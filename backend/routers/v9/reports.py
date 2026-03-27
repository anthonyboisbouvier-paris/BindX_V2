"""
BindX V9 — Report builder router.

Generates professional PDF reports for a phase's molecules using ReportLab.
"""

from __future__ import annotations

import io
import logging
from datetime import datetime, timezone
from typing import List, Optional
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field
from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import selectinload

from auth_v9 import require_v9_user
from database_v9 import get_v9_db
from models_v9 import (
    CampaignORM_V9,
    MoleculeORM_V9,
    MoleculePropertyORM_V9,
    PhaseORM_V9,
    ProjectORM_V9,
    RunORM_V9,
)
from routers.v9.deps import get_phase_owned

logger = logging.getLogger(__name__)

router = APIRouter()

AVAILABLE_SECTIONS = [
    "summary",
    "top_molecules",
    "property_distributions",
    "admet_profiles",
    "safety_overview",
    "pareto_analysis",
    "retrosynthesis",
    "clustering",
    "scaffold",
]


class ReportRequest(BaseModel):
    sections: List[str] = Field(default=AVAILABLE_SECTIONS)
    top_n: int = Field(default=10, ge=1, le=100)
    format: str = Field(default="pdf")
    molecule_ids: Optional[List[str]] = Field(default=None, description="Filter to specific molecule IDs (max 200)")


# ---------------------------------------------------------------------------
# Helpers — flatten properties (same logic as frontend flattenMoleculeProperties)
# ---------------------------------------------------------------------------

_PROP_ALIASES = {
    "hbd": "HBD", "hba": "HBA", "qed": "QED", "tpsa": "TPSA",
    "bbb_permeability": "BBB", "herg_inhibition": "hERG",
    "color_code": "safety_color_code",
    "affinity": "docking_score", "vina_score": "docking_score",
}

_SKIP_KEYS = {
    "flags", "note", "status", "smiles", "confidence_modifier",
    "nearest_tanimoto", "docking_status", "tree", "children",
    "reaction", "reactants", "reactant_names", "conditions",
    "results", "warnings", "pose_molblock", "preparation",
    "breakdown", "weights_used", "scaffold_positions", "scaffold_svg",
    "interactions_detail", "interactions_method", "mapped_functional",
}


def _deep_flatten(obj: dict, out: dict | None = None) -> dict:
    if out is None:
        out = {}
    for k, v in obj.items():
        if k in _SKIP_KEYS:
            continue
        if isinstance(v, list):
            continue
        if isinstance(v, dict):
            _deep_flatten(v, out)
        else:
            alias = _PROP_ALIASES.get(k, k)
            if isinstance(v, str) and v != "":
                try:
                    out[alias] = float(v)
                except ValueError:
                    out[alias] = v
            else:
                out[alias] = v
    return out


def flatten_molecule(mol: MoleculeORM_V9) -> dict:
    """Flatten a molecule ORM object into a flat dict with all properties."""
    flat = {
        "id": str(mol.id),
        "name": mol.name or "",
        "smiles": mol.smiles or "",
        "bookmarked": mol.bookmarked,
    }
    props = {}
    for prop in (mol.properties or []):
        if prop.property_value and isinstance(prop.property_value, dict):
            props[prop.property_name] = prop.property_value
    flat.update(_deep_flatten(props))
    return flat


# ---------------------------------------------------------------------------
# PDF generation with ReportLab
# ---------------------------------------------------------------------------

def _build_pdf(
    phase_name: str,
    project_name: str,
    campaign_name: str,
    molecules: list[dict],
    sections: list[str],
    top_n: int,
    run_count: int,
) -> bytes:
    """Build PDF bytes using ReportLab."""
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import mm, cm
    from reportlab.platypus import (
        SimpleDocTemplate,
        Paragraph,
        Spacer,
        Table,
        TableStyle,
        PageBreak,
    )
    from reportlab.graphics.shapes import Drawing, String
    from reportlab.graphics.charts.barcharts import VerticalBarChart

    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4, topMargin=2 * cm, bottomMargin=2 * cm)
    styles = getSampleStyleSheet()

    # Custom styles
    title_style = ParagraphStyle("BXTitle", parent=styles["Title"], fontSize=24, spaceAfter=6, textColor=colors.HexColor("#1e293b"))
    h1_style = ParagraphStyle("BXH1", parent=styles["Heading1"], fontSize=16, spaceBefore=18, spaceAfter=8, textColor=colors.HexColor("#1e3a5f"))
    h2_style = ParagraphStyle("BXH2", parent=styles["Heading2"], fontSize=13, spaceBefore=12, spaceAfter=6, textColor=colors.HexColor("#334155"))
    body_style = ParagraphStyle("BXBody", parent=styles["Normal"], fontSize=10, leading=14)
    small_style = ParagraphStyle("BXSmall", parent=styles["Normal"], fontSize=8, textColor=colors.HexColor("#64748b"))

    story = []

    # Header color
    header_bg = colors.HexColor("#1e3a5f")
    header_fg = colors.white
    alt_row = colors.HexColor("#f8fafc")

    # Sort by composite_score desc
    sorted_mols = sorted(molecules, key=lambda m: m.get("composite_score") or m.get("weighted_score") or 0, reverse=True)
    top_mols = sorted_mols[:top_n]

    total = len(molecules)
    bookmarked = sum(1 for m in molecules if m.get("bookmarked"))
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    # === COVER PAGE ===
    story.append(Spacer(1, 60))
    story.append(Paragraph("BindX", title_style))
    story.append(Spacer(1, 10))
    story.append(Paragraph(f"Phase Report — {phase_name}", h1_style))
    story.append(Spacer(1, 8))
    story.append(Paragraph(f"Project: {project_name}", body_style))
    story.append(Paragraph(f"Campaign: {campaign_name}", body_style))
    story.append(Paragraph(f"Generated: {now}", body_style))
    story.append(Spacer(1, 20))
    story.append(Paragraph(f"{total} molecules | {bookmarked} bookmarked | {run_count} runs", body_style))
    story.append(PageBreak())

    def _table_style():
        return TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), header_bg),
            ("TEXTCOLOR", (0, 0), (-1, 0), header_fg),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, 0), 8),
            ("FONTSIZE", (0, 1), (-1, -1), 7),
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#e2e8f0")),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, alt_row]),
            ("TOPPADDING", (0, 0), (-1, -1), 4),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ])

    def _fmt(val, decimals=2):
        if val is None:
            return "—"
        if isinstance(val, bool):
            return "Yes" if val else "No"
        if isinstance(val, float):
            return f"{val:.{decimals}f}"
        return str(val)

    # === SUMMARY ===
    if "summary" in sections:
        story.append(Paragraph("1. Summary", h1_style))
        summary_data = [
            ["Metric", "Value"],
            ["Total Molecules", str(total)],
            ["Bookmarked", str(bookmarked)],
            ["Runs Completed", str(run_count)],
            ["Report Date", now],
        ]
        # Add average composite score
        scores = [m.get("composite_score") or m.get("weighted_score") for m in molecules if (m.get("composite_score") or m.get("weighted_score")) is not None]
        if scores:
            avg_score = sum(scores) / len(scores)
            summary_data.append(["Avg Composite Score", f"{avg_score:.3f}"])
        # Lipinski pass rate
        lip_pass = [m for m in molecules if m.get("lipinski_pass") is True]
        if lip_pass or any(m.get("lipinski_pass") is not None for m in molecules):
            lip_total = sum(1 for m in molecules if m.get("lipinski_pass") is not None)
            rate = len(lip_pass) / lip_total * 100 if lip_total else 0
            summary_data.append(["Lipinski Pass Rate", f"{rate:.1f}%"])

        t = Table(summary_data, colWidths=[150, 200])
        t.setStyle(_table_style())
        story.append(t)
        story.append(Spacer(1, 12))

    # === TOP MOLECULES ===
    if "top_molecules" in sections and top_mols:
        story.append(Paragraph(f"2. Top {min(top_n, len(top_mols))} Molecules", h1_style))
        headers = ["#", "Name", "SMILES", "Dock.", "QED", "logP", "Safety", "Score"]
        data = [headers]
        for i, m in enumerate(top_mols, 1):
            smiles = m.get("smiles", "")
            if len(smiles) > 30:
                smiles = smiles[:30] + "..."
            data.append([
                str(i),
                m.get("name", "—")[:20],
                smiles,
                _fmt(m.get("docking_score")),
                _fmt(m.get("QED")),
                _fmt(m.get("logP")),
                m.get("safety_color_code", "—"),
                _fmt(m.get("composite_score") or m.get("weighted_score")),
            ])
        col_widths = [25, 70, 120, 45, 40, 40, 45, 45]
        t = Table(data, colWidths=col_widths)
        t.setStyle(_table_style())
        story.append(t)
        story.append(Spacer(1, 12))

    # === PROPERTY DISTRIBUTIONS ===
    if "property_distributions" in sections:
        story.append(Paragraph("3. Property Distributions", h1_style))
        dist_keys = ["docking_score", "composite_score", "MW", "logP", "QED", "TPSA"]
        for key in dist_keys:
            values = [m.get(key) for m in molecules if m.get(key) is not None and isinstance(m.get(key), (int, float))]
            if len(values) < 3:
                continue
            story.append(Paragraph(f"{key}", h2_style))

            # Compute histogram bins
            min_v, max_v = min(values), max(values)
            if min_v == max_v:
                story.append(Paragraph(f"All values = {min_v:.2f}", body_style))
                continue

            n_bins = min(20, len(values) // 2 + 1)
            bin_width = (max_v - min_v) / n_bins
            bins = [0] * n_bins
            for v in values:
                idx = min(int((v - min_v) / bin_width), n_bins - 1)
                bins[idx] += 1

            # Draw bar chart
            drawing = Drawing(400, 120)
            chart = VerticalBarChart()
            chart.x = 40
            chart.y = 15
            chart.width = 340
            chart.height = 90
            chart.data = [bins]
            chart.categoryAxis.categoryNames = [f"{min_v + i * bin_width:.1f}" for i in range(n_bins)]
            chart.categoryAxis.labels.fontSize = 5
            chart.categoryAxis.labels.angle = 45
            chart.valueAxis.labels.fontSize = 6
            chart.bars[0].fillColor = colors.HexColor("#3b82f6")
            drawing.add(chart)
            story.append(drawing)
            story.append(Spacer(1, 6))

            # Stats
            import statistics
            mean_v = statistics.mean(values)
            std_v = statistics.stdev(values) if len(values) > 1 else 0
            story.append(Paragraph(
                f"n={len(values)}, mean={mean_v:.2f}, std={std_v:.2f}, min={min_v:.2f}, max={max_v:.2f}",
                small_style
            ))
            story.append(Spacer(1, 8))

    # === ADMET PROFILES ===
    if "admet_profiles" in sections and top_mols:
        story.append(Paragraph(f"4. ADMET Profiles (Top {min(top_n, len(top_mols))})", h1_style))
        admet_keys = ["solubility", "BBB", "oral_bioavailability", "plasma_protein_binding", "half_life", "hERG"]
        headers = ["Name"] + [k.replace("_", " ").title() for k in admet_keys]
        data = [headers]
        for m in top_mols:
            row = [m.get("name", "—")[:20]]
            for k in admet_keys:
                row.append(_fmt(m.get(k)))
            data.append(row)
        col_widths = [80] + [60] * len(admet_keys)
        t = Table(data, colWidths=col_widths)
        t.setStyle(_table_style())
        story.append(t)
        story.append(Spacer(1, 12))

    # === SAFETY OVERVIEW ===
    if "safety_overview" in sections:
        story.append(Paragraph("5. Safety Overview", h1_style))
        safety_counts = {}
        for m in molecules:
            sc = m.get("safety_color_code", "unknown")
            safety_counts[sc] = safety_counts.get(sc, 0) + 1

        data = [["Safety Color", "Count", "Percentage"]]
        for color, count in sorted(safety_counts.items(), key=lambda x: -x[1]):
            pct = count / total * 100 if total else 0
            data.append([color, str(count), f"{pct:.1f}%"])
        t = Table(data, colWidths=[120, 80, 80])
        t.setStyle(_table_style())
        story.append(t)

        # PAINS / Brenk alerts
        pains_count = sum(1 for m in molecules if m.get("pains_alert") is True)
        brenk_count = sum(1 for m in molecules if m.get("brenk_alert") is True)
        if pains_count or brenk_count:
            story.append(Spacer(1, 8))
            story.append(Paragraph(f"PAINS alerts: {pains_count} | Brenk alerts: {brenk_count}", body_style))
        story.append(Spacer(1, 12))

    # === PARETO ANALYSIS ===
    if "pareto_analysis" in sections:
        story.append(Paragraph("6. Pareto Analysis", h1_style))
        story.append(Paragraph(
            "Pareto analysis identifies molecules that represent optimal trade-offs "
            "between competing objectives (e.g., docking score vs. ADME properties). "
            "These molecules are not dominated by any other molecule on all criteria.",
            body_style
        ))
        # Count molecules with both docking and composite scores
        pareto_candidates = [m for m in molecules if m.get("docking_score") is not None and (m.get("composite_score") or m.get("weighted_score")) is not None]
        story.append(Paragraph(f"{len(pareto_candidates)} molecules have both docking and composite scores for Pareto analysis.", body_style))
        story.append(Spacer(1, 12))

    # === RETROSYNTHESIS ===
    if "retrosynthesis" in sections and top_mols:
        synth_mols = [m for m in top_mols if m.get("n_synth_steps") is not None]
        if synth_mols:
            story.append(Paragraph(f"7. Retrosynthesis (Top {len(synth_mols)})", h1_style))
            headers = ["Name", "Steps", "Confidence", "Cost", "Reagents"]
            data = [headers]
            for m in synth_mols:
                data.append([
                    m.get("name", "—")[:20],
                    _fmt(m.get("n_synth_steps"), 0),
                    _fmt(m.get("synth_confidence")),
                    m.get("synth_cost_estimate", "—"),
                    "Yes" if m.get("reagents_available") else "No",
                ])
            t = Table(data, colWidths=[90, 50, 70, 70, 60])
            t.setStyle(_table_style())
            story.append(t)
            story.append(Spacer(1, 12))

    # === CLUSTERING ===
    if "clustering" in sections or "scaffold" in sections:
        clustered = [m for m in molecules if m.get("cluster_id") is not None]
        if clustered:
            story.append(Paragraph("8. Clustering Summary", h1_style))
            clusters = {}
            for m in clustered:
                cid = m.get("cluster_id")
                if cid not in clusters:
                    clusters[cid] = {"count": 0, "representatives": 0}
                clusters[cid]["count"] += 1
                if m.get("is_representative"):
                    clusters[cid]["representatives"] += 1

            headers = ["Cluster", "Members", "Representatives"]
            data = [headers]
            for cid in sorted(clusters.keys()):
                info = clusters[cid]
                data.append([str(cid), str(info["count"]), str(info["representatives"])])
            t = Table(data, colWidths=[80, 80, 100])
            t.setStyle(_table_style())
            story.append(t)
            story.append(Spacer(1, 12))

    # Build PDF
    doc.build(story)
    return buf.getvalue()


# ---------------------------------------------------------------------------
# Endpoint
# ---------------------------------------------------------------------------

@router.post("/phases/{phase_id}/report")
async def generate_report(
    phase_id: UUID,
    body: ReportRequest,
    user_id: str = Depends(require_v9_user),
    db: AsyncSession = Depends(get_v9_db),
):
    """Generate a PDF report for a phase."""

    phase = await get_phase_owned(phase_id, user_id, db)

    # Validate sections
    invalid = [s for s in body.sections if s not in AVAILABLE_SECTIONS]
    if invalid:
        raise HTTPException(status_code=422, detail=f"Invalid sections: {invalid}")

    # Validate molecule_ids limit
    if body.molecule_ids and len(body.molecule_ids) > 200:
        raise HTTPException(status_code=422, detail="Maximum 200 molecules per report")

    # Load molecules with properties
    stmt = (
        select(MoleculeORM_V9)
        .where(MoleculeORM_V9.phase_id == phase.id)
        .options(selectinload(MoleculeORM_V9.properties))
    )
    if body.molecule_ids:
        from uuid import UUID as _UUID
        try:
            uuids = [_UUID(mid) for mid in body.molecule_ids]
        except (ValueError, AttributeError):
            raise HTTPException(status_code=422, detail="Invalid molecule IDs")
        stmt = stmt.where(MoleculeORM_V9.id.in_(uuids))
    result = await db.execute(stmt)
    molecules_orm = result.scalars().all()

    if not molecules_orm:
        raise HTTPException(status_code=404, detail="No molecules found")

    molecules = [flatten_molecule(m) for m in molecules_orm]

    # Get run count
    run_stmt = select(RunORM_V9).where(RunORM_V9.phase_id == phase.id, RunORM_V9.status == "completed")
    run_result = await db.execute(run_stmt)
    run_count = len(run_result.scalars().all())

    # Get phase/campaign/project names
    phase_type_labels = {
        "hit_discovery": "Phase A — Hit Discovery",
        "hit_to_lead": "Phase B — Hit-to-Lead",
        "lead_optimization": "Phase C — Lead Optimization",
    }
    phase_name = phase_type_labels.get(phase.type, phase.type)

    # Load campaign and project names
    campaign = phase.campaign
    project = campaign.project if campaign else None
    campaign_name = campaign.name if campaign else "—"
    project_name = project.name if project else "—"

    # Generate PDF
    try:
        pdf_bytes = _build_pdf(
            phase_name=phase_name,
            project_name=project_name,
            campaign_name=campaign_name,
            molecules=molecules,
            sections=body.sections,
            top_n=body.top_n,
            run_count=run_count,
        )
    except Exception as exc:
        logger.error("PDF generation failed: %s", exc, exc_info=True)
        raise HTTPException(status_code=500, detail=f"PDF generation failed: {exc}")

    filename = f"BindX_Report_{phase.type}_{datetime.now(timezone.utc).strftime('%Y%m%d')}.pdf"

    return StreamingResponse(
        io.BytesIO(pdf_bytes),
        media_type="application/pdf",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )
