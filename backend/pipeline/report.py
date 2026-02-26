"""
DockIt pipeline — PDF report and ZIP archive generation.

Uses ReportLab for PDF creation and the standard-library ``zipfile``
module for archive bundling.
"""

from __future__ import annotations

import csv
import io
import logging
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# PDF report
# ---------------------------------------------------------------------------

def generate_pdf_report(
    job_data: dict,
    results: list[dict],
    output_path: Path,
    generated_results: list[dict] | None = None,
) -> Path:
    """Generate a PDF report summarising the docking run.

    Parameters
    ----------
    job_data : dict
        Job metadata (``job_id``, ``uniprot_id``, ``protein_name``, ``mode``).
    results : list[dict]
        Scored docking results (top N).
    output_path : Path
        Where to write the PDF file.
    generated_results : list[dict], optional
        V2: AI-generated molecules results.

    Returns
    -------
    Path
        Absolute path to the generated PDF.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        return _generate_pdf_reportlab(job_data, results, output_path, generated_results)
    except ImportError:
        logger.warning("ReportLab not available; generating plain-text report")
        return _generate_txt_report(job_data, results, output_path, generated_results)
    except Exception as exc:
        logger.error("PDF generation failed: %s", exc)
        return _generate_txt_report(job_data, results, output_path, generated_results)


def _generate_pdf_reportlab(
    job_data: dict,
    results: list[dict],
    output_path: Path,
    generated_results: list[dict] | None = None,
) -> Path:
    """Create PDF using ReportLab.

    Sections
    --------
    1. Header / job metadata
    2. Top 10 results table
    3. SMILES of top compounds
    4. AI-generated molecules table (if any)
    5. **Detailed Candidate Profiles** for the top 5 (ADMET, synthesis,
       off-target, confidence, Pareto)
    6. **Pipeline Summary** (structure source, pocket, cutoffs, clusters)
    7. Legacy ADMET / retrosynthesis summaries
    8. Footer
    """
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.platypus import (
        HRFlowable,
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    doc = SimpleDocTemplate(
        str(output_path),
        pagesize=A4,
        leftMargin=20 * mm,
        rightMargin=20 * mm,
        topMargin=20 * mm,
        bottomMargin=20 * mm,
    )

    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        "DockItTitle",
        parent=styles["Title"],
        fontSize=22,
        textColor=colors.HexColor("#1e3a5f"),
        spaceAfter=12,
    )
    heading_style = ParagraphStyle(
        "DockItHeading",
        parent=styles["Heading2"],
        textColor=colors.HexColor("#1e3a5f"),
        spaceAfter=8,
    )
    subheading_style = ParagraphStyle(
        "DockItSubHeading",
        parent=styles["Heading3"],
        textColor=colors.HexColor("#1e3a5f"),
        fontSize=11,
        spaceAfter=4,
        spaceBefore=6,
    )
    detail_style = ParagraphStyle(
        "DockItDetail",
        parent=styles["Normal"],
        fontSize=8,
        leading=11,
        leftIndent=10,
        spaceAfter=2,
    )
    detail_bold_style = ParagraphStyle(
        "DockItDetailBold",
        parent=detail_style,
        fontName="Helvetica-Bold",
        fontSize=8,
    )

    elements: list[Any] = []

    # -----------------------------------------------------------------------
    # Title
    # -----------------------------------------------------------------------
    elements.append(Paragraph("BindX - Virtual Screening Report", title_style))
    elements.append(Spacer(1, 6 * mm))

    # -----------------------------------------------------------------------
    # Job info
    # -----------------------------------------------------------------------
    job_id = job_data.get("job_id", "N/A")
    uniprot = job_data.get("uniprot_id", "N/A")
    protein = job_data.get("protein_name", uniprot)
    date_str = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")

    mode = job_data.get("mode", "basic")
    n_generated = len(generated_results) if generated_results else 0

    info_text = (
        f"<b>Job ID:</b> {job_id}<br/>"
        f"<b>Target protein:</b> {protein} ({uniprot})<br/>"
        f"<b>Date:</b> {date_str}<br/>"
        f"<b>Mode:</b> {mode.upper()}<br/>"
        f"<b>Known ligands screened:</b> {len(results)}"
    )
    if n_generated:
        info_text += f"<br/><b>AI-generated molecules:</b> {n_generated}"
    elements.append(Paragraph(info_text, styles["Normal"]))
    elements.append(Spacer(1, 8 * mm))

    # -----------------------------------------------------------------------
    # Results table (top 10)
    # -----------------------------------------------------------------------
    elements.append(Paragraph("Top Docking Results", heading_style))

    top_n = results[:10]
    table_data = [["Rank", "Name", "Affinity\n(kcal/mol)", "QED", "LogP", "MW", "Score", "ADMET"]]
    for i, r in enumerate(top_n, 1):
        admet = r.get("admet", {})
        admet_label = "-"
        if isinstance(admet, dict) and admet.get("color_code"):
            admet_label = admet["color_code"].upper()
        table_data.append([
            str(i),
            str(r.get("name", "?"))[:40],
            f"{r.get('affinity', 0):.1f}",
            f"{r.get('qed', 0):.2f}" if r.get("qed") is not None else "-",
            f"{r.get('logP') or r.get('logp', 0):.1f}" if r.get("logP") or r.get("logp") is not None else "-",
            f"{r.get('MW') or r.get('mw', 0):.0f}" if r.get("MW") or r.get("mw") is not None else "-",
            f"{r.get('composite_score', 0):.3f}",
            admet_label,
        ])

    table = Table(table_data, repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1e3a5f")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, 0), 9),
        ("FONTSIZE", (0, 1), (-1, -1), 8),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.whitesmoke, colors.white]),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    elements.append(table)
    elements.append(Spacer(1, 10 * mm))

    # -----------------------------------------------------------------------
    # SMILES list (name truncated to 40 chars)
    # -----------------------------------------------------------------------
    elements.append(Paragraph("SMILES of Top Compounds", heading_style))
    for i, r in enumerate(top_n, 1):
        smi = r.get("smiles", "N/A")
        name = str(r.get("name", "?"))[:40]
        elements.append(Paragraph(
            f"<b>{i}. {name}:</b> <font size=7>{smi}</font>",
            styles["Normal"],
        ))
    elements.append(Spacer(1, 6 * mm))

    # -----------------------------------------------------------------------
    # AI-generated molecules table (V2)
    # -----------------------------------------------------------------------
    if generated_results:
        elements.append(Paragraph("AI-Generated Molecules", heading_style))
        elements.append(Paragraph(
            f"REINVENT4 generated <b>{len(generated_results)}</b> novel molecules optimized for this target.",
            styles["Normal"],
        ))
        elements.append(Spacer(1, 4 * mm))

        gen_data = [["Rank", "Name", "Affinity\n(kcal/mol)", "QED", "Novelty", "ADMET", "Score"]]
        for i, r in enumerate(generated_results[:10], 1):
            admet = r.get("admet", {})
            admet_score = admet.get("composite_score", "-") if isinstance(admet, dict) else "-"
            gen_data.append([
                str(i),
                str(r.get("name", "?"))[:20],
                f"{r.get('affinity', 0):.1f}",
                f"{r.get('qed', 0):.2f}" if r.get("qed") is not None else "-",
                f"{r.get('novelty_score', 0):.2f}" if r.get("novelty_score") is not None else "-",
                f"{admet_score:.2f}" if isinstance(admet_score, (int, float)) else "-",
                f"{r.get('composite_score', 0):.3f}",
            ])

        gen_table = Table(gen_data, repeatRows=1)
        gen_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#7c3aed")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, 0), 9),
            ("FONTSIZE", (0, 1), (-1, -1), 8),
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.Color(0.96, 0.95, 1.0), colors.white]),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("TOPPADDING", (0, 0), (-1, -1), 4),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ]))
        elements.append(gen_table)
        elements.append(Spacer(1, 8 * mm))

    # ===================================================================
    # NEW SECTION: Detailed Candidate Profiles (top 5)
    # ===================================================================
    try:
        top5 = results[:5]
        if top5:
            elements.append(Paragraph("Detailed Candidate Profiles", heading_style))
            elements.append(Paragraph(
                "In-depth analysis of the top 5 candidates, including docking scores, "
                "ADMET properties, synthesis feasibility, off-target safety, and confidence.",
                styles["Normal"],
            ))
            elements.append(Spacer(1, 4 * mm))

            for rank, r in enumerate(top5, 1):
                try:
                    cand_name = str(r.get("name", "Unknown"))[:40]
                    elements.append(Paragraph(
                        f"#{rank} — {cand_name}",
                        subheading_style,
                    ))

                    # -- Core identifiers --
                    smi = r.get("smiles", "N/A")
                    elements.append(Paragraph(
                        f"<b>SMILES:</b> <font size=7>{smi}</font>",
                        detail_style,
                    ))

                    # -- Docking scores --
                    score_parts: list[str] = []
                    if r.get("affinity") is not None:
                        score_parts.append(f"Affinity: {r['affinity']:.1f} kcal/mol")
                    if r.get("cnn_score") is not None:
                        score_parts.append(f"CNN Score: {r['cnn_score']:.3f}")
                    if r.get("cnn_affinity") is not None:
                        score_parts.append(f"CNN Affinity: {r['cnn_affinity']:.2f} pK")
                    if r.get("vina_score") is not None:
                        score_parts.append(f"Vina Score: {r['vina_score']:.1f}")
                    if r.get("composite_score") is not None:
                        score_parts.append(f"Composite: {r['composite_score']:.3f}")
                    if r.get("score_100") is not None and r["score_100"] > 0:
                        score_parts.append(f"Score/100: {r['score_100']}")
                    if score_parts:
                        elements.append(Paragraph(
                            f"<b>Scores:</b> {' | '.join(score_parts)}",
                            detail_style,
                        ))

                    # -- Drug-likeness --
                    dl_parts: list[str] = []
                    if r.get("qed") is not None:
                        dl_parts.append(f"QED: {r['qed']:.2f}")
                    logp_val = r.get("logP") or r.get("logp")
                    if logp_val is not None:
                        dl_parts.append(f"LogP: {logp_val:.1f}")
                    mw_val = r.get("MW") or r.get("mw")
                    if mw_val is not None:
                        dl_parts.append(f"MW: {mw_val:.0f}")
                    if r.get("tpsa") is not None:
                        dl_parts.append(f"TPSA: {r['tpsa']:.1f}")
                    if r.get("hbd") is not None:
                        dl_parts.append(f"HBD: {r['hbd']}")
                    if r.get("hba") is not None:
                        dl_parts.append(f"HBA: {r['hba']}")
                    if r.get("rotatable_bonds") is not None:
                        dl_parts.append(f"RotBonds: {r['rotatable_bonds']}")
                    if r.get("sa_score") is not None:
                        dl_parts.append(f"SA: {r['sa_score']:.1f}")
                    if dl_parts:
                        elements.append(Paragraph(
                            f"<b>Drug-likeness:</b> {' | '.join(dl_parts)}",
                            detail_style,
                        ))

                    # -- Pareto ranking --
                    pareto_parts: list[str] = []
                    if r.get("pareto_rank") is not None:
                        pareto_parts.append(f"Rank: {r['pareto_rank']}")
                    if r.get("pareto_front") is not None:
                        pareto_parts.append(
                            "On Pareto front" if r["pareto_front"] else "Not on Pareto front"
                        )
                    if r.get("pareto_objectives") and isinstance(r["pareto_objectives"], dict):
                        obj_strs = [
                            f"{k}: {v:.3f}" if isinstance(v, float) else f"{k}: {v}"
                            for k, v in r["pareto_objectives"].items()
                        ]
                        pareto_parts.append(f"Objectives({', '.join(obj_strs)})")
                    if pareto_parts:
                        elements.append(Paragraph(
                            f"<b>Pareto:</b> {' | '.join(pareto_parts)}",
                            detail_style,
                        ))

                    # -- Cluster info --
                    if r.get("cluster_id") is not None:
                        cl_parts = [f"Cluster #{r['cluster_id']}"]
                        if r.get("cluster_size") is not None:
                            cl_parts.append(f"size {r['cluster_size']}")
                        if r.get("is_representative"):
                            cl_parts.append("(representative)")
                        elements.append(Paragraph(
                            f"<b>Cluster:</b> {', '.join(cl_parts)}",
                            detail_style,
                        ))

                    # -- Elimination status --
                    if r.get("eliminated"):
                        reason = r.get("elimination_reason", "unspecified")
                        elements.append(Paragraph(
                            f'<b>Status:</b> <font color="red">ELIMINATED</font> — {reason}',
                            detail_style,
                        ))

                    # -- ADMET details --
                    admet = r.get("admet")
                    if isinstance(admet, dict) and admet:
                        admet_color = admet.get("color_code", "yellow").upper()
                        admet_comp = admet.get("composite_score")
                        admet_header = f"<b>ADMET:</b> {admet_color}"
                        if admet_comp is not None:
                            admet_header += f" (score: {admet_comp:.2f})"
                        flags = admet.get("flags", [])
                        if flags:
                            admet_header += f" — Flags: {', '.join(str(f) for f in flags)}"
                        else:
                            admet_header += " — No alerts"
                        elements.append(Paragraph(admet_header, detail_style))

                        # Individual ADMET fields (non-None values)
                        admet_fields = [
                            ("oral_bioavailability", "Oral Bioavailability"),
                            ("intestinal_permeability", "Intestinal Permeability"),
                            ("solubility", "Solubility"),
                            ("plasma_protein_binding", "Plasma Protein Binding"),
                            ("bbb_permeability", "BBB Permeability"),
                            ("cyp1a2_inhibitor", "CYP1A2 Inhibitor"),
                            ("cyp2c9_inhibitor", "CYP2C9 Inhibitor"),
                            ("cyp2c19_inhibitor", "CYP2C19 Inhibitor"),
                            ("cyp2d6_inhibitor", "CYP2D6 Inhibitor"),
                            ("cyp3a4_inhibitor", "CYP3A4 Inhibitor"),
                            ("clearance", "Clearance"),
                            ("half_life", "Half-life"),
                            ("herg_inhibition", "hERG Inhibition"),
                            ("ames_mutagenicity", "Ames Mutagenicity"),
                            ("hepatotoxicity", "Hepatotoxicity"),
                            ("skin_sensitization", "Skin Sensitization"),
                            ("carcinogenicity", "Carcinogenicity"),
                        ]
                        field_strs: list[str] = []
                        for key, label in admet_fields:
                            val = admet.get(key)
                            if val is not None:
                                field_strs.append(f"{label}: {val:.2f}")
                        if field_strs:
                            # Display as a compact multi-column line
                            elements.append(Paragraph(
                                "<b>ADMET Details:</b> " + " | ".join(field_strs),
                                detail_style,
                            ))

                    # -- Synthesis route --
                    synth = r.get("synthesis_route")
                    if isinstance(synth, dict) and synth:
                        n_steps = synth.get("n_steps", 0)
                        conf = synth.get("confidence", 0)
                        avail = synth.get("all_reagents_available", False)
                        est_cost = synth.get("estimated_cost", "unknown")
                        avail_txt = "all reagents available" if avail else "some reagents unavailable"
                        synth_line = (
                            f"<b>Synthesis:</b> {n_steps} steps, "
                            f"confidence {conf:.0%}, "
                            f"cost: {est_cost}, {avail_txt}"
                        )
                        elements.append(Paragraph(synth_line, detail_style))

                        # V6.3 cost_estimate
                        cost_est = synth.get("cost_estimate")
                        if isinstance(cost_est, dict) and cost_est:
                            cost_parts: list[str] = []
                            if cost_est.get("total_cost_usd") is not None:
                                cost_parts.append(f"Total: ${cost_est['total_cost_usd']:.0f}")
                            if cost_est.get("reagent_cost") is not None:
                                cost_parts.append(f"Reagents: ${cost_est['reagent_cost']:.0f}")
                            if cost_est.get("labor_cost") is not None:
                                cost_parts.append(f"Labor: ${cost_est['labor_cost']:.0f}")
                            if cost_est.get("currency"):
                                cost_parts.append(f"Currency: {cost_est['currency']}")
                            if cost_parts:
                                elements.append(Paragraph(
                                    f"<b>Cost Estimate:</b> {' | '.join(cost_parts)}",
                                    detail_style,
                                ))

                    # -- Off-target results --
                    off_target = r.get("off_target")
                    if isinstance(off_target, dict) and off_target:
                        sel_score = off_target.get("selectivity_score")
                        n_safe = off_target.get("n_safe")
                        n_total = off_target.get("n_total")
                        warnings = off_target.get("warnings", [])
                        ot_parts: list[str] = []
                        if sel_score is not None:
                            ot_parts.append(f"Selectivity: {sel_score:.2f}")
                        if n_safe is not None and n_total is not None:
                            ot_parts.append(f"Safe: {n_safe}/{n_total}")
                        ot_line = "<b>Off-target:</b> "
                        if ot_parts:
                            ot_line += " | ".join(ot_parts)
                        if warnings:
                            ot_line += f' — <font color="red">Warnings: {", ".join(str(w) for w in warnings)}</font>'
                        elif ot_parts:
                            ot_line += ' — <font color="green">No warnings</font>'
                        if ot_parts or warnings:
                            elements.append(Paragraph(ot_line, detail_style))

                    # -- Combined off-target (V6.3) --
                    comb_ot = r.get("combined_off_target")
                    if isinstance(comb_ot, dict) and comb_ot:
                        cot_parts: list[str] = []
                        if comb_ot.get("selectivity_score") is not None:
                            cot_parts.append(f"Combined Selectivity: {comb_ot['selectivity_score']:.2f}")
                        if cot_parts:
                            elements.append(Paragraph(
                                f"<b>Combined Off-target:</b> {' | '.join(cot_parts)}",
                                detail_style,
                            ))

                    # -- hERG specialized (V6.3) --
                    herg = r.get("herg_specialized")
                    if isinstance(herg, dict) and herg:
                        herg_parts: list[str] = []
                        if herg.get("predicted_ic50") is not None:
                            herg_parts.append(f"Predicted IC50: {herg['predicted_ic50']:.2f}")
                        if herg.get("risk_level"):
                            herg_parts.append(f"Risk: {herg['risk_level']}")
                        if herg_parts:
                            elements.append(Paragraph(
                                f"<b>hERG:</b> {' | '.join(herg_parts)}",
                                detail_style,
                            ))

                    # -- Confidence score --
                    conf_data = r.get("confidence")
                    if isinstance(conf_data, dict) and conf_data:
                        overall = conf_data.get("overall")
                        conf_line = "<b>Confidence:</b> "
                        if overall is not None:
                            conf_line += f"Overall: {overall:.2f}"
                        # List sub-scores
                        sub_parts: list[str] = []
                        for ckey, cval in conf_data.items():
                            if ckey == "overall":
                                continue
                            if isinstance(cval, (int, float)):
                                sub_parts.append(f"{ckey}: {cval:.2f}")
                        if sub_parts:
                            conf_line += " (" + ", ".join(sub_parts) + ")"
                        elements.append(Paragraph(conf_line, detail_style))

                    # Separator between candidates
                    elements.append(Spacer(1, 2 * mm))
                    if rank < len(top5):
                        elements.append(HRFlowable(
                            width="100%", thickness=0.5,
                            color=colors.lightgrey, spaceAfter=2,
                        ))

                except Exception as cand_exc:
                    logger.warning(
                        "PDF: failed rendering candidate #%d (%s): %s",
                        rank, r.get("name", "?"), cand_exc,
                    )
                    elements.append(Paragraph(
                        f"<i>Could not render full profile for candidate #{rank}.</i>",
                        detail_style,
                    ))

            elements.append(Spacer(1, 8 * mm))
    except Exception as sec_exc:
        logger.warning("PDF: Detailed Candidate Profiles section failed: %s", sec_exc)

    # ===================================================================
    # NEW SECTION: Pipeline Summary
    # ===================================================================
    try:
        pipeline_summary = job_data.get("pipeline_summary", {})
        if isinstance(pipeline_summary, dict) and pipeline_summary:
            elements.append(Paragraph("Pipeline Summary", heading_style))

            ps_lines: list[str] = []

            # Structure source
            src = pipeline_summary.get("structure_source")
            if src:
                ps_lines.append(f"<b>Structure source:</b> {src}")
            pdb_id = pipeline_summary.get("pdb_id")
            if pdb_id:
                ps_lines.append(f"<b>PDB ID:</b> {pdb_id}")
            resolution = pipeline_summary.get("resolution")
            if resolution is not None:
                ps_lines.append(f"<b>Resolution:</b> {resolution}")

            # Pocket detection
            pocket_method = pipeline_summary.get("pocket_method") or pipeline_summary.get("pocket_detection")
            if pocket_method:
                ps_lines.append(f"<b>Pocket detection:</b> {pocket_method}")
            n_pockets = pipeline_summary.get("n_pockets")
            if n_pockets is not None:
                ps_lines.append(f"<b>Pockets found:</b> {n_pockets}")

            # Docking engine
            engine = pipeline_summary.get("docking_engine")
            if engine:
                ps_lines.append(f"<b>Docking engine:</b> {engine}")

            # Molecules screened
            n_screened = pipeline_summary.get("n_molecules_screened") or pipeline_summary.get("total_screened")
            if n_screened is not None:
                ps_lines.append(f"<b>Molecules screened:</b> {n_screened}")

            # Hard cutoffs
            n_passed = pipeline_summary.get("n_passed") or pipeline_summary.get("passed_cutoffs")
            n_eliminated = pipeline_summary.get("n_eliminated") or pipeline_summary.get("eliminated_cutoffs")
            if n_passed is not None or n_eliminated is not None:
                cutoff_str = "<b>Hard cutoffs:</b>"
                if n_passed is not None:
                    cutoff_str += f" passed={n_passed}"
                if n_eliminated is not None:
                    cutoff_str += f" eliminated={n_eliminated}"
                ps_lines.append(cutoff_str)

            # Clusters
            n_clusters = pipeline_summary.get("n_clusters")
            if n_clusters is not None:
                ps_lines.append(f"<b>Clusters:</b> {n_clusters}")

            # Duration
            duration = pipeline_summary.get("duration_seconds") or pipeline_summary.get("total_time")
            if duration is not None:
                if isinstance(duration, (int, float)):
                    minutes = int(duration) // 60
                    seconds = int(duration) % 60
                    ps_lines.append(f"<b>Total time:</b> {minutes}m {seconds}s")
                else:
                    ps_lines.append(f"<b>Total time:</b> {duration}")

            # Render all collected lines
            for line in ps_lines:
                elements.append(Paragraph(line, detail_style))

            # Also dump any extra keys not already handled as a fallback
            handled_keys = {
                "structure_source", "pdb_id", "resolution",
                "pocket_method", "pocket_detection", "n_pockets",
                "docking_engine",
                "n_molecules_screened", "total_screened",
                "n_passed", "passed_cutoffs",
                "n_eliminated", "eliminated_cutoffs",
                "n_clusters",
                "duration_seconds", "total_time",
            }
            extra_keys = set(pipeline_summary.keys()) - handled_keys
            for ek in sorted(extra_keys):
                val = pipeline_summary[ek]
                if val is not None and val != "" and val != []:
                    elements.append(Paragraph(
                        f"<b>{ek.replace('_', ' ').title()}:</b> {val}",
                        detail_style,
                    ))

            elements.append(Spacer(1, 8 * mm))
    except Exception as ps_exc:
        logger.warning("PDF: Pipeline Summary section failed: %s", ps_exc)

    # -----------------------------------------------------------------------
    # ADMET summary for top molecules (legacy section, kept for compat)
    # -----------------------------------------------------------------------
    try:
        all_with_admet = [r for r in (results[:5] + (generated_results or [])[:5]) if r.get("admet")]
        if all_with_admet:
            elements.append(Paragraph("ADMET Profiles - Top Candidates", heading_style))
            for r in all_with_admet[:5]:
                admet = r.get("admet", {})
                name = r.get("name", "?")
                color = admet.get("color_code", "yellow") if isinstance(admet, dict) else "yellow"
                flags = admet.get("flags", []) if isinstance(admet, dict) else []
                score = admet.get("composite_score", 0) if isinstance(admet, dict) else 0
                flag_text = ", ".join(str(f) for f in flags) if flags else "No alerts"
                elements.append(Paragraph(
                    f"<b>{name}</b> — ADMET Score: {score:.2f} ({color.upper()}) — {flag_text}",
                    styles["Normal"],
                ))
            elements.append(Spacer(1, 8 * mm))
    except Exception as admet_exc:
        logger.warning("PDF: ADMET summary section failed: %s", admet_exc)

    # -----------------------------------------------------------------------
    # Retrosynthesis summary (legacy section, kept for compat)
    # -----------------------------------------------------------------------
    try:
        all_with_synth = [r for r in (results + (generated_results or [])) if r.get("synthesis_route")]
        if all_with_synth:
            elements.append(Paragraph("Retrosynthesis Routes", heading_style))
            for r in all_with_synth[:5]:
                synth = r.get("synthesis_route", {})
                name = r.get("name", "?")
                n_steps = synth.get("n_steps", 0) if isinstance(synth, dict) else 0
                conf = synth.get("confidence", 0) if isinstance(synth, dict) else 0
                available = synth.get("all_reagents_available", False) if isinstance(synth, dict) else False
                cost = synth.get("estimated_cost", "unknown") if isinstance(synth, dict) else "unknown"
                avail_text = "All reagents available" if available else "Some reagents may need custom synthesis"
                elements.append(Paragraph(
                    f"<b>{name}</b> — {n_steps} steps, confidence {conf:.0%}, cost: {cost}. {avail_text}",
                    styles["Normal"],
                ))
            elements.append(Spacer(1, 8 * mm))
    except Exception as synth_exc:
        logger.warning("PDF: Retrosynthesis summary section failed: %s", synth_exc)

    # -----------------------------------------------------------------------
    # Footer
    # -----------------------------------------------------------------------
    elements.append(Paragraph(
        "<i>Generated by BindX V6 - The Canva of molecular docking</i>",
        styles["Normal"],
    ))

    doc.build(elements)
    logger.info("PDF report generated: %s", output_path)
    return output_path


def _generate_txt_report(
    job_data: dict,
    results: list[dict],
    output_path: Path,
    generated_results: list[dict] | None = None,
) -> Path:
    """Plain-text fallback when ReportLab is unavailable."""
    txt_path = output_path.with_suffix(".txt")
    lines = [
        "=" * 60,
        "DockIt - Virtual Screening Report",
        "=" * 60,
        f"Job ID:   {job_data.get('job_id', 'N/A')}",
        f"Target:   {job_data.get('protein_name', job_data.get('uniprot_id', 'N/A'))}",
        f"UniProt:  {job_data.get('uniprot_id', 'N/A')}",
        f"Date:     {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')}",
        f"Ligands:  {len(results)}",
        "",
        "-" * 60,
        f"{'Rank':<5}{'Name':<25}{'Affinity':<10}{'QED':<8}{'Score':<8}",
        "-" * 60,
    ]
    for i, r in enumerate(results[:10], 1):
        lines.append(
            f"{i:<5}{str(r.get('name','?'))[:24]:<25}"
            f"{r.get('affinity', 0):>8.1f}  "
            f"{r.get('qed', 0):>6.2f}  "
            f"{r.get('composite_score', 0):>6.3f}"
        )
    lines.append("-" * 60)

    if generated_results:
        lines.append("")
        lines.append("=" * 60)
        lines.append("AI-GENERATED MOLECULES")
        lines.append("=" * 60)
        lines.append(f"{'Rank':<5}{'Name':<25}{'Affinity':<10}{'QED':<8}{'Score':<8}")
        lines.append("-" * 60)
        for i, r in enumerate(generated_results[:10], 1):
            lines.append(
                f"{i:<5}{str(r.get('name','?'))[:24]:<25}"
                f"{r.get('affinity', 0):>8.1f}  "
                f"{r.get('qed', 0):>6.2f}  "
                f"{r.get('composite_score', 0):>6.3f}"
            )
        lines.append("-" * 60)

    lines.append("\nGenerated by DockIt V2")

    txt_path.write_text("\n".join(lines))
    logger.info("Text report generated: %s", txt_path)
    return txt_path


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

def generate_csv(results: list[dict], output_path: Path) -> Path:
    """Write results to a CSV file.

    Parameters
    ----------
    results : list[dict]
        Scored docking results.
    output_path : Path
        Where to write the CSV.

    Returns
    -------
    Path
        Path to the CSV file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "rank", "name", "smiles", "affinity_kcal", "composite_score",
        "qed", "logP", "MW", "tpsa", "hbd", "hba", "rotatable_bonds", "source",
        "docking_method", "novelty_score", "admet_score", "admet_color",
    ]

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for i, r in enumerate(results, 1):
            admet = r.get("admet", {})
            row = {
                "rank": i,
                "name": r.get("name", ""),
                "smiles": r.get("smiles", ""),
                "affinity_kcal": r.get("affinity", ""),
                "composite_score": r.get("composite_score", ""),
                "qed": r.get("qed", ""),
                "logP": r.get("logP", ""),
                "MW": r.get("MW", ""),
                "tpsa": r.get("tpsa", ""),
                "hbd": r.get("hbd", ""),
                "hba": r.get("hba", ""),
                "rotatable_bonds": r.get("rotatable_bonds", ""),
                "source": r.get("source", ""),
                "docking_method": r.get("docking_method", ""),
                "novelty_score": r.get("novelty_score", ""),
                "admet_score": admet.get("composite_score", "") if isinstance(admet, dict) else "",
                "admet_color": admet.get("color_code", "") if isinstance(admet, dict) else "",
            }
            writer.writerow(row)

    logger.info("CSV written: %s", output_path)
    return output_path


# ---------------------------------------------------------------------------
# ZIP archive
# ---------------------------------------------------------------------------

def generate_zip_archive(work_dir: Path, output_path: Path) -> Path:
    """Bundle all job output files into a ZIP archive.

    Includes PDB, PDBQT poses, PDF/TXT report, CSV, and ligand files.

    Parameters
    ----------
    work_dir : Path
        Job working directory containing all output files.
    output_path : Path
        Where to write the ZIP file.

    Returns
    -------
    Path
        Path to the ZIP archive.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    extensions = {".pdb", ".pdbqt", ".pdf", ".txt", ".csv", ".sdf"}

    with zipfile.ZipFile(output_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for fpath in sorted(work_dir.rglob("*")):
            if fpath.is_file() and fpath.suffix.lower() in extensions:
                arcname = fpath.relative_to(work_dir)
                zf.write(fpath, arcname)
                logger.debug("Added to ZIP: %s", arcname)

    size_mb = output_path.stat().st_size / (1024 * 1024)
    logger.info("ZIP archive generated: %s (%.1f MB)", output_path, size_mb)
    return output_path
