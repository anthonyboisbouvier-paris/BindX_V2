"""
DockIt — FastAPI application.

Provides the REST API for creating docking jobs, polling status,
and downloading results.
V2: advanced mode with generation, ADMET, retrosynthesis.
V3: rapid/standard/deep modes, sequence input, auto strategy, score_100,
    affinity_stars, toxicity, pedagogical tips, pipeline_summary.
V7: authentication (JWT), user accounts, project management.
"""

from __future__ import annotations

import asyncio
import json
import logging
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeoutError
from contextlib import asynccontextmanager
from pathlib import Path
from typing import AsyncGenerator, Optional
from uuid import uuid4

import requests as http_requests
from fastapi import Depends, FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse

from auth import create_token, decode_token, hash_password, verify_password
from database import (
    create_job, get_job, init_db, list_jobs,
    create_user, get_user, get_user_by_email,
    create_project, get_project, list_projects, update_project,
    count_project_jobs, list_project_jobs,
    create_assessment, get_assessment, get_assessment_by_hash, update_assessment_agent,
)
from models import (
    ADMETResult, DockingResult, JobCreate, JobResults, JobStatus,
    OffTargetResult, OptimizationRequest, OptimizationStatus,
    PoseQuality, ScaffoldAnalysisResponse,
    SynthesisRoute, SynthesisStep,
    UserCreate, UserLogin, UserResponse,
    ProjectCreate, ProjectUpdate, ProjectResponse,
    TargetAssessmentRequest, TargetAssessmentResult,
    AgentQuery, AgentResponse,
)
from pipeline.structure import query_rcsb_pdb, fetch_structure_from_sequence

# V9 routers
from routers.v9 import health as v9_health_router
from routers.v9 import projects as v9_projects_router
from routers.v9 import campaigns as v9_campaigns_router
from routers.v9 import phases as v9_phases_router
from routers.v9 import runs as v9_runs_router
from routers.v9 import molecules as v9_molecules_router

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
)
logger = logging.getLogger("dockit.api")

DATA_DIR = Path(os.environ.get("DOCKIT_DATA_DIR", "/data")) / "jobs"

# Thread pool for non-blocking Celery dispatch
_dispatch_pool = ThreadPoolExecutor(max_workers=2)

# Thread pool for preview-target blocking I/O (API calls, pocket detection)
_preview_pool = ThreadPoolExecutor(max_workers=4)


# ---------------------------------------------------------------------------
# V7: Auth dependency
# ---------------------------------------------------------------------------

async def get_current_user_id(request: Request) -> Optional[str]:
    """Extract user_id from JWT if present. Returns None for anonymous requests."""
    auth_header = request.headers.get("Authorization", "")
    if not auth_header.startswith("Bearer "):
        return None
    token = auth_header[7:]
    user_id = decode_token(token)
    return user_id


# ---------------------------------------------------------------------------
# V3 helpers
# ---------------------------------------------------------------------------

def _compute_score_100(composite_score: float) -> int:
    """Map a composite score (0.0-1.0 range) to an integer 0-100.

    Parameters
    ----------
    composite_score : float
        Raw composite score, typically between 0 and 1.

    Returns
    -------
    int
        Score mapped to 0-100, clamped.
    """
    return max(0, min(100, int(round(composite_score * 100))))


def _compute_affinity_stars(affinity: float) -> int:
    """Map a binding affinity (kcal/mol, negative is better) to 1-5 stars.

    Parameters
    ----------
    affinity : float
        Binding affinity in kcal/mol (more negative = stronger binding).

    Returns
    -------
    int
        Star rating from 1 (weak) to 5 (very strong).
    """
    if affinity < -10.0:
        return 5
    elif affinity < -8.0:
        return 4
    elif affinity < -6.0:
        return 3
    elif affinity < -4.0:
        return 2
    else:
        return 1


def _get_toxicity_info(admet_data: Optional[dict]) -> tuple[str, str]:
    """Derive toxicity level and color from ADMET data.

    Parameters
    ----------
    admet_data : dict or None
        ADMET result dictionary.

    Returns
    -------
    tuple[str, str]
        (toxicity_level, toxicity_color) where level is low/medium/high/unknown
        and color is a hex color code.
    """
    if not admet_data or not isinstance(admet_data, dict):
        return "unknown", "#9ca3af"

    # Check the color_code or composite_score from ADMET
    color_code = admet_data.get("color_code", "yellow")
    composite = admet_data.get("composite_score", 0.5)
    flags = admet_data.get("flags", [])

    # Derive from flags count and composite
    if color_code == "green" or (composite >= 0.7 and len(flags) == 0):
        return "low", "#22c55e"
    elif color_code == "red" or composite < 0.3 or len(flags) >= 3:
        return "high", "#ef4444"
    else:
        return "medium", "#f59e0b"


def _get_synthesis_info(synthesis_route: Optional[dict]) -> tuple[Optional[int], Optional[str], Optional[str]]:
    """Extract synthesis feasibility from a synthesis route dict.

    Parameters
    ----------
    synthesis_route : dict or None
        Synthesis route dictionary.

    Returns
    -------
    tuple[Optional[int], Optional[str], Optional[str]]
        (n_steps, feasibility, color) where feasibility is easy/moderate/hard/unknown.
    """
    if not synthesis_route or not isinstance(synthesis_route, dict):
        return None, None, None

    n_steps = synthesis_route.get("n_steps", 0)
    confidence = synthesis_route.get("confidence", 0.0)
    available = synthesis_route.get("all_reagents_available", False)

    if n_steps <= 3 and confidence >= 0.7 and available:
        return n_steps, "easy", "#22c55e"
    elif n_steps <= 6 and confidence >= 0.4:
        return n_steps, "moderate", "#f59e0b"
    elif n_steps > 0:
        return n_steps, "hard", "#ef4444"
    else:
        return None, None, None


def _get_pedagogical_tip(step: str) -> str:
    """Return an educational text explaining what is happening at each pipeline step.

    Parameters
    ----------
    step : str
        Current pipeline step name or keyword.

    Returns
    -------
    str
        A short pedagogical explanation in English.
    """
    tips = {
        "structure": (
            "The 3D structure of the protein is essential for molecular docking. "
            "It reveals the exact shape of the active site where candidate molecules "
            "will bind. We use AlphaFold (AI prediction) or experimental structures (X-ray, cryo-EM)."
        ),
        "pockets": (
            "Binding pockets are cavities on the protein surface where small molecules "
            "can bind. The fpocket algorithm analyzes the protein geometry to identify "
            "these cavities. The best pocket is selected for docking."
        ),
        "prepare": (
            "Receptor preparation converts the PDB structure to PDBQT format, "
            "which is required by AutoDock Vina. This step adds partial charges "
            "and atom types used in the affinity calculation."
        ),
        "ligands": (
            "Ligands are the candidate small molecules that will be tested against the protein. "
            "They come from databases such as ChEMBL (known molecules) or ZINC (drug-like molecules). "
            "The more studied a target is, the more reference compounds are available."
        ),
        "docking": (
            "Molecular docking simulates the interaction between each ligand and the protein. "
            "AutoDock Vina explores thousands of conformations to find the optimal position "
            "of the ligand in the pocket. The affinity (in kcal/mol) indicates the binding strength -- the more negative, the better."
        ),
        "scoring": (
            "Scoring combines binding affinity with physicochemical properties "
            "(Lipinski, QED, logP) to produce a composite score. A good drug candidate "
            "must have strong affinity AND good drug-like properties."
        ),
        "admet": (
            "ADMET stands for Absorption, Distribution, Metabolism, Excretion, and Toxicity. "
            "These predictions evaluate whether a molecule will be well absorbed by the body, "
            "whether it will reach its target, and whether it poses toxicity risks."
        ),
        "generation": (
            "Molecule generation uses artificial intelligence (REINVENT4) "
            "to create new molecules optimized for the target. "
            "This is useful when few known molecules exist for a protein."
        ),
        "retrosynthesis": (
            "Retrosynthesis plans how to synthesize a molecule in the laboratory, "
            "starting from the final product and working backwards to commercially available reagents. "
            "Fewer steps and more available reagents make the synthesis more feasible."
        ),
        "report": (
            "The final report compiles all results: top molecules, scores, "
            "properties, and visualizations. It can be downloaded as a PDF to be "
            "shared with collaborators or used to guide laboratory testing."
        ),
    }

    # Match the step to a tip key
    step_lower = step.lower()
    for key, tip in tips.items():
        if key in step_lower:
            return tip

    return (
        "The molecular docking pipeline systematically analyzes the interaction "
        "between a target protein and candidate molecules to identify "
        "potential drug compounds."
    )


# ---------------------------------------------------------------------------
# V3: mode to pipeline params mapping
# ---------------------------------------------------------------------------

def _map_mode_to_params(mode: str, body: JobCreate) -> dict:
    """Map request body to pipeline parameters.

    V11: mode is kept for backward compat but the key control is now
    `enable_dmpk`. The frontend sends compounds count, engine choice,
    and DMPK toggle directly.
    """
    enable_dmpk = body.enable_dmpk if hasattr(body, 'enable_dmpk') else True

    # V12: granular analysis flags (backward compat: enable_dmpk=False disables all)
    master = enable_dmpk
    enable_admet = getattr(body, 'enable_admet', master) if master else False
    enable_synthesis = getattr(body, 'enable_synthesis', master) if master else False
    enable_selectivity = getattr(body, 'enable_selectivity', master) if master else False
    enable_herg = getattr(body, 'enable_herg', master) if master else False
    enable_safety = getattr(body, 'enable_safety', master) if master else False
    box_size = getattr(body, 'box_size', None)

    base_params: dict = {
        "uniprot_id": body.uniprot_id or "",
        "sequence": body.sequence,
        "smiles_list": body.smiles_list,
        "mode": mode,
        "docking_engine": body.docking_engine,
        "notification_email": body.notification_email,
        "enable_dmpk": enable_dmpk,
        "enable_admet": enable_admet,
        "enable_synthesis": enable_synthesis,
        "enable_selectivity": enable_selectivity,
        "enable_herg": enable_herg,
        "enable_safety": enable_safety,
        "box_size": box_size,
        "max_ligands": body.max_ligands or 50,
        "use_chembl": body.use_chembl if body.use_chembl is not None else True,
        "use_zinc": body.use_zinc if body.use_zinc is not None else True,
        "enable_generation": False,
        "enable_diffdock": body.enable_diffdock if body.enable_diffdock is not None else False,
        "enable_retrosynthesis": enable_synthesis,
        "n_generated_molecules": 0,
        "auto_strategy": True,
    }

    # Deep screening for very large compound sets (>200)
    if (body.max_ligands or 0) > 200 or mode == "deep":
        base_params.update({
            "max_ligands": body.max_ligands or 5000,
            "enable_generation": enable_dmpk,
            "n_generated_molecules": body.n_generated_molecules or 200 if enable_dmpk else 0,
        })

    return base_params


# ---------------------------------------------------------------------------
# Lifespan (startup / shutdown)
# ---------------------------------------------------------------------------

@asynccontextmanager
async def lifespan(app: FastAPI) -> AsyncGenerator[None, None]:
    """Application lifespan: initialise DB on startup."""
    logger.info("DockIt API starting up (V3)")
    init_db()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    yield
    logger.info("DockIt API shutting down")
    _dispatch_pool.shutdown(wait=False)
    _preview_pool.shutdown(wait=False)


# ---------------------------------------------------------------------------
# FastAPI app
# ---------------------------------------------------------------------------

app = FastAPI(
    title="DockIt API",
    description="Virtual screening / molecular docking pipeline — V7.0 (Auth + Projects)",
    version="7.0.0",
    lifespan=lifespan,
)

# CORS — allow all origins for development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# V9 routers
app.include_router(v9_health_router.router, prefix="/api/v9", tags=["v9"])
app.include_router(v9_projects_router.router, prefix="/api/v9", tags=["v9-projects"])
app.include_router(v9_campaigns_router.router, prefix="/api/v9", tags=["v9-campaigns"])
app.include_router(v9_phases_router.router, prefix="/api/v9", tags=["v9-phases"])
app.include_router(v9_runs_router.router, prefix="/api/v9", tags=["v9-runs"])
app.include_router(v9_molecules_router.router, prefix="/api/v9", tags=["v9-molecules"])


# ---------------------------------------------------------------------------
# Helper: dispatch Celery task with timeout
# ---------------------------------------------------------------------------

def _dispatch_celery_task(job_id: str, task_params: dict, task_name: str = "run_pipeline") -> None:
    """Send the pipeline task to Celery, with a short timeout to avoid blocking.

    Parameters
    ----------
    job_id : str
        Job UUID.
    task_params : dict
        Parameters to pass to the Celery task.
    task_name : str
        Either "run_pipeline" or "run_deep_screening".
    """
    if task_name == "run_deep_screening":
        from tasks import run_deep_screening
        run_deep_screening.apply_async(
            args=[job_id, task_params],
            retry=False,
            retry_policy={"max_retries": 0},
        )
    else:
        from tasks import run_pipeline
        run_pipeline.apply_async(
            args=[job_id, task_params],
            retry=False,
            retry_policy={"max_retries": 0},
        )


# ---------------------------------------------------------------------------
# Helper: parse results JSON into DockingResult models (V3 enriched)
# ---------------------------------------------------------------------------

def _extract_admet_domain(admet_data: Optional[dict]) -> Optional[dict]:
    """Extract the applicability domain from ADMET data if present."""
    if not admet_data or not isinstance(admet_data, dict):
        return None
    ad = admet_data.get("applicability_domain")
    if isinstance(ad, dict):
        return ad
    return None


def _lazy_svg(smiles: Optional[str]) -> Optional[str]:
    """Generate 2D SVG on-the-fly for results that lack a pre-rendered SVG."""
    if not smiles:
        return None
    try:
        from pipeline.scoring import generate_2d_svg
        return generate_2d_svg(smiles)
    except Exception:
        return None


def _infer_docking_status(r: dict) -> str:
    """Infer docking status from result data."""
    engine = r.get("docking_engine") or r.get("docking_method") or ""
    has_pose = bool(r.get("pose_pdbqt_path") or r.get("pose_sdf_path"))
    if has_pose and engine.lower() in ("gnina", "vina", "autodock_vina"):
        return "docked"
    if engine.lower() in ("mock",):
        return "docked"  # mock still counts as having been through docking
    return "not_docked"


def _parse_results(raw_list: list[dict]) -> list[DockingResult]:
    """Convert raw JSON dicts into DockingResult pydantic models with V3 fields."""
    results: list[DockingResult] = []
    for r in raw_list:
        # Parse ADMET sub-object
        admet_data = r.get("admet")
        admet_obj = None
        if admet_data and isinstance(admet_data, dict):
            admet_obj = ADMETResult(
                oral_bioavailability=admet_data.get("absorption", {}).get("oral_bioavailability") if isinstance(admet_data.get("absorption"), dict) else admet_data.get("oral_bioavailability"),
                herg_inhibition=admet_data.get("toxicity", {}).get("herg_inhibition") if isinstance(admet_data.get("toxicity"), dict) else admet_data.get("herg_inhibition"),
                hepatotoxicity=admet_data.get("toxicity", {}).get("hepatotoxicity") if isinstance(admet_data.get("toxicity"), dict) else admet_data.get("hepatotoxicity"),
                ames_mutagenicity=admet_data.get("toxicity", {}).get("ames_mutagenicity") if isinstance(admet_data.get("toxicity"), dict) else admet_data.get("ames_mutagenicity"),
                bbb_permeability=admet_data.get("distribution", {}).get("bbb_permeability") if isinstance(admet_data.get("distribution"), dict) else admet_data.get("bbb_permeability"),
                plasma_protein_binding=admet_data.get("distribution", {}).get("plasma_protein_binding") if isinstance(admet_data.get("distribution"), dict) else admet_data.get("plasma_protein_binding"),
                composite_score=admet_data.get("composite_score", 0.5),
                flags=admet_data.get("flags", []),
                color_code=admet_data.get("color_code", "yellow"),
            )

        # Parse synthesis route
        synth_data = r.get("synthesis_route")
        synth_obj = None
        if synth_data and isinstance(synth_data, dict):
            steps = []
            for s in synth_data.get("steps", []):
                steps.append(SynthesisStep(
                    reaction=s.get("reaction", ""),
                    reactants=s.get("reactants", []),
                    reactant_names=s.get("reactant_names", []),
                    conditions=s.get("conditions", ""),
                    confidence=s.get("confidence", 0.0),
                ))
            synth_obj = SynthesisRoute(
                n_steps=synth_data.get("n_steps", 0),
                confidence=synth_data.get("confidence", 0.0),
                steps=steps,
                all_reagents_available=synth_data.get("all_reagents_available", False),
                estimated_cost=synth_data.get("estimated_cost", "unknown"),
                tree=synth_data.get("tree"),
                # V6.3 cost fields
                cost_estimate=synth_data.get("cost_estimate"),
                reagent_availability=synth_data.get("reagent_availability"),
            )

        # V3: compute derived fields
        affinity = r.get("affinity", 0.0)
        composite = r.get("composite_score", 0.0)
        score_100 = r.get("score_100", _compute_score_100(composite))
        affinity_stars = r.get("affinity_stars", _compute_affinity_stars(affinity))
        toxicity_level, toxicity_color = _get_toxicity_info(admet_data)
        synth_steps, synth_feasibility, synth_color = _get_synthesis_info(synth_data)

        results.append(DockingResult(
            name=r.get("name", "Unknown"),
            smiles=r.get("smiles", ""),
            affinity=affinity,
            logp=r.get("logP") or r.get("logp"),
            mw=r.get("MW") or r.get("mw"),
            tpsa=r.get("tpsa"),
            qed=r.get("qed"),
            hbd=r.get("hbd"),
            hba=r.get("hba"),
            rotatable_bonds=r.get("rotatable_bonds"),
            composite_score=composite,
            svg=r.get("svg_2d") or r.get("svg") or _lazy_svg(r.get("smiles")),
            pose_pdbqt=r.get("pose_pdbqt_path") or r.get("pose_sdf_path"),
            source=r.get("source"),
            admet=admet_obj,
            synthesis_route=synth_obj,
            docking_method=r.get("docking_method"),
            novelty_score=r.get("novelty_score"),
            # V3 fields
            score_100=score_100,
            affinity_stars=affinity_stars,
            toxicity_level=r.get("toxicity_level", toxicity_level),
            toxicity_color=r.get("toxicity_color", toxicity_color),
            synthesis_steps=r.get("synthesis_steps", synth_steps),
            synthesis_feasibility=r.get("synthesis_feasibility", synth_feasibility),
            synthesis_color=r.get("synthesis_color", synth_color),
            # V5 fields
            off_target=r.get("off_target_results"),
            confidence=r.get("confidence"),
            # V5bis fields
            vina_score=r.get("vina_score"),
            cnn_score=r.get("cnn_score"),
            cnn_affinity=r.get("cnn_affinity"),
            cnn_vs=r.get("cnn_vs") or (
                round((r.get("cnn_score") or 0) * (r.get("cnn_affinity") or 0), 4)
                if r.get("cnn_score") and r.get("cnn_affinity") else None
            ),
            consensus_rank=r.get("consensus_rank"),
            consensus_robust=r.get("consensus_robust"),
            interactions=r.get("interactions"),
            interaction_quality=r.get("interaction_quality"),
            cluster_id=r.get("cluster_id"),
            cluster_size=r.get("cluster_size"),
            is_representative=r.get("is_representative"),
            eliminated=r.get("eliminated"),
            elimination_reason=r.get("elimination_reason"),
            admet_domain=_extract_admet_domain(admet_data),
            pains_alert=r.get("pains_alert"),
            sa_score=r.get("sa_score"),
            # V6.1 fields
            consensus_detail=r.get("consensus_detail"),
            # V6.2 fields
            pareto_rank=r.get("pareto_rank"),
            pareto_front=r.get("pareto_front"),
            pareto_objectives=r.get("pareto_objectives"),
            # V6.3 fields
            combined_off_target=r.get("combined_off_target"),
            herg_specialized=r.get("herg_specialized"),
            # V11 fields: docking transparency
            docking_engine=r.get("docking_engine"),
            docking_status=r.get("docking_status") or _infer_docking_status(r),
            # V12 fields: pose quality
            pose_quality=r.get("pose_quality"),
        ))
    return results


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------

@app.get("/api/health")
async def health_check() -> dict:
    """Health-check endpoint."""
    return {"status": "ok", "service": "dockit-api", "version": "7.0.0"}


# ---------------------------------------------------------------------------
# V7: Auth endpoints
# ---------------------------------------------------------------------------

@app.post("/api/auth/register", response_model=dict)
async def register(payload: UserCreate):
    """Create a new user account."""
    existing = get_user_by_email(payload.email)
    if existing:
        raise HTTPException(status_code=409, detail="Email already registered")

    user_id = str(uuid4())
    pw_hash = hash_password(payload.password)
    user = create_user(user_id, payload.email, payload.username, pw_hash)
    token = create_token(user_id)

    return {
        "user": {
            "id": user.id,
            "email": user.email,
            "username": user.username,
        },
        "token": token,
    }


@app.post("/api/auth/login", response_model=dict)
async def login(payload: UserLogin):
    """Login and receive a JWT token."""
    user = get_user_by_email(payload.email)
    if not user or not verify_password(payload.password, user.password_hash):
        raise HTTPException(status_code=401, detail="Invalid email or password")

    token = create_token(user.id)
    return {
        "user": {
            "id": user.id,
            "email": user.email,
            "username": user.username,
        },
        "token": token,
    }


@app.get("/api/auth/me", response_model=dict)
async def get_me(user_id: Optional[str] = Depends(get_current_user_id)):
    """Get current user info from JWT."""
    if not user_id:
        raise HTTPException(status_code=401, detail="Not authenticated")

    user = get_user(user_id)
    if not user:
        raise HTTPException(status_code=404, detail="User not found")

    return {
        "id": user.id,
        "email": user.email,
        "username": user.username,
        "created_at": user.created_at.isoformat() if user.created_at else None,
    }


# ---------------------------------------------------------------------------
# V7: Project endpoints
# ---------------------------------------------------------------------------

@app.post("/api/projects", response_model=dict)
async def create_project_endpoint(
    payload: ProjectCreate,
    user_id: Optional[str] = Depends(get_current_user_id),
):
    """Create a new project (requires auth)."""
    if not user_id:
        raise HTTPException(status_code=401, detail="Authentication required to create projects")

    import json as _json
    project_id = str(uuid4())
    preview_json = _json.dumps(payload.target_preview_json) if payload.target_preview_json else None

    project = create_project(
        project_id=project_id,
        user_id=user_id,
        name=payload.name,
        uniprot_id=payload.uniprot_id,
        sequence=payload.sequence,
        description=payload.description,
        target_preview_json=preview_json,
    )

    return {
        "id": project.id,
        "name": project.name,
        "uniprot_id": project.uniprot_id,
        "description": project.description,
        "created_at": project.created_at.isoformat() if project.created_at else None,
    }


@app.get("/api/projects", response_model=list)
async def list_projects_endpoint(
    user_id: Optional[str] = Depends(get_current_user_id),
):
    """List current user's projects (requires auth)."""
    if not user_id:
        raise HTTPException(status_code=401, detail="Authentication required")

    import json as _json
    projects = list_projects(user_id)
    result = []
    for p in projects:
        job_count = count_project_jobs(p.id)
        preview = None
        if p.target_preview_json:
            try:
                preview = _json.loads(p.target_preview_json)
            except Exception:
                preview = None
        result.append({
            "id": p.id,
            "name": p.name,
            "uniprot_id": p.uniprot_id,
            "description": p.description,
            "target_preview_json": preview,
            "created_at": p.created_at.isoformat() if p.created_at else None,
            "job_count": job_count,
        })
    return result


def _job_best_score(job) -> float | None:
    """Extract best composite_score from a completed job's results_json."""
    if job.status != "completed" or not job.results_json:
        return None
    try:
        results = json.loads(job.results_json)
        scores = [r.get("composite_score", 0.0) for r in results if isinstance(r, dict)]
        return round(max(scores), 4) if scores else None
    except Exception:
        return None


@app.get("/api/projects/{project_id}", response_model=dict)
async def get_project_endpoint(
    project_id: str,
    user_id: Optional[str] = Depends(get_current_user_id),
):
    """Get project detail with its jobs list."""
    project = get_project(project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")

    # Only owner can see project details
    if project.user_id != user_id:
        raise HTTPException(status_code=403, detail="Access denied")

    import json as _json
    jobs = list_project_jobs(project_id)
    preview = None
    if project.target_preview_json:
        try:
            preview = _json.loads(project.target_preview_json)
        except Exception:
            preview = None

    return {
        "id": project.id,
        "name": project.name,
        "uniprot_id": project.uniprot_id,
        "sequence": project.sequence,
        "description": project.description,
        "target_preview_json": preview,
        "created_at": project.created_at.isoformat() if project.created_at else None,
        "jobs": [{
            "id": j.id,
            "status": j.status,
            "mode": j.mode,
            "progress": j.progress,
            "current_step": j.current_step,
            "protein_name": j.protein_name,
            "created_at": j.created_at.isoformat() if j.created_at else None,
            "completed_at": j.completed_at.isoformat() if j.completed_at else None,
            "enable_generation": bool(j.enable_generation) if hasattr(j, 'enable_generation') else False,
            "has_custom_smiles": bool(j.smiles_list) if hasattr(j, 'smiles_list') else False,
            "best_score": _job_best_score(j),
            "error_message": j.error_message if j.status == "failed" else None,
        } for j in jobs],
    }


@app.put("/api/projects/{project_id}", response_model=dict)
async def update_project_endpoint(
    project_id: str,
    payload: ProjectUpdate,
    user_id: Optional[str] = Depends(get_current_user_id),
):
    """Update a project."""
    if not user_id:
        raise HTTPException(status_code=401, detail="Authentication required")

    project = get_project(project_id)
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    if project.user_id != user_id:
        raise HTTPException(status_code=403, detail="Access denied")

    import json as _json
    updates = {}
    if payload.name is not None:
        updates["name"] = payload.name
    if payload.description is not None:
        updates["description"] = payload.description
    if payload.target_preview_json is not None:
        updates["target_preview_json"] = _json.dumps(payload.target_preview_json)

    if updates:
        update_project(project_id, **updates)

    return {"status": "updated", "id": project_id}


# ---------------------------------------------------------------------------
# Job endpoints
# ---------------------------------------------------------------------------

@app.get("/api/jobs")
async def list_all_jobs(
    limit: int = 50,
    project_id: Optional[str] = None,
    user_id: Optional[str] = Depends(get_current_user_id),
) -> list[dict]:
    """List recent jobs, ordered by creation date descending.

    Parameters
    ----------
    limit : int
        Maximum number of jobs to return (default 50, max 200).
    project_id : str, optional
        Filter by project ID (V7).
    user_id : str, optional
        Injected from JWT. Used to filter by user when project_id is given.

    Returns
    -------
    list[dict]
        List of job summaries with ``job_id``, ``uniprot_id``, ``mode``,
        ``status``, and ``created_at``.
    """
    limit = max(1, min(limit, 200))
    jobs = list_jobs(limit=limit, project_id=project_id, user_id=None)
    return [
        {
            "job_id": job.id,
            "uniprot_id": job.uniprot_id or None,
            "mode": job.mode or "rapid",
            "status": job.status,
            "created_at": job.created_at.isoformat() if job.created_at else None,
        }
        for job in jobs
    ]


# ---------------------------------------------------------------------------
# Preview target helpers (blocking I/O, run in thread pool)
# ---------------------------------------------------------------------------

def _fetch_protein_name(uniprot_id: str) -> Optional[str]:
    """Fetch the recommended protein name from the UniProt REST API.

    Parameters
    ----------
    uniprot_id : str
        UniProt accession (e.g. ``"P00533"``).

    Returns
    -------
    str or None
        The recommended full name, or None on failure.
    """
    try:
        resp = http_requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json",
            timeout=(5, 15),
        )
        resp.raise_for_status()
        data = resp.json()
        # Navigate: proteinDescription -> recommendedName -> fullName -> value
        protein_desc = data.get("proteinDescription", {})
        rec_name = protein_desc.get("recommendedName", {})
        full_name = rec_name.get("fullName", {})
        name = full_name.get("value")
        if name:
            return name
        # Fallback: submittedName
        submitted = protein_desc.get("submittedName", [])
        if submitted and isinstance(submitted, list):
            return submitted[0].get("fullName", {}).get("value")
        return None
    except Exception as exc:
        logger.warning("Failed to fetch protein name for %s: %s", uniprot_id, exc)
        return None


def _check_alphafold_availability(uniprot_id: str) -> bool:
    """Check whether AlphaFold DB has a prediction for this UniProt ID.

    Parameters
    ----------
    uniprot_id : str
        UniProt accession.

    Returns
    -------
    bool
        True if AlphaFold has a prediction available.
    """
    try:
        resp = http_requests.head(
            f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}",
            timeout=(5, 10),
        )
        return resp.status_code == 200
    except Exception:
        return False


def _fetch_chembl_info(uniprot_id: str) -> dict:
    """Fetch ChEMBL target info and activity counts for a UniProt accession.

    Parameters
    ----------
    uniprot_id : str
        UniProt accession (e.g. ``"P00533"``).

    Returns
    -------
    dict
        Dictionary with ``has_data``, ``n_actives``, ``n_with_ic50``,
        ``target_chembl_id``, and ``description`` keys.
    """
    result: dict = {
        "has_data": False,
        "n_actives": 0,
        "n_with_ic50": 0,
        "target_chembl_id": None,
        "description": "No ChEMBL data available for this target.",
    }

    try:
        # Step 1: Get target_chembl_id from UniProt accession
        target_resp = http_requests.get(
            "https://www.ebi.ac.uk/chembl/api/data/target.json",
            params={
                "target_components__accession": uniprot_id,
                "limit": 1,
                "format": "json",
            },
            timeout=(5, 15),
        )
        target_resp.raise_for_status()
        target_data = target_resp.json()

        targets = target_data.get("targets", [])
        if not targets:
            return result

        target_chembl_id = targets[0].get("target_chembl_id")
        if not target_chembl_id:
            return result

        result["target_chembl_id"] = target_chembl_id

        # Step 2: Count total activities for this target
        activity_resp = http_requests.get(
            "https://www.ebi.ac.uk/chembl/api/data/activity.json",
            params={
                "target_chembl_id": target_chembl_id,
                "limit": 0,
                "format": "json",
            },
            timeout=(5, 15),
        )
        activity_resp.raise_for_status()
        activity_data = activity_resp.json()

        total_activities = activity_data.get("page_meta", {}).get("total_count", 0)

        # Step 3: Count activities with IC50 data
        ic50_resp = http_requests.get(
            "https://www.ebi.ac.uk/chembl/api/data/activity.json",
            params={
                "target_chembl_id": target_chembl_id,
                "standard_type": "IC50",
                "limit": 0,
                "format": "json",
            },
            timeout=(5, 15),
        )
        ic50_resp.raise_for_status()
        ic50_data = ic50_resp.json()

        n_ic50 = ic50_data.get("page_meta", {}).get("total_count", 0)

        result["has_data"] = total_activities > 0
        result["n_actives"] = total_activities
        result["n_with_ic50"] = n_ic50

        if total_activities > 0:
            result["description"] = (
                f"{total_activities} known active compounds found in ChEMBL, "
                f"including {n_ic50} with IC50 data. "
                f"ChEMBL is recommended as the primary screening library for this target."
            )
        else:
            result["description"] = (
                "No known active compounds found in ChEMBL for this target. "
                "ZINC drug-like library or de novo generation is recommended."
            )

    except Exception as exc:
        logger.warning("ChEMBL lookup failed for %s: %s", uniprot_id, exc)
        result["description"] = f"ChEMBL lookup failed: {exc}"

    return result


def _build_structures_list(pdb_info: Optional[dict], uniprot_id: str) -> list[dict]:
    """Return all available structure options for the preview, in priority order.

    Returns a list with up to 3 entries: PDB (if found), AlphaFold (if available),
    ESMFold (always, as on-demand fallback). The first entry is marked recommended.
    """
    structures: list[dict] = []

    if pdb_info is not None:
        pdb_id = pdb_info.get("pdb_id", "unknown")
        resolution = pdb_info.get("resolution")
        method = pdb_info.get("method", "UNKNOWN")
        ligand_id = pdb_info.get("ligand_id")
        ligand_name = pdb_info.get("ligand_name")
        source = "pdb_holo" if ligand_id else "pdb_experimental"
        label = "PDB Holo" if ligand_id else "PDB Experimental"
        confidence = 0.98 if ligand_id else 0.90
        desc = f"Experimental {method} structure ({pdb_id})"
        if resolution:
            desc += f" at {resolution} Å"
        if ligand_id:
            lig = f"{ligand_id} ({ligand_name})" if ligand_name else ligand_id
            desc += f" — co-crystallized ligand {lig}"
        structures.append({
            "source": source,
            "label": label,
            "pdb_id": pdb_id,
            "resolution": resolution,
            "method": method,
            "ligand_id": ligand_id,
            "ligand_name": ligand_name,
            "confidence": confidence,
            "recommended": True,
            "explanation": desc,
            "download_url": pdb_info.get("download_url"),
        })

    has_alphafold = _check_alphafold_availability(uniprot_id)
    if has_alphafold:
        structures.append({
            "source": "alphafold",
            "label": "AlphaFold",
            "pdb_id": None,
            "resolution": None,
            "method": "AI prediction",
            "ligand_id": None,
            "ligand_name": None,
            "confidence": 0.85,
            "recommended": len(structures) == 0,
            "explanation": f"AlphaFold DB predicted structure for {uniprot_id}. Good for targets without experimental structure.",
            "download_url": f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb",
        })

    # ESMFold always available as on-demand option
    structures.append({
        "source": "esmfold",
        "label": "ESMFold",
        "pdb_id": None,
        "resolution": None,
        "method": "ESMFold (on-demand)",
        "ligand_id": None,
        "ligand_name": None,
        "confidence": 0.60,
        "recommended": len(structures) == 0,
        "explanation": "On-demand ESMFold structure prediction. Lower quality than experimental structures, but always available.",
        "download_url": None,
    })

    return structures


def _build_structure_info(pdb_info: Optional[dict], uniprot_id: str) -> dict:
    """Build the ``structure`` section of the preview response.

    Parameters
    ----------
    pdb_info : dict or None
        Result from ``query_rcsb_pdb()`` or None if no PDB found.
    uniprot_id : str
        UniProt accession for AlphaFold fallback check.

    Returns
    -------
    dict
        Structure information with ``source``, ``explanation``, and optional
        PDB metadata.
    """
    if pdb_info is not None:
        pdb_id = pdb_info.get("pdb_id", "unknown")
        resolution = pdb_info.get("resolution")
        method = pdb_info.get("method", "UNKNOWN")
        ligand_id = pdb_info.get("ligand_id")
        ligand_name = pdb_info.get("ligand_name")

        explanation_parts = [
            f"Experimental {method} structure from PDB ({pdb_id})"
        ]
        if resolution:
            explanation_parts[0] += f" at {resolution} A resolution"
        if ligand_id:
            lig_desc = ligand_id
            if ligand_name:
                lig_desc = f"{ligand_id} ({ligand_name})"
            explanation_parts.append(
                f"with co-crystallized ligand {lig_desc}"
            )
        explanation_parts.append(
            "This is the highest quality structure available."
        )

        return {
            "source": "pdb_experimental",
            "pdb_id": pdb_id,
            "resolution": resolution,
            "method": method,
            "ligand_id": ligand_id,
            "ligand_name": ligand_name,
            "explanation": ". ".join(explanation_parts) + ".",
        }

    # No PDB -- check AlphaFold
    has_alphafold = _check_alphafold_availability(uniprot_id)
    if has_alphafold:
        return {
            "source": "alphafold",
            "pdb_id": None,
            "resolution": None,
            "method": "AI prediction",
            "ligand_id": None,
            "ligand_name": None,
            "explanation": (
                f"No experimental structure found in PDB for {uniprot_id}. "
                "AlphaFold DB has a predicted structure available. "
                "AI-predicted structures are suitable for pocket detection "
                "but may be less accurate for docking near flexible loops."
            ),
        }

    return {
        "source": "none",
        "pdb_id": None,
        "resolution": None,
        "method": None,
        "ligand_id": None,
        "ligand_name": None,
        "explanation": (
            f"No experimental or AlphaFold structure found for {uniprot_id}. "
            "ESMFold on-demand folding will be attempted when the pipeline runs. "
            "This may take longer and the predicted structure may have lower confidence."
        ),
    }


def _build_pockets_info(
    pdb_info: Optional[dict],
    uniprot_id: str,
) -> list[dict]:
    """Detect pockets if a PDB structure can be quickly obtained.

    This downloads the PDB to a temporary directory, runs pocket detection,
    and returns formatted pocket information.

    Parameters
    ----------
    pdb_info : dict or None
        Result from ``query_rcsb_pdb()`` with ``download_url``.
    uniprot_id : str
        UniProt accession.

    Returns
    -------
    list[dict]
        List of pocket dicts with ``rank``, ``method``, ``probability``,
        ``residues_count``, ``center``, ``selected``, and ``explanation``.
    """
    if pdb_info is None:
        return []

    from pipeline.pockets import detect_pockets

    try:
        # Download PDB to a temporary directory
        download_url = pdb_info.get("download_url", "")
        if not download_url:
            return []

        with tempfile.TemporaryDirectory(prefix="dockit_preview_") as tmp_dir:
            tmp_path = Path(tmp_dir)
            pdb_path = tmp_path / f"{uniprot_id}.pdb"

            resp = http_requests.get(download_url, timeout=(5, 30))
            resp.raise_for_status()
            content = resp.text
            if "ATOM" not in content or len(content) < 100:
                logger.warning("Preview PDB download has no ATOM records")
                return []

            pdb_path.write_text(content)

            # Run pocket detection
            ligand_id = pdb_info.get("ligand_id")
            raw_pockets = detect_pockets(
                pdb_path=pdb_path,
                work_dir=tmp_path,
                ligand_id=ligand_id,
            )

            # Format pockets for the preview response
            formatted: list[dict] = []
            for i, pocket in enumerate(raw_pockets[:5]):  # Limit to top 5
                method = pocket.get("method", "unknown")
                probability = pocket.get("probability", 0.0)
                residues = pocket.get("residues", [])
                center = pocket.get("center", (0.0, 0.0, 0.0))

                # Build explanation based on method
                if method == "co-crystallized_ligand":
                    lig = pocket.get("ligand_id", ligand_id or "unknown")
                    explanation = (
                        f"Pocket defined by co-crystallized ligand {lig}. "
                        "This is an experimentally validated binding site."
                    )
                elif method == "p2rank":
                    explanation = (
                        f"Pocket detected by P2Rank ML model with "
                        f"{probability:.0%} probability."
                    )
                elif method == "p2rank_mock":
                    explanation = (
                        f"Pocket detected by P2Rank with "
                        f"{probability:.0%} probability."
                    )
                elif method == "fpocket":
                    explanation = (
                        f"Pocket detected by fpocket geometric algorithm "
                        f"(druggability score: {probability:.2f})."
                    )
                else:
                    explanation = (
                        f"Pocket detected by {method} method "
                        f"(score: {probability:.2f})."
                    )

                formatted.append({
                    "rank": i + 1,
                    "method": method,
                    "probability": round(probability, 4),
                    "residues_count": len(residues) if isinstance(residues, list) else 0,
                    "residues": residues[:50] if isinstance(residues, list) else [],
                    "center": list(center) if isinstance(center, tuple) else center,
                    "selected": (i == 0),  # First pocket is the selected one
                    "explanation": explanation,
                    "volume": round(pocket.get("volume", 0.0), 1),
                })

            return formatted

    except Exception as exc:
        logger.warning(
            "Pocket detection failed during preview for %s: %s",
            uniprot_id, exc,
        )
        return []


# ---------------------------------------------------------------------------
# Preview target endpoint
# ---------------------------------------------------------------------------

@app.post("/api/preview-target")
async def preview_target(body: dict) -> dict:
    """Preview protein structure and target information before running the pipeline.

    This endpoint fetches protein name, structure info, pocket detection, and
    ChEMBL compound availability without creating a job. It is designed to be
    fast and return partial results if some lookups fail.

    Parameters
    ----------
    body : dict
        Must contain ``"uniprot_id"`` (str).

    Returns
    -------
    dict
        Preview data with ``uniprot_id``, ``protein_name``, ``structure``,
        ``pockets``, and ``chembl_info`` sections.

    Raises
    ------
    HTTPException
        400 if ``uniprot_id`` is missing or empty.
    """
    uniprot_id = body.get("uniprot_id", "").strip().upper()
    if not uniprot_id:
        raise HTTPException(
            status_code=400,
            detail="uniprot_id is required",
        )

    loop = asyncio.get_event_loop()

    # Launch all blocking I/O tasks concurrently in the thread pool
    protein_name_future = loop.run_in_executor(
        _preview_pool, _fetch_protein_name, uniprot_id,
    )
    pdb_info_future = loop.run_in_executor(
        _preview_pool, query_rcsb_pdb, uniprot_id,
    )
    chembl_future = loop.run_in_executor(
        _preview_pool, _fetch_chembl_info, uniprot_id,
    )

    # Await the independent lookups concurrently
    protein_name, pdb_info, chembl_info = await asyncio.gather(
        protein_name_future,
        pdb_info_future,
        chembl_future,
        return_exceptions=True,
    )

    # Handle exceptions from gather -- replace with safe defaults
    if isinstance(protein_name, BaseException):
        logger.warning("Protein name lookup raised: %s", protein_name)
        protein_name = None
    if isinstance(pdb_info, BaseException):
        logger.warning("RCSB PDB lookup raised: %s", pdb_info)
        pdb_info = None
    if isinstance(chembl_info, BaseException):
        logger.warning("ChEMBL lookup raised: %s", chembl_info)
        chembl_info = {
            "has_data": False,
            "n_actives": 0,
            "n_with_ic50": 0,
            "target_chembl_id": None,
            "description": f"ChEMBL lookup failed: {chembl_info}",
        }

    # Build all available structure options (PDB, AlphaFold, ESMFold)
    structures_list = _build_structures_list(pdb_info, uniprot_id)
    # Keep backward-compat single `structure` field = first (recommended) option
    structure_info = structures_list[0] if structures_list else _build_structure_info(pdb_info, uniprot_id)

    # Pocket detection requires downloading the PDB -- run in executor
    pockets_info: list[dict] = []
    try:
        pockets_info = await loop.run_in_executor(
            _preview_pool, _build_pockets_info, pdb_info, uniprot_id,
        )
    except Exception as exc:
        logger.warning("Pocket detection raised during preview: %s", exc)

    return {
        "uniprot_id": uniprot_id,
        "protein_name": protein_name or uniprot_id,
        "structure": structure_info,
        "structures": structures_list,
        "pockets": pockets_info,
        "chembl_info": chembl_info,
    }


@app.post("/api/preview-sequence")
async def preview_sequence(body: dict) -> dict:
    """Fold a raw amino-acid sequence with ESMFold and detect binding pockets."""
    sequence = body.get("sequence", "").strip().upper()
    if not sequence or len(sequence) < 10:
        raise HTTPException(400, "sequence must be at least 10 amino acids")

    loop = asyncio.get_event_loop()
    work_dir = Path("/tmp/dockit_preview_sequences")
    work_dir.mkdir(parents=True, exist_ok=True)

    # 1. Fold with ESMFold (blocking, run in thread pool)
    try:
        pdb_path, source = await loop.run_in_executor(
            _preview_pool, fetch_structure_from_sequence, sequence, work_dir
        )
    except Exception as exc:
        logger.error("ESMFold failed for sequence preview: %s", exc)
        raise HTTPException(500, f"Structure prediction failed: {exc}")

    # 2. Read PDB content to send to frontend (for 3D viewer)
    pdb_data = pdb_path.read_text() if pdb_path and pdb_path.exists() else ""

    # 3. Pocket detection on the folded PDB (run directly on local path)
    pockets_info: list[dict] = []
    try:
        def _detect_pockets_local() -> list[dict]:
            from pipeline.pockets import detect_pockets
            raw = detect_pockets(pdb_path=pdb_path, work_dir=pdb_path.parent, ligand_id=None)
            formatted: list[dict] = []
            for i, pocket in enumerate(raw[:5]):
                method = pocket.get("method", "unknown")
                probability = float(pocket.get("probability", 0.0))
                residues = pocket.get("residues", [])
                center = pocket.get("center", (0.0, 0.0, 0.0))
                formatted.append({
                    "rank": i + 1,
                    "method": method,
                    "probability": round(probability, 4),
                    "n_residues": len(residues) if isinstance(residues, list) else 0,
                    "volume": pocket.get("volume"),
                    "center": list(center) if isinstance(center, tuple) else center,
                    "selected": (i == 0),
                    "explanation": f"Pocket detected by {method} (score: {probability:.2f}).",
                })
            return formatted

        pockets_info = await loop.run_in_executor(_preview_pool, _detect_pockets_local)
    except Exception as exc:
        logger.warning("Pocket detection failed for sequence preview: %s", exc)

    structure = {
        "source": source,
        "label": "ESMFold" if source == "esmfold" else "Mock",
        "pdb_id": None,
        "resolution": None,
        "method": "ESMFold (on-demand)",
        "ligand_id": None,
        "ligand_name": None,
        "confidence": 0.60,
        "recommended": True,
        "explanation": f"Structure predicted from sequence ({len(sequence)} aa) via ESMFold.",
        "pdb_data": pdb_data,
    }

    return {
        "uniprot_id": None,
        "protein_name": f"Custom sequence ({len(sequence)} aa)",
        "sequence": sequence,
        "structure": structure,
        "structures": [structure],
        "pockets": pockets_info,
        "chembl_info": {"has_data": False, "n_actives": 0},
    }


@app.post("/api/jobs")
async def create_docking_job(
    body: JobCreate,
    user_id: Optional[str] = Depends(get_current_user_id),
) -> dict:
    """Create a new docking job and dispatch it to the Celery worker.

    V3: accepts sequence input, mode rapid/standard/deep, auto-determines
    ligand strategy. Deep mode dispatches run_deep_screening task.
    V7: optionally associates job with a project and user.
    """
    job_id = str(uuid4())

    # Normalize mode: map V2 names to V3
    mode = body.mode
    if mode == "basic":
        mode = "rapid"
    elif mode == "advanced":
        mode = "standard"

    # Build pipeline parameters based on mode
    task_params = _map_mode_to_params(mode, body)

    # V8: inject target_config to skip structure+pockets recomputation
    target_config_raw: Optional[str] = body.target_config_json
    if not target_config_raw and body.project_id:
        from database import get_project as _get_project
        proj = _get_project(body.project_id)
        if proj and getattr(proj, "target_preview_json", None):
            target_config_raw = proj.target_preview_json
    if target_config_raw:
        try:
            task_params["target_config"] = (
                json.loads(target_config_raw)
                if isinstance(target_config_raw, str)
                else target_config_raw
            )
        except Exception:
            pass  # malformed JSON — pipeline will re-compute

    # Determine effective boolean flags for DB storage
    use_chembl = task_params.get("use_chembl", True)
    use_zinc = task_params.get("use_zinc", False)
    max_ligands = task_params.get("max_ligands", 50)
    enable_gen = task_params.get("enable_generation", False)
    enable_dd = task_params.get("enable_diffdock", False)
    enable_retro = task_params.get("enable_retrosynthesis", False)
    n_gen = task_params.get("n_generated_molecules", 100)

    # Persist in DB
    create_job(
        job_id=job_id,
        uniprot_id=body.uniprot_id,
        use_chembl=use_chembl,
        use_zinc=use_zinc,
        max_ligands=max_ligands,
        smiles_list=body.smiles_list,
        mode=mode,
        enable_generation=enable_gen,
        enable_diffdock=enable_dd,
        enable_retrosynthesis=enable_retro,
        n_generated_molecules=n_gen,
        # V3 fields
        sequence=body.sequence,
        notification_email=body.notification_email,
        docking_engine=body.docking_engine,
        # V7 fields
        project_id=body.project_id,
        user_id=user_id,
    )

    # Choose task type based on mode
    task_name = "run_deep_screening" if mode == "deep" else "run_pipeline"

    try:
        future = _dispatch_pool.submit(_dispatch_celery_task, job_id, task_params, task_name)
        future.result(timeout=5)  # Wait at most 5 seconds for Redis
        display_id = body.uniprot_id or f"seq[{len(body.sequence or '')}aa]"
        logger.info("Job %s dispatched for %s (mode=%s, task=%s)", job_id, display_id, mode, task_name)
    except FuturesTimeoutError:
        logger.warning("Celery dispatch timed out for job %s (Redis may be offline)", job_id)
    except Exception as exc:
        logger.warning("Celery dispatch failed for job %s: %s", job_id, exc)

    return {"job_id": job_id, "mode": mode}


@app.get("/api/jobs/{job_id}")
async def get_job_status(job_id: str) -> JobStatus:
    """Return the current status and progress of a job (V3 enriched)."""
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    # V3: parse pipeline summary from JSON if stored
    pipeline_summary = None
    pipeline_summary_json = getattr(job, "pipeline_summary_json", None)
    if pipeline_summary_json:
        try:
            pipeline_summary = json.loads(pipeline_summary_json)
        except (json.JSONDecodeError, TypeError):
            pass

    # V3: build results summary from stored results
    results_summary = None
    if job.status == "completed" and job.results_json:
        try:
            raw_results = json.loads(job.results_json)
            results_summary = {
                "total_results": len(raw_results),
                "best_affinity": min((r.get("affinity", 0.0) for r in raw_results), default=0.0),
                "best_score": max((r.get("composite_score", 0.0) for r in raw_results), default=0.0),
            }
        except (json.JSONDecodeError, TypeError):
            pass

    # V3: pedagogical tip based on current step
    current_step = job.current_step or "Queued"
    pedagogical_tip = _get_pedagogical_tip(current_step)

    job_mode = getattr(job, "mode", "rapid") or "rapid"
    # Normalize V2 mode names
    if job_mode == "basic":
        job_mode = "rapid"
    elif job_mode == "advanced":
        job_mode = "standard"

    return JobStatus(
        job_id=job.id,
        status=job.status,
        progress=job.progress,
        current_step=current_step,
        error_message=job.error_message,
        created_at=job.created_at.isoformat() if job.created_at else None,
        completed_at=job.completed_at.isoformat() if job.completed_at else None,
        mode=job_mode,
        # V3 fields
        step_details=getattr(job, "current_step", None),
        strategy_message=getattr(job, "strategy_message", None),
        pedagogical_tip=pedagogical_tip,
        structure_source=getattr(job, "structure_source", None),
        pipeline_summary=pipeline_summary,
        results_summary=results_summary,
    )


@app.get("/api/jobs/{job_id}/results")
async def get_job_results(job_id: str) -> JobResults:
    """Return the scored docking results for a completed job (V3 enriched)."""
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    if job.status != "completed":
        raise HTTPException(
            status_code=400,
            detail=f"Job is {job.status}; results are only available for completed jobs",
        )

    # Parse known results from JSON
    results: list[DockingResult] = []
    if job.results_json:
        try:
            raw = json.loads(job.results_json)
            results = _parse_results(raw)
        except (json.JSONDecodeError, TypeError, KeyError) as exc:
            logger.error("Failed to deserialise results for job %s: %s", job_id, exc)

    # Parse generated molecules (V2)
    generated: list[DockingResult] = []
    generated_json = getattr(job, "generated_json", None)
    if generated_json:
        try:
            raw_gen = json.loads(generated_json)
            generated = _parse_results(raw_gen)
        except (json.JSONDecodeError, TypeError, KeyError) as exc:
            logger.error("Failed to deserialise generated results for job %s: %s", job_id, exc)

    job_mode = getattr(job, "mode", "rapid") or "rapid"
    if job_mode == "basic":
        job_mode = "rapid"
    elif job_mode == "advanced":
        job_mode = "standard"

    # V3: parse pipeline summary
    pipeline_summary = None
    pipeline_summary_json = getattr(job, "pipeline_summary_json", None)
    if pipeline_summary_json:
        try:
            pipeline_summary = json.loads(pipeline_summary_json)
        except (json.JSONDecodeError, TypeError):
            pass

    return JobResults(
        job_id=job_id,
        protein_name=job.protein_name or job.uniprot_id or "Unknown",
        uniprot_id=job.uniprot_id or None,
        results=results,
        generated_molecules=generated,
        pdb_file_url=f"/api/jobs/{job_id}/protein" if job.pdb_path else None,
        report_pdf_url=f"/api/jobs/{job_id}/report" if job.report_path else None,
        zip_url=f"/api/jobs/{job_id}/download" if job.zip_path else None,
        mode=job_mode,
        pipeline_summary=pipeline_summary,
    )


@app.get("/api/jobs/{job_id}/synthesis/{mol_index}")
async def get_synthesis_route(job_id: str, mol_index: int) -> dict:
    """Return detailed synthesis route for a specific molecule."""
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    # Combine known + generated results
    all_results = []
    if job.results_json:
        try:
            all_results.extend(json.loads(job.results_json))
        except json.JSONDecodeError:
            pass
    generated_json = getattr(job, "generated_json", None)
    if generated_json:
        try:
            all_results.extend(json.loads(generated_json))
        except json.JSONDecodeError:
            pass

    if mol_index < 0 or mol_index >= len(all_results):
        raise HTTPException(
            status_code=404,
            detail=f"Molecule index {mol_index} out of range (0-{len(all_results)-1})",
        )

    mol = all_results[mol_index]
    route = mol.get("synthesis_route")
    if not route:
        raise HTTPException(status_code=404, detail="No synthesis route available for this molecule")

    return {
        "molecule_name": mol.get("name", "Unknown"),
        "smiles": mol.get("smiles", ""),
        "synthesis_route": route,
    }


@app.get("/api/jobs/{job_id}/report")
async def download_report(job_id: str) -> FileResponse:
    """Download the PDF (or text) report for a completed job."""
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    if not job.report_path:
        raise HTTPException(status_code=404, detail="Report not available")

    report_path = Path(job.report_path)
    if not report_path.exists():
        # Try .txt fallback
        txt_path = report_path.with_suffix(".txt")
        if txt_path.exists():
            report_path = txt_path
        else:
            raise HTTPException(status_code=404, detail="Report file not found on disk")

    media_type = "application/pdf" if report_path.suffix == ".pdf" else "text/plain"
    return FileResponse(
        path=str(report_path),
        media_type=media_type,
        filename=f"dockit_report_{job_id[:8]}{report_path.suffix}",
    )


@app.get("/api/jobs/{job_id}/download")
async def download_zip(job_id: str) -> FileResponse:
    """Download the ZIP archive of all results."""
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    if not job.zip_path:
        raise HTTPException(status_code=404, detail="ZIP archive not available")

    zip_path = Path(job.zip_path)
    if not zip_path.exists():
        raise HTTPException(status_code=404, detail="ZIP file not found on disk")

    return FileResponse(
        path=str(zip_path),
        media_type="application/zip",
        filename=f"dockit_results_{job_id[:8]}.zip",
    )


@app.get("/api/jobs/{job_id}/protein")
async def download_protein(job_id: str) -> FileResponse:
    """Download the PDB file of the target protein."""
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    if not job.pdb_path:
        raise HTTPException(status_code=404, detail="PDB file not available")

    pdb_path = Path(job.pdb_path)
    if not pdb_path.exists():
        raise HTTPException(status_code=404, detail="PDB file not found on disk")

    return FileResponse(
        path=str(pdb_path),
        media_type="chemical/x-pdb",
        filename=f"{job.uniprot_id or 'protein'}.pdb",
    )


@app.get("/api/jobs/{job_id}/pose/{ligand_index}")
async def download_pose(job_id: str, ligand_index: int) -> FileResponse:
    """Download the docked PDBQT pose for a specific ligand (by index)."""
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    if not job.results_json:
        raise HTTPException(status_code=404, detail="No results available")

    try:
        results = json.loads(job.results_json)
    except json.JSONDecodeError:
        raise HTTPException(status_code=500, detail="Corrupt results data")

    if ligand_index < 0 or ligand_index >= len(results):
        raise HTTPException(
            status_code=404,
            detail=f"Ligand index {ligand_index} out of range (0-{len(results)-1})",
        )

    pose_path_str = (
        results[ligand_index].get("pose_pdbqt_path")
        or results[ligand_index].get("pose_sdf_path")
    )
    if not pose_path_str:
        raise HTTPException(status_code=404, detail="Pose file not available for this ligand")

    pose_path = Path(pose_path_str)
    if not pose_path.exists():
        raise HTTPException(status_code=404, detail="Pose file not found on disk")

    # Detect actual format from content (file extension may be wrong,
    # e.g. .pdbqt extension but SDF data from GNINA)
    try:
        content = pose_path.read_text()
    except Exception:
        content = ""

    if "V2000" in content or "V3000" in content or "$$$$" in content:
        media_type = "chemical/x-mdl-sdfile"
        filename_ext = ".sdf"
    elif "BRANCH" in content or "ENDROOT" in content or "TORSDOF" in content:
        media_type = "chemical/x-pdbqt"
        filename_ext = ".pdbqt"
    elif "ATOM" in content or "HETATM" in content:
        media_type = "chemical/x-pdb"
        filename_ext = ".pdb"
    else:
        media_type = "application/octet-stream"
        filename_ext = pose_path.suffix

    ligand_name = results[ligand_index].get("name", f"ligand_{ligand_index}")
    return FileResponse(
        path=str(pose_path),
        media_type=media_type,
        filename=f"{ligand_name}_pose{filename_ext}",
    )


# ---------------------------------------------------------------------------
# V5 Endpoints: Audit log, Lead optimization
# ---------------------------------------------------------------------------

# In-memory store for optimization results (keyed by optimization_id)
# TODO: persist to DB — currently lost on process restart
_optimization_results: dict[str, dict] = {}


@app.get("/api/jobs/{job_id}/audit_log")
async def get_audit_log(job_id: str) -> dict:
    """Return the audit log JSON for a completed (or running) job.

    The audit log is stored as ``audit_log.json`` in the job's
    working directory.
    """
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    audit_path = DATA_DIR / job_id / "audit_log.json"
    if not audit_path.exists():
        raise HTTPException(
            status_code=404,
            detail="Audit log not available for this job",
        )

    try:
        with open(audit_path, "r", encoding="utf-8") as f:
            audit_data = json.load(f)
        return audit_data
    except (json.JSONDecodeError, IOError) as exc:
        logger.error("Failed to read audit log for job %s: %s", job_id, exc)
        raise HTTPException(
            status_code=500,
            detail="Failed to read audit log file",
        )


@app.post("/api/molecule/analyze-scaffold")
async def analyze_scaffold_endpoint(body: dict) -> ScaffoldAnalysisResponse:
    """Analyze a molecule's scaffold, R-groups, and return annotated SVG.

    Parameters
    ----------
    body : dict
        Must contain ``smiles`` key.

    Returns
    -------
    ScaffoldAnalysisResponse
    """
    smiles = body.get("smiles", "")
    if not smiles:
        raise HTTPException(status_code=422, detail="smiles field is required")

    try:
        from pipeline.scaffold_analysis import analyze_scaffold
        result = analyze_scaffold(smiles)
        return ScaffoldAnalysisResponse(**result)
    except Exception as exc:
        logger.error("Scaffold analysis failed for %s: %s", smiles[:60], exc)
        raise HTTPException(status_code=500, detail=f"Scaffold analysis failed: {exc}")


@app.post("/api/jobs/{job_id}/optimize")
async def start_optimization(job_id: str, body: OptimizationRequest) -> dict:
    """Start a lead optimization task for a molecule from a completed job.

    Dispatches the optimization to a Celery task and returns an
    optimization_id for polling.

    Parameters
    ----------
    job_id : str
        UUID of the parent docking job (must be completed).
    body : OptimizationRequest
        Optimization parameters: smiles, molecule_name, weights,
        n_iterations, variants_per_iter.

    Returns
    -------
    dict
        ``{"optimization_id": str, "status": "queued"}``
    """
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    if job.status != "completed":
        raise HTTPException(
            status_code=400,
            detail=f"Job is {job.status}; optimization requires a completed job",
        )

    opt_id = str(uuid4())

    # Build optimization parameters
    opt_params: dict = {
        "job_id": job_id,
        "optimization_id": opt_id,
        "smiles": body.smiles,
        "molecule_name": body.molecule_name,
        "weights": body.weights,
        "n_iterations": body.n_iterations,
        "variants_per_iter": body.variants_per_iter,
        "structural_rules": body.modification_rules.model_dump() if body.modification_rules else None,
        "docking_engine": body.docking_engine,
        "dock_top_n": body.dock_top_n,
        "exhaustiveness": body.exhaustiveness,
        # V12: granular analysis flags
        "enable_admet": body.enable_admet,
        "enable_synthesis": body.enable_synthesis,
        "enable_selectivity": body.enable_selectivity,
        "enable_herg": body.enable_herg,
        "enable_safety": body.enable_safety,
        "box_size": body.box_size,
    }

    # Store initial status
    _optimization_results[opt_id] = {
        "optimization_id": opt_id,
        "job_id": job_id,
        "status": "queued",
        "progress": 0,
        "current_iteration": 0,
        "total_iterations": body.n_iterations,
        "best_score": None,
        "result": None,
        "error_message": None,
    }

    # Try to dispatch via Celery; fall back to synchronous execution
    try:
        future = _dispatch_pool.submit(
            _dispatch_optimization_task, opt_id, opt_params,
        )
        future.result(timeout=5)
        logger.info(
            "Optimization %s dispatched for job %s, molecule %s",
            opt_id, job_id, body.molecule_name,
        )
    except FuturesTimeoutError:
        logger.warning(
            "Optimization dispatch timed out for %s (will try sync)", opt_id,
        )
        # Run synchronously as fallback
        _run_optimization_sync(opt_id, opt_params)
    except Exception as exc:
        logger.warning(
            "Optimization dispatch failed for %s: %s (will try sync)",
            opt_id, exc,
        )
        _run_optimization_sync(opt_id, opt_params)

    return {"optimization_id": opt_id, "status": "queued"}


@app.get("/api/jobs/{job_id}/optimization/{opt_id}")
async def get_optimization_status(job_id: str, opt_id: str) -> OptimizationStatus:
    """Return the status and results of a lead optimization task.

    Parameters
    ----------
    job_id : str
        UUID of the parent docking job.
    opt_id : str
        UUID of the optimization task.

    Returns
    -------
    OptimizationStatus
        Current status, progress, and results (if completed).
    """
    opt_data = _optimization_results.get(opt_id)

    # If not in memory, try to recover from disk (e.g. after API restart)
    if opt_data is None:
        result_path = DATA_DIR / job_id / "optimization" / opt_id / "result.json"
        if result_path.exists():
            try:
                with open(result_path, "r", encoding="utf-8") as f:
                    result = json.load(f)
                # Try to recover created_job_id from meta file
                created_job_id = None
                meta_path = DATA_DIR / job_id / "optimization" / opt_id / "meta.json"
                if meta_path.exists():
                    try:
                        meta = json.loads(meta_path.read_text())
                        created_job_id = meta.get("created_job_id")
                    except Exception:
                        pass
                opt_data = {
                    "optimization_id": opt_id,
                    "job_id": job_id,
                    "status": "completed",
                    "progress": 100,
                    "current_iteration": result.get("total_tested", 10),
                    "total_iterations": len(result.get("iterations", [])) or 10,
                    "best_score": result.get("final_lead", {}).get("score")
                                 or result.get("best_molecule", {}).get("score"),
                    "result": result,
                    "error_message": None,
                    "created_job_id": created_job_id,
                }
                _optimization_results[opt_id] = opt_data
            except (json.JSONDecodeError, IOError):
                pass

    if opt_data is None:
        raise HTTPException(
            status_code=404,
            detail=f"Optimization {opt_id} not found",
        )

    if opt_data.get("job_id") != job_id:
        raise HTTPException(
            status_code=404,
            detail=f"Optimization {opt_id} does not belong to job {job_id}",
        )

    # If still queued/running, check if Celery worker wrote result.json
    if opt_data.get("status") in ("queued", "running"):
        result_path = DATA_DIR / job_id / "optimization" / opt_id / "result.json"
        if result_path.exists():
            try:
                with open(result_path, "r", encoding="utf-8") as f:
                    result = json.load(f)
                # Also recover created_job_id from meta.json
                created_job_id = opt_data.get("created_job_id")
                if not created_job_id:
                    meta_path = DATA_DIR / job_id / "optimization" / opt_id / "meta.json"
                    if meta_path.exists():
                        try:
                            meta = json.loads(meta_path.read_text())
                            created_job_id = meta.get("created_job_id")
                        except Exception:
                            pass
                opt_data.update({
                    "status": "completed",
                    "progress": 100,
                    "result": result,
                    "best_score": result.get("final_lead", {}).get("score")
                                 or result.get("best_molecule", {}).get("score"),
                    "current_iteration": opt_data.get("total_iterations", 10),
                    "created_job_id": created_job_id,
                })
            except (json.JSONDecodeError, IOError):
                pass

    return OptimizationStatus(
        optimization_id=opt_data["optimization_id"],
        job_id=opt_data["job_id"],
        status=opt_data.get("status", "unknown"),
        progress=opt_data.get("progress", 0),
        current_iteration=opt_data.get("current_iteration", 0),
        total_iterations=opt_data.get("total_iterations", 10),
        best_score=opt_data.get("best_score"),
        result=opt_data.get("result"),
        error_message=opt_data.get("error_message"),
        created_job_id=opt_data.get("created_job_id"),
    )


# ---------------------------------------------------------------------------
# Optimization dispatch helpers
# ---------------------------------------------------------------------------

def _dispatch_optimization_task(opt_id: str, params: dict) -> None:
    """Send the optimization task to Celery.

    Parameters
    ----------
    opt_id : str
        Optimization UUID.
    params : dict
        Optimization parameters.
    """
    try:
        from tasks import run_lead_optimization
        run_lead_optimization.apply_async(
            args=[opt_id, params],
            retry=False,
            retry_policy={"max_retries": 0},
        )
    except Exception as exc:
        logger.warning(
            "Celery dispatch for optimization %s failed: %s (falling back to sync)",
            opt_id, exc,
        )
        _run_optimization_sync(opt_id, params)


def _create_job_from_optimization(opt_id: str, params: dict, result: dict) -> Optional[str]:
    """Create a new job in the project from optimization results.

    Delegates to the same function in tasks.py for consistency (single
    implementation of the enrichment pipeline integration).

    Parameters
    ----------
    opt_id : str
        Optimization UUID.
    params : dict
        Original optimization parameters (contains ``job_id``).
    result : dict
        Optimization result dict with ``final_lead``, ``best_molecule``,
        ``iterations``, etc.

    Returns
    -------
    str or None
        The new job_id, or None on failure.
    """
    try:
        from tasks import _create_job_from_optimization as _tasks_create
        return _tasks_create(opt_id, params, result)
    except Exception as exc:
        logger.warning("Failed to create job from optimization %s: %s", opt_id, exc)
        return None


def _run_optimization_sync(opt_id: str, params: dict) -> None:
    """Run optimization synchronously (fallback when Celery is unavailable).

    Parameters
    ----------
    opt_id : str
        Optimization UUID.
    params : dict
        Optimization parameters.
    """
    import traceback

    try:
        _optimization_results[opt_id]["status"] = "running"

        from pipeline.lead_optimization import run_optimization

        job_id = params["job_id"]
        work_dir = DATA_DIR / job_id

        # Get pocket center and receptor path from the job's pipeline summary
        job = get_job(job_id)
        pocket_center = [22.0, 0.5, 18.0]  # fallback default
        target_pdbqt = None
        if job and job.pipeline_summary_json:
            try:
                summary = json.loads(job.pipeline_summary_json)
                pocket_center = summary.get("best_pocket_center", pocket_center)
                # V12: resolve receptor path from pipeline summary
                rp = summary.get("receptor_path", "")
                if rp and Path(rp).exists():
                    target_pdbqt = rp
            except (json.JSONDecodeError, TypeError):
                pass

        # V12: resolve receptor path (same logic as celery task)
        if not target_pdbqt:
            receptor_candidates = sorted(work_dir.glob("*receptor*.pdbqt"))
            if receptor_candidates:
                target_pdbqt = str(receptor_candidates[0])
            else:
                target_pdbqt = str(work_dir / "receptor.pdbqt")
                logger.warning("No receptor found for opt %s, using default: %s", opt_id, target_pdbqt)

        def progress_cb(iteration: int, score: float, message: str) -> None:
            n_iters = params.get("n_iterations", 10)
            progress = int(iteration / n_iters * 100)
            _optimization_results[opt_id].update({
                "progress": progress,
                "current_iteration": iteration,
                "best_score": round(score, 4),
            })

        result = run_optimization(
            starting_smiles=params["smiles"],
            starting_name=params.get("molecule_name", "molecule"),
            target_pdbqt=target_pdbqt,
            pocket_center=pocket_center,
            weights=params.get("weights"),
            n_iterations=params.get("n_iterations", 10),
            variants_per_iter=params.get("variants_per_iter", 50),
            work_dir=work_dir / "optimization" / opt_id,
            progress_callback=progress_cb,
            structural_rules=params.get("structural_rules"),
            docking_engine=params.get("docking_engine", "auto"),
            pocket_size=params.get("pocket_size", [25.0, 25.0, 25.0]),
            exhaustiveness=params.get("exhaustiveness", 8),
            dock_top_n=params.get("dock_top_n", 20),
        )

        # Auto-create a new job in the project with optimization results
        new_job_id = _create_job_from_optimization(opt_id, params, result)

        _optimization_results[opt_id].update({
            "status": "completed",
            "progress": 100,
            "result": result,
            "best_score": result["final_lead"]["score"],
            "created_job_id": new_job_id,
        })

        logger.info(
            "Optimization %s completed: score %.3f -> %.3f (new job: %s)",
            opt_id,
            result["starting_molecule"]["score"],
            result["final_lead"]["score"],
            new_job_id,
        )

    except Exception as exc:
        error_msg = f"{type(exc).__name__}: {exc}"
        logger.error(
            "Optimization %s failed: %s\n%s",
            opt_id, error_msg, traceback.format_exc(),
        )
        _optimization_results[opt_id].update({
            "status": "failed",
            "error_message": error_msg,
        })


# ---------------------------------------------------------------------------
# BindX: Target Assessment endpoints
# ---------------------------------------------------------------------------

_assessment_pool = ThreadPoolExecutor(max_workers=2, thread_name_prefix="assessment")


@app.post("/api/target-assessment")
async def run_target_assessment(req: TargetAssessmentRequest) -> JSONResponse:
    """Run a full target assessment (5 scores + recommendation + optional AI agents).

    Returns the assessment result with composite score and GO/CAUTION/NO-GO.
    Cached by input_hash for idempotence.
    """
    import hashlib as _hashlib

    # Compute input hash for cache lookup
    input_hash = _hashlib.md5(
        json.dumps({
            "uniprot_id": req.uniprot_id,
            "disease_context": req.disease_context,
            "modality": req.modality,
        }, sort_keys=True).encode()
    ).hexdigest()

    # Check cache (skip if force_refresh or if agents were requested but previous cache has no agent data)
    cached = get_assessment_by_hash(input_hash)
    if cached and not req.force_refresh:
        result = json.loads(cached.assessment_json)
        agent_data = json.loads(cached.agent_responses_json) if cached.agent_responses_json else None

        # If agents requested but cached result has unavailable agents, re-run agents only
        needs_agent_rerun = (
            req.include_agents
            and (not agent_data or not agent_data.get("available", False))
        )
        if not needs_agent_rerun:
            logger.info("Returning cached assessment for %s (hash=%s)", req.uniprot_id, input_hash)
            return JSONResponse({
                "id": cached.id,
                "cached": True,
                **result,
                "agent_analysis": agent_data,
            })

        # Re-run agent with cached assessment scores
        logger.info("Cached assessment found but agent unavailable — re-running agent for %s", req.uniprot_id)
        try:
            from pipeline.agents.target_agent import TargetAssessmentAgent
            agent = TargetAssessmentAgent()
            agent_context = {
                "uniprot_id": req.uniprot_id,
                "disease_context": req.disease_context,
                "scores": result.get("scores", {}),
                "composite_score": result.get("composite_score"),
                "recommendation": result.get("recommendation"),
                "flags": result.get("flags", []),
            }
            loop = asyncio.get_event_loop()
            agent_data = await loop.run_in_executor(
                _assessment_pool,
                lambda: agent.query_sync(agent_context),
            )
            # Update DB with fresh agent data
            try:
                update_assessment_agent(cached.id, json.dumps(agent_data))
            except Exception:
                pass
        except Exception as agent_exc:
            logger.warning("Agent re-run failed: %s", agent_exc)
            agent_data = {"available": False, "fallback": str(agent_exc)}

        return JSONResponse({
            "id": cached.id,
            "cached": True,
            **result,
            "agent_analysis": agent_data,
        })

    # Get preview data from project if available
    preview_data = None
    if req.project_id:
        project = get_project(req.project_id)
        if project and project.target_preview_json:
            try:
                preview_data = json.loads(project.target_preview_json) if isinstance(
                    project.target_preview_json, str
                ) else project.target_preview_json
            except (json.JSONDecodeError, TypeError):
                pass

    # Run assessment in thread pool (involves blocking HTTP calls)
    loop = asyncio.get_event_loop()
    try:
        from pipeline.target_assessment import assess_target
        assessment = await loop.run_in_executor(
            _assessment_pool,
            lambda: assess_target(
                uniprot_id=req.uniprot_id,
                disease_context=req.disease_context,
                modality=req.modality,
                preview_data=preview_data,
                weights=req.weights,
            ),
        )
    except Exception as exc:
        logger.error("Target assessment failed for %s: %s", req.uniprot_id, exc)
        raise HTTPException(500, f"Assessment failed: {exc}")

    # Optionally run Agent 1 (Target Assessment Agent)
    agent_data = None
    if req.include_agents:
        try:
            from pipeline.agents.target_agent import TargetAssessmentAgent
            agent = TargetAssessmentAgent()
            agent_context = {
                "uniprot_id": req.uniprot_id,
                "disease_context": req.disease_context,
                "scores": assessment.get("scores", {}),
                "composite_score": assessment.get("composite_score"),
                "recommendation": assessment.get("recommendation"),
                "flags": assessment.get("flags", []),
            }
            agent_result = await loop.run_in_executor(
                _assessment_pool,
                lambda: agent.query_sync(agent_context),
            )
            agent_data = agent_result
        except Exception as agent_exc:
            logger.warning("Agent 1 failed: %s", agent_exc)
            agent_data = {"available": False, "fallback": str(agent_exc)}

    # Store in DB
    assessment_id = str(uuid4())
    try:
        create_assessment(
            assessment_id=assessment_id,
            uniprot_id=req.uniprot_id,
            assessment_json=json.dumps(assessment),
            project_id=req.project_id,
            disease_context=req.disease_context,
            agent_responses_json=json.dumps(agent_data) if agent_data else None,
            input_hash=input_hash,
        )
    except Exception as db_exc:
        logger.warning("Failed to store assessment: %s", db_exc)

    return JSONResponse({
        "id": assessment_id,
        "cached": False,
        **assessment,
        "agent_analysis": agent_data,
    })


@app.get("/api/target-assessment/{assessment_id}")
async def get_target_assessment(assessment_id: str) -> JSONResponse:
    """Retrieve a stored target assessment by ID."""
    row = get_assessment(assessment_id)
    if not row:
        raise HTTPException(404, "Assessment not found")

    result = json.loads(row.assessment_json) if row.assessment_json else {}
    agent_data = json.loads(row.agent_responses_json) if row.agent_responses_json else None

    return JSONResponse({
        "id": row.id,
        "uniprot_id": row.uniprot_id,
        "disease_context": row.disease_context,
        "created_at": row.created_at.isoformat() if row.created_at else None,
        **result,
        "agent_analysis": agent_data,
    })


# ---------------------------------------------------------------------------
# BindX: Agent ad-hoc query endpoint
# ---------------------------------------------------------------------------

AGENT_REGISTRY = {
    "target": "pipeline.agents.target_agent.TargetAssessmentAgent",
    "run_analysis": "pipeline.agents.run_analysis_agent.RunAnalysisAgent",
    "candidate": "pipeline.agents.candidate_agent.CandidateEvaluationAgent",
    "optimization": "pipeline.agents.optimization_agent.OptimizationStrategyAgent",
}


@app.post("/api/agent/{agent_name}/query")
async def query_agent(agent_name: str, req: AgentQuery) -> JSONResponse:
    """Query a specific AI agent with arbitrary context.

    Available agents: target, run_analysis, candidate, optimization.
    If project_id is provided, enriches context with project target info.
    """
    if agent_name not in AGENT_REGISTRY:
        raise HTTPException(
            404,
            f"Unknown agent '{agent_name}'. Available: {list(AGENT_REGISTRY.keys())}",
        )

    # Enrich context with project info if project_id provided
    enriched_context = dict(req.context)
    if req.project_id:
        try:
            from database import get_db
            from models import ProjectORM
            with get_db() as session:
                proj = session.query(ProjectORM).filter_by(id=req.project_id).first()
                if proj:
                    if proj.uniprot_id and "uniprot_id" not in enriched_context:
                        enriched_context["uniprot_id"] = proj.uniprot_id
                    # Extract gene name and target name from target_preview_json
                    if proj.target_preview_json:
                        try:
                            tp = json.loads(proj.target_preview_json) if isinstance(proj.target_preview_json, str) else proj.target_preview_json
                            if tp.get("protein_name") and "target_name" not in enriched_context:
                                enriched_context["target_name"] = tp["protein_name"]
                            # Try to get gene name from assessment or UniProt
                            assessment = tp.get("assessment", {})
                            if assessment.get("provenance", {}).get("gene_name"):
                                enriched_context["gene_name"] = assessment["provenance"]["gene_name"]
                        except (json.JSONDecodeError, TypeError):
                            pass
        except Exception as proj_exc:
            logger.debug("Could not enrich context from project %s: %s", req.project_id, proj_exc)

    module_path, class_name = AGENT_REGISTRY[agent_name].rsplit(".", 1)

    try:
        import importlib
        mod = importlib.import_module(module_path)
        agent_class = getattr(mod, class_name)
        agent = agent_class()

        loop = asyncio.get_running_loop()
        result = await loop.run_in_executor(
            _assessment_pool,
            lambda: agent.query_sync(enriched_context),
        )

        return JSONResponse(result)

    except Exception as exc:
        logger.error("Agent '%s' query failed: %s", agent_name, exc)
        return JSONResponse({
            "available": False,
            "agent_name": agent_name,
            "fallback": "Agent analysis failed. Please try again.",
        })


def _compute_score_stats(molecules: list) -> dict:
    """Compute score statistics for a list of molecules."""
    if not molecules:
        return {}
    affinities = [m.get("affinity", 0) for m in molecules if m.get("affinity") is not None]
    scores = [m.get("composite_score", 0) for m in molecules if m.get("composite_score") is not None]
    qeds = [m.get("qed", 0) for m in molecules if m.get("qed") is not None]
    return {
        "n_molecules": len(molecules),
        "affinity": {"min": min(affinities) if affinities else None, "max": max(affinities) if affinities else None,
                     "mean": sum(affinities) / len(affinities) if affinities else None},
        "composite_score": {"min": min(scores) if scores else None, "max": max(scores) if scores else None,
                           "mean": sum(scores) / len(scores) if scores else None},
        "qed": {"min": min(qeds) if qeds else None, "max": max(qeds) if qeds else None,
                "mean": sum(qeds) / len(qeds) if qeds else None},
    }


def _compute_hit_classes(molecules: list) -> dict:
    """Classify hits by score thresholds."""
    thresholds = {"excellent": 0.7, "good": 0.5, "moderate": 0.3}
    counts = {k: 0 for k in thresholds}
    counts["weak"] = 0
    for m in molecules:
        score = m.get("composite_score", 0)
        if score >= thresholds["excellent"]:
            counts["excellent"] += 1
        elif score >= thresholds["good"]:
            counts["good"] += 1
        elif score >= thresholds["moderate"]:
            counts["moderate"] += 1
        else:
            counts["weak"] += 1
    return counts


@app.post("/api/jobs/{job_id}/agent-analysis")
async def trigger_run_analysis(job_id: str) -> JSONResponse:
    """Trigger Agent 2 (Run Analysis) on a completed job's results.

    Builds context from the job's results.json and runs the Run Analysis Agent.
    """
    from database import get_db
    from models import JobORM

    with get_db() as session:
        job = session.query(JobORM).filter_by(id=job_id).first()
        if not job:
            raise HTTPException(404, f"Job {job_id} not found")
        if job.status != "completed":
            raise HTTPException(400, f"Job {job_id} is not completed (status: {job.status})")
        if not job.results_json:
            raise HTTPException(400, f"Job {job_id} has no results")

        molecules = json.loads(job.results_json)
        score_stats = _compute_score_stats(molecules)
        hit_classes = _compute_hit_classes(molecules)

        # Build context for Agent 2
        context = {
            "job_id": job_id,
            "mode": job.mode or "rapid",
            "docking_engine": job.docking_engine or "auto",
            "n_ligands_screened": len(molecules),
            "uniprot_id": job.uniprot_id,
            "protein_name": job.protein_name,
            "score_statistics": score_stats,
            "hit_classification": hit_classes,
            "top_molecules": [
                {
                    "name": m.get("name", ""),
                    "smiles": m.get("smiles", ""),
                    "affinity": m.get("affinity"),
                    "composite_score": m.get("composite_score"),
                    "qed": m.get("qed"),
                    "logp": m.get("logp"),
                    "mw": m.get("mw"),
                    "source": m.get("source"),
                    "cnn_score": m.get("cnn_score"),
                    "toxicity_level": m.get("toxicity_level"),
                    "admet_composite": m.get("admet", {}).get("composite_score") if isinstance(m.get("admet"), dict) else None,
                }
                for m in sorted(molecules, key=lambda x: x.get("composite_score", 0), reverse=True)[:10]
            ],
        }

    try:
        from pipeline.agents.run_analysis_agent import RunAnalysisAgent
        agent = RunAnalysisAgent()

        loop = asyncio.get_running_loop()
        result = await loop.run_in_executor(
            _assessment_pool,
            lambda: agent.query_sync(context),
        )

        return JSONResponse(result)

    except Exception as exc:
        logger.error("Run analysis agent failed for job %s: %s", job_id, exc)
        return JSONResponse({
            "available": False,
            "agent_name": "Run Analysis Agent",
            "fallback": f"Agent error: {exc}",
        })


# ---------------------------------------------------------------------------
# V9 — Time estimation endpoint
# ---------------------------------------------------------------------------

@app.post("/api/estimate-time")
async def estimate_time(req: Request) -> JSONResponse:
    """Estimate pipeline runtime based on run configuration parameters.

    Accepts JSON body with keys matching pipeline params:
    n_ligands, mode, docking_engine, enable_generation, etc.
    """
    try:
        body = await req.json()
    except Exception:
        body = {}

    from pipeline.time_estimation import estimate_pipeline_time, quick_estimate

    n_ligands = int(body.get("n_ligands", 50))
    mode = body.get("mode", "rapid")
    docking_engine = body.get("docking_engine", "auto")

    # Full estimate
    full = estimate_pipeline_time(
        n_ligands=n_ligands,
        mode=mode,
        docking_engine=docking_engine,
        enable_generation=body.get("enable_generation", mode != "rapid"),
        n_generated_molecules=int(body.get("n_generated_molecules", 100)),
        enable_retrosynthesis=body.get("enable_retrosynthesis", mode != "rapid"),
        enable_diffdock=body.get("enable_diffdock", False),
        structure_source=body.get("structure_source", "alphafold"),
        use_pubchem=body.get("use_pubchem", True),
        use_enamine=body.get("use_enamine", n_ligands < 100),
        use_chembl=body.get("use_chembl", True),
        include_agents=body.get("include_agents", False),
        protein_length=int(body.get("protein_length", 400)),
    )

    # Quick summary
    quick = quick_estimate(n_ligands=n_ligands, mode=mode, docking_engine=docking_engine)

    return JSONResponse({
        "estimate": full,
        "quick": quick,
    })


# ---------------------------------------------------------------------------
# V9 — Chemical database info endpoint
# ---------------------------------------------------------------------------

@app.get("/api/databases")
async def list_databases() -> JSONResponse:
    """List available chemical databases and their status."""
    return JSONResponse({
        "databases": [
            {
                "id": "chembl",
                "name": "ChEMBL",
                "description": "Bioactive compounds with drug-like properties from EBI",
                "type": "bioactivity",
                "url": "https://www.ebi.ac.uk/chembl/",
                "status": "active",
                "avg_query_time_s": 10,
            },
            {
                "id": "pubchem",
                "name": "PubChem",
                "description": "NCBI chemical compound database with bioassay data",
                "type": "bioactivity",
                "url": "https://pubchem.ncbi.nlm.nih.gov/",
                "status": "active",
                "avg_query_time_s": 12,
            },
            {
                "id": "zinc",
                "name": "ZINC20",
                "description": "Drug-like molecules for virtual screening",
                "type": "drug-like",
                "url": "https://zinc20.docking.org/",
                "status": "active",
                "avg_query_time_s": 1,
            },
            {
                "id": "enamine_real",
                "name": "Enamine REAL",
                "description": "Make-on-demand virtual library (37B+ compounds)",
                "type": "virtual_library",
                "url": "https://enamine.net/compound-collections/real-compounds",
                "status": "active",
                "avg_query_time_s": 2,
            },
        ],
    })


# ---------------------------------------------------------------------------
# V9 — Scientific references & expected results
# ---------------------------------------------------------------------------

@app.get("/api/references")
async def get_references(category: Optional[str] = None) -> JSONResponse:
    """Return scientific references for all tools/methods used in DockIt."""
    from pipeline.references import get_references, get_categories
    refs = get_references(category)
    cats = get_categories()
    return JSONResponse({"references": refs, "categories": cats})


@app.get("/api/expected-results")
async def get_expected_results() -> JSONResponse:
    """Return expected results and benchmarks for each pipeline mode."""
    from pipeline.references import get_expected_results
    return JSONResponse(get_expected_results())


# ---------------------------------------------------------------------------
# V9 — Pipeline analytics / monitoring
# ---------------------------------------------------------------------------

@app.get("/api/analytics")
async def get_analytics() -> JSONResponse:
    """Return pipeline analytics: run counts, success rates, avg times by mode.

    Useful for monitoring dashboards (n8n, Grafana) and performance tracking.
    """
    from database import get_db
    from models import JobORM

    with get_db() as db:
        all_jobs = db.query(JobORM).all()
        total = len(all_jobs)
        completed = sum(1 for j in all_jobs if j.status == "completed")
        failed = sum(1 for j in all_jobs if j.status == "failed")
        running = sum(1 for j in all_jobs if j.status == "running")

        # Mode breakdown
        mode_counts = {}
        for j in all_jobs:
            m = j.mode or "rapid"
            mode_counts[m] = mode_counts.get(m, 0) + 1

        # Success rate
        success_rate = (completed / total * 100) if total > 0 else 0

        # Recent failures (last 5)
        recent_failures = [
            {
                "job_id": j.id,
                "error": j.error_message or "Unknown",
                "mode": j.mode,
                "created_at": j.created_at.isoformat() if j.created_at else None,
            }
            for j in sorted(all_jobs, key=lambda x: x.created_at or "", reverse=True)
            if j.status == "failed"
        ][:5]

        return JSONResponse({
            "total_jobs": total,
            "completed": completed,
            "failed": failed,
            "running": running,
            "success_rate_pct": round(success_rate, 1),
            "mode_breakdown": mode_counts,
            "recent_failures": recent_failures,
        })
