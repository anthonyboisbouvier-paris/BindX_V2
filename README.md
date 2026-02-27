# BindX V2 — In Silico Drug Discovery Platform

> Architecture V9 : Project > Campaign > Phase > Run

## Stack

- **Backend** : FastAPI + Celery + Redis + Supabase (PostgreSQL)
- **Frontend** : React 18 + Vite + Tailwind + Recharts
- **Auth** : Supabase Auth (JWT, RLS)
- **Docking GPU** : GNINA via RunPod serverless

## Quick Start

```bash
cp .env.example .env          # Edit with Supabase credentials
docker-compose up --build
```

- Frontend : http://localhost:3000
- Backend API : http://localhost:8000/docs

## Project Structure

```
BindX_V2/
├── backend/
│   ├── main.py                 # FastAPI app (V8 endpoints + V9 router mount)
│   ├── models_v9.py            # 9 SQLAlchemy ORM models (V9)
│   ├── database_v9.py          # Async Supabase connection
│   ├── auth_v9.py              # Supabase JWT verification
│   ├── schemas_v9.py           # Pydantic request/response schemas
│   ├── tasks_v9.py             # Celery tasks (import/calculation/generation)
│   ├── celery_app.py           # Celery + Redis config
│   ├── routers/v9/             # 6 API routers (health, projects, campaigns, phases, runs, molecules)
│   └── pipeline/               # 35 computation modules (docking, scoring, ADMET, agents...)
├── frontend/src/
│   ├── pages/                  # ProjectListPage, ProjectHome, PhaseDashboard, Login, Register
│   ├── components/             # SidebarLayout, MoleculeTable, RunCreator, RunProgress, ErrorBoundary...
│   ├── contexts/               # AuthContext (Supabase), WorkspaceContext, ToastContext
│   ├── api.js                  # V9 API client (runs, molecules, bookmarks)
│   ├── mock/data.js            # Mock data for dev
│   └── lib/supabase.js         # Supabase client
├── infrastructure/
│   └── supabase/migrations/    # v9_schema.sql, v9_rls.sql
├── docs/
│   └── CDC_V9_FINAL.md         # Cahier des charges complet
├── docker-compose.yml          # redis + backend + celery_worker + frontend
└── CLAUDE.md                   # Context for Claude Code AI assistant
```

## Architecture V9

### Data Model
- **Project** : 1 cible therapeutique (UniProt/PDB/FASTA)
- **Campaign** : 1 pocket selectionnee + regles de scoring
- **Phase** : Hit Discovery > Hit-to-Lead > Lead Optimization
- **Run** : unite atomique de calcul

### Run Types
| Type | Description |
|------|-------------|
| **import** | SDF/SMILES ingestion + pre-filtres |
| **calculation** | Docking, ADMET, scoring, enrichment, clustering, off-target, confidence, retrosynthesis, safety |
| **generation** | De novo (REINVENT4) |

### Phase Dashboard
- 1 molecule = 1 row (deduplicated by canonical SMILES)
- Runs add columns progressively
- Bookmark + Freeze workflow between phases

## Specs

Full specification : [docs/CDC_V9_FINAL.md](docs/CDC_V9_FINAL.md)
