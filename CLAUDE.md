# BindX V2 - Plateforme de Drug Discovery In Silico

## MISSION
Tu es un developpeur senior autonome. Tu developpes BindX V9, une refonte structurelle.
Specs completes : docs/CDC_V9_FINAL.md

## ARCHITECTURE V9

### Modele conceptuel
Project (1 cible) > Campaign (cible+pocket+regles) > Phase (A/B/C) > Run (unite calcul)

### Stack
- Backend : FastAPI + Celery + Redis + Supabase (PostgreSQL)
- Frontend : React + Vite + Tailwind
- Auth : Supabase Auth
- Storage : Supabase Storage
- Docking GPU : GNINA via RunPod serverless

### Principes cles
- Dashboard Phase = vue centrale (1 mol = 1 ligne, runs ajoutent colonnes)
- 1 run = 1 type calcul (import, docking, admet, scoring, enrichment, generation, clustering)
- Input run = lignes cochees dans le dashboard
- Bookmark = tag modifiable, Freeze = verrouille pour phase suivante (reversible avec warning)
- 1 run simultane par user V9
- Legacy V8 = lecture seule separee

## STRUCTURE PROJET
```
BindX_V2/
├── docs/                    # Specs, CDC, methodes
├── backend/                 # FastAPI + pipeline + Celery
│   ├── main.py             # API endpoints
│   ├── models.py           # Pydantic + SQLAlchemy
│   ├── database.py         # DB connection
│   ├── auth.py             # JWT auth
│   ├── tasks.py            # Celery tasks
│   ├── celery_app.py       # Celery config
│   ├── pipeline/           # 27 modules (docking, scoring, ADMET, agents...)
│   ├── tests/              # Unit tests
│   └── benchmarks/         # Kinase benchmark suite
├── frontend/               # React + Vite
│   └── src/
│       ├── pages/          # ProjectListPage, ProjectHome, PhaseDashboard, PharmacoDB
│       ├── components/     # 30+ composants (V9 + V8 reutilises)
│       ├── contexts/       # WorkspaceContext, ToastContext, AuthContext
│       ├── mock/           # Mock data (dev)
│       └── lib/            # Supabase client
├── infrastructure/         # Docker, Supabase, deploy
│   ├── docker/nginx/
│   └── supabase/migrations/
├── tools/
│   ├── claude/             # Claude Code workflows
│   └── pharmaco_db/        # MegaPharmaDB scripts
├── scripts/                # dev.sh, automation
└── data/                   # gitignored — structures, cache
```

## REGLES AUTONOMIE
- Ne jamais demander confirmation sauf suppression donnees ou push prod
- Bloque > 3 tentatives : simplifier (mock/stub) + noter TODO
- Tester unitairement chaque module
- docker-compose up --build pour integration

## EPICS V9 (7 semaines)
1. Fondations Supabase (S1-S2)
2. Modele Project/Campaign/Phase (S2-S3)
3. Runs par type (S3-S4)
4. Dashboard Phase (S4-S5)
5. Navigation UX (S5-S6)
6. Polish + Legacy (S6-S7)

## DEMARRAGE RAPIDE
```
docker-compose up --build
# Frontend: http://localhost:3000
# Backend: http://localhost:8000/docs
```
