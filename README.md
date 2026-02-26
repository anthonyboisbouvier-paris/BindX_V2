# 🧬 BindX V2 - Drug Discovery Platform

> **Version 9** - Architecture modulaire Project → Campaign → Phase → Run

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Node.js 18+](https://img.shields.io/badge/Node.js-18+-green.svg)](https://nodejs.org/)
[![Docker](https://img.shields.io/badge/Docker-Required-blue.svg)](https://www.docker.com/)

## 🚀 Quick Start (Ubuntu/WSL)

```bash
git clone https://github.com/anthonyboisbouvier-paris/BindX_V2.git
cd BindX_V2
./scripts/dev.sh
```

**Prérequis :** Docker, Docker Compose, Node.js 18+, Python 3.10+

## 📚 Documentation

- 📋 [CDC V9 Final](docs/CDC_V9_FINAL.md) - Cahier des charges complet  
- 🏗️ [Architecture](docs/ARCHITECTURE.md) - Vue technique  
- 🔄 [Migration V8→V9](docs/MIGRATION_V8_TO_V9.md) - Guide migration  
- 👨‍💻 [Development](docs/DEVELOPMENT.md) - Setup développement  
- 🤖 [Claude Code](tools/claude/) - Workflows IA  

## 🎯 Différence avec DockIt V8

| Aspect | DockIt V8 | BindX V9 |
|--------|-----------|----------|
| Architecture | Job monolithique | Project → Campaign → Phase → Run |
| Base de données | SQLite | Supabase (PostgreSQL) |
| Runs | 1 job = tout | Runs modulaires (import/calcul/génération) |
| Interface | Job-centrique | Phase-dashboard centrique |
| Cible | 1 par job | Target Preview + poches multiples |
| Calculs | Bundle fixe | Multi-sélection (docking+ADMET+scoring...) |

## 🧪 Features V9

### Architecture modulaire
- **Project** : cible thérapeutique + target preview automatique
- **Campaign** : stratégie d'exploration (pocket sélectionnée)  
- **Phase** : Hit Discovery → Hit-to-Lead → Lead Optimization
- **Run** : unité atomique (import/calcul/génération)

### Types de runs
- **Import** : SDF/SMILES, bases connectées, pré-filtres Lipinski/PAINS
- **Calcul** : multi-sélection modulaire
  - 🎯 Docking (Vina, GNINA CPU/GPU, DiffDock)
  - 💊 ADMET (logP, solubility, BBB, hERG, CYP, bioavailability)
  - 📊 Scoring composite pondéré
  - 🔗 Enrichment (ProLIF interactions, clustering)
  - 🎨 Clustering (diversité, scaffolds, Tanimoto)
  - ⚠️ Off-target sélectivité  
  - ✅ Confidence (PAINS, applicability domain)
  - ⚗️ Retrosynthesis (faisabilité, coût, réactifs)
  - 🛡️ Safety (hERG, AMES, hepatotox, carcinogenicity)
- **Génération** : de novo REINVENT4 (mode batch + molécule individuelle)

### Dashboard cumulatif par phase
- **1 molécule = 1 ligne** dédupliquée par SMILES canonique
- **Runs ajoutent des colonnes** progressivement
- **Popups détaillées** : safety, synthesis, confidence
- **Viewer 3D intégré** : drawer 40% + plein écran, 3Dmol.js
- **Bookmarks & freeze** pour pipeline phase → phase

### Agent IA campagne
- **Cross-phases analysis** : attrition, enrichment factor
- **Recommandations intelligentes** : prochains runs, paramètres
- **Safety alerts** : red flags ADMET, scaffold analysis
- **Chat contextuel** avec accès toutes les données

## 🛠️ Tech Stack

- **Frontend** : React 18, Vite, TailwindCSS, 3Dmol.js, Recharts
- **Backend** : FastAPI, SQLAlchemy, Pydantic V2, Celery, Redis  
- **Database** : Supabase (PostgreSQL + Auth + Storage + RLS)
- **Containers** : Docker, Docker Compose, multi-stage builds
- **Pipeline** : AutoDock Vina, GNINA, GNINA GPU (RunPod), DiffDock
- **Cheminformatics** : RDKit, OpenEye, fpocket, ProLIF, MDAnalysis
- **ML** : REINVENT4, scikit-learn, PyTorch, Transformers
- **Deployment** : Nginx, SSL auto-renewal, domain routing

## 📦 Project Structure

```
BindX_V2/
├── docs/                    # 📚 Documentation complète
├── backend/                 # 🐍 API FastAPI + pipeline
├── frontend/                # ⚛️ Interface React
├── infrastructure/          # 🏗️ Docker, Supabase, deploy
├── tools/claude/           # 🤖 Claude Code workflows  
├── scripts/                # 🚀 Automation
└── data/                   # 📊 Structures, cache (gitignored)
```

## 🚀 Development Quickstart

### 1. Environment Setup
```bash
# Clone & setup
git clone https://github.com/anthonyboisbouvier-paris/BindX_V2.git
cd BindX_V2

# Configure environment  
cp .env.example .env
# Edit .env with your Supabase credentials
```

### 2. Docker Development
```bash
# Start all services
./scripts/dev.sh

# Or manually
docker-compose up -d
```

### 3. Access Points
- 🌐 **Frontend** : http://localhost:3000
- 🐍 **Backend** : http://localhost:8000  
- 📚 **API Docs** : http://localhost:8000/docs
- 🗄️ **Database** : Supabase Dashboard

### 4. Claude Code Integration
```bash
# AI-powered development
cd ~/BindX_V2
claude

# Available shortcuts (see tools/claude/shortcuts.sh)
bx-dev      # Start development
bx-logs     # View container logs  
bx-test     # Run tests
bx-reset    # Restart containers
```

## 🧬 Scientific Pipeline

### Supported Calculations
- **Structure Processing** : PDB cleanup, pocket detection (fpocket), validation
- **Molecular Docking** : Vina (CPU), GNINA (CPU/GPU), DiffDock (ML-based)
- **ADMET Prediction** : Solubility, permeability, metabolism, toxicity
- **Drug-likeness** : Lipinski Ro5, QED, PAINS filtering
- **De Novo Generation** : REINVENT4 goal-directed optimization
- **Interaction Analysis** : ProLIF protein-ligand fingerprints
- **Selectivity Screening** : Off-target binding prediction
- **Synthetic Accessibility** : Retrosynthetic route planning

### Validation & Quality
- **Confidence Scoring** : Applicability domain, model uncertainty
- **Safety Profiling** : hERG, AMES, hepatotoxicity, carcinogenicity  
- **PAINS Detection** : Pan-assay interference compounds
- **Structural Alerts** : Reactive/toxic substructures

## 📊 Data Sources

- **ChEMBL** : Bioactivity data, known targets
- **PubChem** : Chemical structures, properties
- **UniProt** : Protein sequences, annotations
- **PDB** : 3D structures (experimental)
- **AlphaFold** : 3D structures (predicted)
- **DrugBank** : Approved drugs, targets

## 🔗 Integration

### Docking Engines
- **AutoDock Vina** : Fast, reliable CPU docking
- **GNINA** : CNN-enhanced scoring (CPU/GPU)
- **DiffDock** : Diffusion model for complex poses

### Cloud Computing  
- **RunPod Serverless** : GPU acceleration (10x speedup)
- **Auto-fallback** : CPU → GPU → Mock for development

### Authentication
- **Supabase Auth** : Email/password, magic links
- **Row Level Security** : Multi-tenant data isolation
- **SSO Ready** : SAML/OAuth2 (V10)

## 🔄 Migration from DockIt V8

BindX V9 is a complete rewrite with **no automatic migration**. See [MIGRATION.md](docs/MIGRATION_V8_TO_V9.md) for:

- Data export strategies
- Workflow mapping V8→V9  
- Feature comparison matrix
- Step-by-step migration guide

**Legacy support** : V8 projects accessible in read-only mode at `/legacy`.

## 🤝 Contributing

1. Fork the repository
2. Create feature branch (`git checkout -b feature/amazing-feature`)
3. Commit changes (`git commit -m 'Add amazing feature'`)
4. Push to branch (`git push origin feature/amazing-feature`)
5. Open Pull Request

See [DEVELOPMENT.md](docs/DEVELOPMENT.md) for detailed contributor guide.

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **OpenEye Scientific** : OMEGA conformer generation
- **Schrödinger** : Molecular modeling inspiration
- **ChEMBL Team** : Bioactivity database
- **RDKit Community** : Cheminformatics toolkit
- **AutoDock Team** : Docking algorithms

## 📞 Support

- 🐛 **Issues** : GitHub Issues tracker
- 💬 **Discussions** : GitHub Discussions  
- 📖 **Documentation** : `/docs` folder
- 🤖 **AI Assistant** : Claude Code integration

---

**🔬 Built for computational chemists, by computational chemists.**

**From molecular hypothesis to lead compound in record time.**
