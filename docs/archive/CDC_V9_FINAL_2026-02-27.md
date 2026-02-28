# BindX V9 — Cahier des Charges Produit

## Refonte structurelle : Project > Campaign > Phase > Run

**Version** : 9.0 — CDC Final  
**Date** : 25 février 2026  
**Auteur** : Anthony Boisbouvier + Claude Opus 4.6

---

## 1. Vision produit

BindX est une plateforme open-source de drug discovery in silico. La V9 transforme le pipeline monolithique actuel (1 job = tout) en une architecture modulaire pensée pour le workflow réel d'un chimiste computationnel :

- **Project** : une cible thérapeutique principale (champ libre)
- **Campaign** : une stratégie d'exploration (cible + pocket + objectif)
- **Phase** : une étape du funnel (Hit Discovery > Hit-to-Lead > Lead Optimization)
- **Run** : une unité de calcul atomique qui enrichit le dashboard de la phase (ou un import de composés)

**Principe fondamental** : le dashboard de la phase est la vue centrale. Chaque run ajoute des colonnes (propriétés calculées) aux molécules, ou des lignes lorsque c'est un import de nouveaux composés. L'utilisateur travaille toujours dans un seul tableau qui s'enrichit progressivement.

---

## 2. Personas

### Chimiste computationnel (utilisateur principal)
Expert drug discovery. Lance des campagnes de screening virtuel, sélectionne des hits, optimise des leads. Veut un outil rapide, pas de friction, des résultats exploitables immédiatement.

### Bioinformaticien
Prépare les structures protéiques, configure les poches de binding. Veut du contrôle sur les paramètres techniques.

### Chef de projet pharma
Suit l'avancement des campagnes, compare les phases, exporte les résultats pour les comités.

---

## 3. Modèle conceptuel

### 3.1 PROJECT

Conteneur organisationnel. 1 projet = 1 cible thérapeutique principale.

| Champ | Description |
|-------|-------------|
| name | Nom du projet (ex: EGFR Inhibitors) |
| description | Contexte scientifique |
| status | active, completed, archived |

#### Rubrique Target (création du projet)

Au moment de la création du projet, l'utilisateur renseigne sa cible via l'une des méthodes suivantes :

- **UniProt ID** (ex: P00533) — déclenche la récupération automatique de la structure
- **Séquence FASTA** — collée directement, validée (commence par M, AA valides, min 50 résidus)
- **PDB ID** — structure expérimentale directe

Le backend effectue automatiquement le **Target Preview** :

| Information | Source |
|-------------|--------|
| Structure source | PDB expérimental, AlphaFold, ESMFold (détection automatique avec fallback) |
| Résolution, méthode | RCSB PDB |
| Ligand co-cristallisé | PDB |
| Poches détectées | fpocket (classées par druggability) |
| Actifs connus | ChEMBL (nombre + IC50 médiane) |
| Prédiction désordre (IDR) | IUPred / AlphaFold pLDDT |

Ces informations sont affichées dans une **TargetPreviewCard** persistante dans le projet et réutilisées lors de la création de campagnes.

V10 : support multi-cibles par projet (polypharmacologie).

### 3.2 CAMPAIGN

Programme scientifique concret. Définit la stratégie d'exploration.

| Champ | Description |
|-------|-------------|
| project_id | FK vers Project |
| name | Nom (ex: ATP pocket - ChEMBL screening) |
| target_config | Cible + pocket sélectionnée (parmi celles détectées dans le Target Preview) |

- 1 campagne par défaut à la création du projet
- UI pour en ajouter (pas prioritaire V9)
- La campagne choisit UNE cible + UNE pocket du projet

### 3.3 PHASE

Étape du funnel de drug discovery. L'utilisateur crée chaque phase quand il en a besoin.

| Type | Objectif | Input typique |
|------|----------|---------------|
| A — Hit Discovery | Trouver des hits dans une librairie | Librairie externe (SDF/SMILES) ou librairie standard connectée à BindX |
| B — Hit-to-Lead | Étudier les hits, identifier des patterns | Molécules sélectionnées de Phase A |
| C — Lead Optimization | Affiner les leads via génération de novo, ADMET, sélectivité | Molécules sélectionnées de Phase B |

Chaque phase possède :
- Un **dashboard unique** : tableau cumulatif de toutes les molécules + propriétés calculées par les runs
- Un **état** : active, frozen, completed
- **Colonnes affichées** : toutes les colonnes avec résultats sont affichées par défaut ; l'utilisateur peut masquer celles qu'il ne souhaite pas voir

#### Freeze / Unfreeze
- **Bookmark** : tag interne, modifiable à tout moment
- **Freeze** : verrouille la sélection bookmarkée, disponible comme input phase suivante
- **Unfreeze** : possible avec warning si la phase suivante a déjà des runs

#### Experimental Results (placeholder V9)
Rubrique visible dans la sidebar sous les phases. Servira au fine-tuning futur : l'utilisateur pourra importer des résultats expérimentaux (IC50 mesurées, etc.) qui seront intégrés via des runs dans les phases quand disponible.

V10 : phases multiples du même type, intégration experimental results.

### 3.4 RUN

Unité de calcul atomique. 3 types de run.

#### Types de run

**Run Import**

| Source | Description |
|--------|-------------|
| Fichier SDF/SMILES | Upload externe |
| Sélection interne | Molécules sélectionnées de la phase précédente (bookmarked + frozen) |
| Base connectée | Librairie standard connectée à BindX |

Options de pré-filtres disponibles à l'import (Lipinski, MW range, etc.).

Colonnes ajoutées : SMILES, canonical_smiles, name, source.

**Run Calcul**

L'utilisateur sélectionne un ou plusieurs types de calcul dans un seul run via des checkboxes :

| Calcul | Description | Colonnes ajoutées |
|--------|-------------|-------------------|
| Docking | Docking moléculaire | docking_score, CNNscore, CNNaffinity, poses |
| ADMET | Propriétés ADMET complètes | logP, solubility, BBB, hERG_inhibition, metabolic_stability, oral_bioavailability, plasma_protein_binding, CYP inhibitions |
| Scoring | Score composite pondéré | composite_score |
| Enrichment | ProLIF interactions, clustering | interactions, interaction_count, cluster_id, scaffold |
| Clustering | Diversité, scaffolds, Tanimoto | cluster_id, scaffold_smiles, tanimoto_to_centroid |
| Off-target | Sélectivité off-target | selectivity_score, off_target_hits, selectivity_ratio |
| Confidence | PAINS, applicability domain, convergence | confidence_score, pains_alert, applicability_domain, confidence_flags |
| Retrosynthesis | Faisabilité de synthèse | n_synth_steps, synth_confidence, synth_cost_estimate, reagents_available |
| Safety | Profil de sécurité complet | herg_risk, ames_mutagenicity, hepatotoxicity, skin_sensitization, carcinogenicity, safety_color_code |

Chaque calcul coché ajoute ses colonnes au dashboard. Le backend exécute les calculs séquentiellement dans l'ordre logique (ex: docking avant scoring).

**Configuration du docking** (quand coché) :
- Moteur : Vina (CPU classique), GNINA CPU (CNN-scored), GNINA GPU (RunPod serverless), DiffDock (deep learning, meilleur pour poses flexibles)
- Box size : automatique depuis la pocket ou personnalisée
- Exhaustiveness : paramétrable

**Configuration du scoring** (quand coché) :
- Éditeur visuel des poids (ScoringWeightsEditor)
- Poids par défaut modifiables par run

**Run Génération**

Spécifique aux Phases B et C. Ouvre une **interface dédiée** avec les composés sélectionnés :

**Mode batch** : sélectionner un groupe de molécules et configurer :
- Nombre d'itérations
- Nombre de molécules générées par itération
- Critères génériques (QED min, Lipinski, PAINS filter)
- Tous les critères existants de la version actuelle (REINVENT4 ou mock RDKit scaffold hopping)

**Mode molécule** : sélectionner une molécule individuelle et configurer :
- Critères au niveau des R-groupes spécifiques (comme DockIt V8)
- Contrôle fin des modifications structurelles

Output :
- N molécules générées avec TOUS les calculs déjà dans le dashboard cochés par défaut (ex: docking + ADMET + scoring si déjà calculés)
- Métadonnées : `generation_level` (0=originale, 1=gen1, 2=gen2...), `parent_smiles` pour tracer la lignée
- Tag **"IA generated"** visible dans le dashboard
- Lien parent → enfant intelligent dans le dashboard (ex: icône arbre, expand/collapse, ou colonne "parent" cliquable)

Chaque itération ajoute des lignes dans le dashboard de la phase.

#### Input d'un run
- **Import** : fichier SDF/SMILES externe OU sélection de molécules internes OU base connectée, avec option de pré-filtres
- **Calcul** : les lignes cochées dans le dashboard de la phase (sélection manuelle)
- **Génération** : les lignes cochées dans le dashboard
- Pas de multi-source : un run prend soit un import externe, soit des cochées du dashboard
- **Select all filtered** : bouton sélectionner toutes les molécules filtrées (pas juste les visibles)

#### Cycle de vie d'un run
`created` → `queued` → `running` → `completed` / `failed`

- Pas de suppression : archive seulement
- Plusieurs runs simultanés si l'implémentation le permet sans risque, sinon 1 run simultané par utilisateur avec file d'attente visible
- **Estimation du temps** affichée avant lancement (basée sur nb molécules × type de calcul)
- **Progress tracking** : barre de progression + étape courante + détails (ex: "Docking 142/500 — GNINA GPU")

#### Cache
`Hash(smiles + run_type + params)` → résultat. Réutilisation sans relancer si déjà calculé.

---

## 4. Dashboard Phase — La vue centrale

Le dashboard est LE lieu de travail de l'utilisateur. Tableau enrichi progressivement par les runs.

### 4.1 Principes
- 1 molécule = 1 ligne (dédupliquée par SMILES canonique)
- Chaque run calcul ajoute des colonnes, pas des lignes (sauf import et génération)
- Toutes les colonnes avec résultats sont affichées par défaut
- L'utilisateur peut masquer/réafficher n'importe quelle colonne

### 4.2 Fonctionnalités du dashboard
- **Tri** : toutes colonnes, multi-sort
- **Filtres** : range sliders numériques, texte catégoriel
- **Sélection** : checkbox par ligne + select all filtered + select bookmarked
- **Bookmark** : toggle par molécule
- **Actions sur sélection** : New Run, Export, Bookmark All
- **Vue 3D** : clic molécule → Viewer3D drawer latéral
- **Safety popup** : clic sur safety_color_code → popup détaillée du profil de sécurité (hERG, AMES, hepatotox, etc.) avec code couleur vert/jaune/rouge
- **Retrosynthesis popup** : clic sur synth_confidence → popup avec arbre rétrosynthétique, coût estimé, disponibilité réactifs (SynthesisTree)
- **Confidence popup** : clic sur confidence_score → popup avec breakdown des flags (PAINS, AD, convergence)

### 4.3 Viewer 3D
- Drawer latéral 40%, table compressée à 60%
- Toggle plein écran (overlay)
- Réutilise Viewer3D existant (3Dmol.js)
- Affiche : protéine + ligand posé + interactions ProLIF + surface pocket
- **Modes de coloration** : standard, hydrophobicité (échelle Kyte-Doolittle), B-factor, secondaire
- **Toggles** : ligand on/off, surface on/off, interactions on/off, pocket résidus on/off, labels on/off
- **Toutes les options de la version actuelle conservées**

### 4.4 InfoTips (pédagogie)
Chaque concept, colonne, paramètre et action dispose d'un **point info (ℹ️)** cliquable qui affiche une explication pédagogique contextuelle. Couvre :
- Que signifie chaque colonne (ex: "CNNscore = score de confiance du réseau de neurones GNINA, de 0 à 1")
- Pourquoi un paramètre est important
- Comment interpréter un résultat
- Bonnes pratiques (ex: "Un logP > 5 réduit généralement la biodisponibilité orale")

Réutilise le composant `InfoTip` existant. Les tips sont stockés en JSON et facilement éditables.

---

## 5. Scoring composite

### Run scoring
Un run scoring calcule un `composite_score` basé sur les colonnes déjà disponibles dans le dashboard.

Configuration par run via l'éditeur visuel `ScoringWeightsEditor` :
- Poids ajustables par glissière pour chaque colonne numérique
- Normalisation automatique (min-max par colonne)
- Prévisualisation du classement avant lancement

Pas de poids par défaut au niveau campagne — chaque run scoring est explicitement configuré.

---

## 6. Agents IA

Opèrent au niveau campagne. Accès lecture : toutes phases, tous runs, toutes molécules.

Capacités :
- Analyser tendances cross-phases (attrition, enrichment factor)
- Recommander prochains runs
- Identifier scaffolds prometteurs
- Suggérer paramètres scoring
- Alerter red flags ADMET / safety
- Répondre aux questions de l'utilisateur sur ses données (chat)
- Citer des faits vérifiables (ChEMBL, propriétés calculées)

Interface : panel chat latéral (CampaignAgentPanel), accessible depuis n'importe quelle phase de la campagne.

---

## 7. Optimisation multi-objectifs (Pareto)

Niveau phase : toutes les molécules de la phase.

- Front de Pareto sur N objectifs sélectionnables
- Visualisations : scatter plot, parallel coordinates
- Tag "Pareto-optimal" sur les molécules identifiées
- Conserve le composant ParetoFront + OptimizationView existants, scope = phase

---

## 8. Export

Scope : sélection cochée.

| Format | Contenu |
|--------|---------|
| CSV | Colonnes visibles, molécules cochées |
| SDF | Structures 3D + propriétés |
| PDF | Rapport synthétique campagne (réutilise le moteur de rapport existant) |

Raccourcis : Export bookmarked, Export Pareto.

Le rapport PDF inclut : résumé campagne, statistiques par phase, top molécules, graphiques Pareto, profils ADMET/safety des leads. Réutilise `ReportPreview` + `report.py`.

---

## 9. Navigation & UX

### Sidebar gauche

```
Project: EGFR Inhibitors
  Target: P00533 — EGFR (PDB: 1M17, AlphaFold)
  Campaign: ATP Pocket ChEMBL
    Phase A — Hit Discovery [active]
    Phase B — Hit-to-Lead [frozen]
    Phase C — Lead Optimization [locked]
  Experimental Results [placeholder]
```

### Flow principal exemple
1. Créer projet (nom + cible : UniProt ID, séquence FASTA, ou PDB ID)
2. Target Preview automatique (structure, poches, ChEMBL actifs)
3. Campagne par défaut créée (choisir pocket parmi celles détectées)
4. Créer Phase A
5. Run import : charger librairie SDF ou source connectée
6. Run calcul : cocher docking + ADMET → lance les deux
7. Run calcul : scoring composite
8. Dashboard : trier, filtrer, bookmarker hits
9. Freeze Phase A → input Phase B
10. Créer Phase B, run génération sur hits

### RunCreator (modal)
Remplace InputForm. Modal depuis dashboard :
- **Étape 1** : Choix type de run (Import / Calcul / Génération)
- **Étape 2** : Configuration
  - Import : source + pré-filtres
  - Calcul : checkboxes des calculs + paramètres par calcul coché (moteur docking, poids scoring, etc.)
  - Génération : mode batch ou molécule, paramètres
- **Étape 3** : Confirmation (nb molécules, estimation temps, calculs sélectionnés)
- Lancement → bandeau "Run in progress" + progress bar + étape courante

### File d'attente
Run en cours → message clair + run queued grisé dans historique.

### Notifications
Email optionnel à la complétion des runs longs (> 10 min estimées). L'utilisateur renseigne son email dans les settings du projet.

---

## 10. Infrastructure

### 10.1 Base de données — Supabase Cloud (PostgreSQL)

```sql
CREATE TABLE projects (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  user_id UUID REFERENCES auth.users(id),
  name TEXT NOT NULL,
  description TEXT,
  -- Target info
  target_input_type TEXT CHECK (target_input_type IN ('uniprot', 'sequence', 'pdb')),
  target_input_value TEXT NOT NULL,
  target_name TEXT,
  target_pdb_id TEXT,
  target_preview JSONB,
  -- Structure source info
  structure_source TEXT, -- 'experimental', 'alphafold', 'esmfold'
  structure_resolution FLOAT,
  structure_method TEXT,
  cocrystal_ligand TEXT,
  -- Pockets détectées
  pockets_detected JSONB, -- [{rank, score, residues, volume, druggability, idr_warning}]
  -- ChEMBL info
  chembl_actives_count INTEGER,
  chembl_median_ic50 FLOAT,
  --
  status TEXT DEFAULT 'active',
  notification_email TEXT,
  created_at TIMESTAMPTZ DEFAULT now(),
  updated_at TIMESTAMPTZ DEFAULT now()
);

CREATE TABLE campaigns (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  project_id UUID REFERENCES projects(id) ON DELETE CASCADE,
  name TEXT NOT NULL,
  pocket_config JSONB, -- pocket sélectionnée parmi pockets_detected
  created_at TIMESTAMPTZ DEFAULT now()
);

CREATE TABLE phases (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  campaign_id UUID REFERENCES campaigns(id) ON DELETE CASCADE,
  type TEXT NOT NULL CHECK (type IN ('hit_discovery', 'hit_to_lead', 'lead_optimization')),
  status TEXT DEFAULT 'active' CHECK (status IN ('active', 'frozen', 'completed')),
  frozen_at TIMESTAMPTZ,
  created_at TIMESTAMPTZ DEFAULT now()
);

CREATE TABLE runs (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  phase_id UUID REFERENCES phases(id) ON DELETE CASCADE,
  type TEXT NOT NULL CHECK (type IN ('import', 'calculation', 'generation')),
  -- Pour run calcul : liste des calculs cochés
  calculation_types TEXT[], -- ex: {'docking','admet','scoring','confidence'}
  status TEXT DEFAULT 'created' CHECK (status IN ('created', 'queued', 'running', 'completed', 'failed')),
  config JSONB NOT NULL, -- paramètres par type de calcul
  input_molecule_ids UUID[],
  input_source TEXT, -- 'file', 'internal', 'connected_library'
  input_file_path TEXT,
  progress INTEGER DEFAULT 0,
  current_step TEXT,
  estimated_duration_seconds INTEGER,
  started_at TIMESTAMPTZ,
  completed_at TIMESTAMPTZ,
  error_message TEXT,
  archived BOOLEAN DEFAULT false,
  created_at TIMESTAMPTZ DEFAULT now()
);

CREATE TABLE molecules (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  phase_id UUID REFERENCES phases(id) ON DELETE CASCADE,
  smiles TEXT NOT NULL,
  canonical_smiles TEXT NOT NULL,
  name TEXT,
  source_run_id UUID REFERENCES runs(id),
  bookmarked BOOLEAN DEFAULT false,
  generation_level INTEGER DEFAULT 0,
  parent_molecule_id UUID REFERENCES molecules(id),
  ai_generated BOOLEAN DEFAULT false,
  created_at TIMESTAMPTZ DEFAULT now(),
  UNIQUE(phase_id, canonical_smiles)
);

CREATE TABLE molecule_properties (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  molecule_id UUID REFERENCES molecules(id) ON DELETE CASCADE,
  run_id UUID REFERENCES runs(id),
  property_name TEXT NOT NULL,
  property_value JSONB NOT NULL,
  created_at TIMESTAMPTZ DEFAULT now(),
  UNIQUE(molecule_id, property_name, run_id)
);

CREATE TABLE calculation_cache (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  cache_key TEXT UNIQUE NOT NULL, -- hash(smiles + calc_type + params)
  result JSONB NOT NULL,
  created_at TIMESTAMPTZ DEFAULT now()
);

CREATE TABLE artifacts (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  run_id UUID REFERENCES runs(id) ON DELETE CASCADE,
  type TEXT NOT NULL, -- 'pose_pdb', 'report_pdf', 'synth_tree_json', etc.
  storage_path TEXT NOT NULL,
  created_at TIMESTAMPTZ DEFAULT now()
);

CREATE TABLE audit_log (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  user_id UUID REFERENCES auth.users(id),
  action TEXT NOT NULL,
  entity_type TEXT,
  entity_id UUID,
  details JSONB,
  created_at TIMESTAMPTZ DEFAULT now()
);

-- Index
CREATE INDEX idx_molecules_phase ON molecules(phase_id);
CREATE INDEX idx_molecules_canonical ON molecules(phase_id, canonical_smiles);
CREATE INDEX idx_molecules_parent ON molecules(parent_molecule_id);
CREATE INDEX idx_mol_props_molecule ON molecule_properties(molecule_id);
CREATE INDEX idx_mol_props_name ON molecule_properties(property_name);
CREATE INDEX idx_runs_phase ON runs(phase_id);
CREATE INDEX idx_cache_key ON calculation_cache(cache_key);
CREATE INDEX idx_audit_created ON audit_log(created_at);
```

Performance EAV : index + vue matérialisée pivot + objectif 100K mols/phase. Tout doit être extrêmement scalable.

### 10.2 Auth — Supabase Auth
Email/password, magic link, RLS toutes tables. SSO V10.

### 10.3 Stockage — Supabase Storage
Buckets : structures (PDB files), artifacts (poses, rapports), exports.

### 10.4 Backend — Celery + Redis
Conservé. Progress tracking via task state. Chaque type de calcul dans un run calcul multi-coché est une sous-tâche Celery chaînée.

### 10.5 Moteurs de docking
Choix par run calcul (quand docking coché) :
- **Vina** : CPU, classique, rapide
- **GNINA CPU** : CNN-scored, plus précis
- **GNINA GPU** : RunPod serverless, ~10x plus rapide
- **DiffDock** : deep learning, meilleur pour poses flexibles

Si RUNPOD_API_KEY absente et GPU sélectionné → erreur explicative claire (pas de fallback silencieux).

---

## 11. Legacy V8

Projets V8 en lecture seule, vue legacy séparée.
Route `/legacy/projects/:id`. Pas de migration auto. V9 = clean slate.

---

## 12. Mapping features existantes → V9

| Feature V8 | Composant actuel | V9 |
|------------|-----------------|-----|
| Job pipeline monolithique | tasks.py | Runs par type (import/calcul/génération) |
| InputForm | InputForm.jsx | RunCreator modal (3 étapes) |
| Target preview | InputForm.jsx (TargetPreviewCard) | Rubrique Target du projet (TargetPreviewCard) |
| Structure auto-detection | structure.py | Conservé dans Target (PDB/AlphaFold/ESMFold) |
| Sequence input | models.py (JobCreate) | Input type dans Target (UniProt/Séquence/PDB) |
| ResultsDashboard | ResultsDashboard.jsx | Dashboard Phase cumulatif |
| Viewer3D | Viewer3D.jsx | Drawer 40% + plein écran, toutes options conservées |
| MoleculeCard | MoleculeCard.jsx | Drawer 3D (MoleculeDrawer) |
| HitSelector | HitSelector.jsx | Bookmark + Freeze |
| ParetoFront | ParetoFront.jsx | Conservé, scope phase |
| OptimizationView | OptimizationView.jsx | Conservé, scope phase |
| ProLIF interactions | InteractionAnalysis.jsx | Run calcul (enrichment coché) |
| InteractionDiagram | InteractionDiagram.jsx | Conservé dans Viewer3D et popups |
| Scaffolds | ScaffoldAnalyzer.jsx | Run calcul (clustering coché) |
| ADMET | ADMETRadar.jsx + admet.py | Run calcul (ADMET coché), radar conservé dans popup |
| Safety report | SafetyReport.jsx | Popup depuis dashboard (clic safety_color_code) |
| Confidence scoring | ConfidenceBreakdown.jsx + confidence.py | Run calcul (confidence coché), popup depuis dashboard |
| Retrosynthesis | SynthesisTree.jsx + retrosynthesis.py | Run calcul (retrosynthesis coché), popup depuis dashboard |
| Off-target selectivity | off_target.py | Run calcul (off-target coché) |
| Scoring weights | ScoringWeightsEditor.jsx + scoring.py | RunCreator quand scoring coché |
| Filter pharma | filter_pharma.py | Pré-filtre disponible dans run import |
| Filter shape | filter_shape.py | Pré-filtre disponible dans run import |
| DiffDock | docking_diffdock.py | Option moteur de docking dans run calcul |
| GNINA GPU | docking_gpu.py | Option moteur de docking dans run calcul |
| Massive screening | screening_massive.py | Run calcul sur grande librairie (pas de mode spécial, juste plus long) |
| Génération de novo | generation.py | Run génération (mode batch + mode molécule) |
| Lead optimization | lead_optimization.py | Run génération itérative |
| Target assessment | TargetAssessment.jsx + target_assessment.py | Niveau campagne (TargetInfoCard dans sidebar) |
| Agent IA | AgentChat.jsx + CampaignAgentPanel.jsx | Panel chat campagne, enrichi |
| AnalysisToggles | AnalysisToggles.jsx | Checkboxes dans RunCreator (run calcul) |
| PedagogicalTip | PedagogicalTip.jsx | InfoTip sur chaque concept |
| Methodology | MethodologyPage.jsx | Page conservée (références scientifiques) |
| Estimation temps | TimeEstimate.jsx + time_estimation.py | Étape 3 du RunCreator |
| Notifications email | notifications.py | Settings projet, runs longs |
| Audit log | audit_log.py | Conservé, table audit_log |
| Report PDF | ReportPreview.jsx + report.py | Export PDF dans ExportModal |
| Export | ExportModal.jsx | Conservé, scope = sélection cochée |
| Activity timeline | ActivityTimeline.jsx | RunHistory (historique des runs) |
| Pipeline summary | PipelineSummary.jsx | Remplacé par RunHistory + progression |
| Download buttons | DownloadButtons.jsx | Intégré dans ExportModal |
| ColumnSelector | ColumnSelector.jsx | Conservé dans dashboard (masquer/afficher colonnes) |
| FilterBar | FilterBar.jsx | Conservé dans dashboard |
| Badge | Badge.jsx | Conservé (composant UI réutilisable) |
| Toast | Toast.jsx | Conservé (notifications UI) |
| ProgressBar | ProgressBar.jsx | Conservé (progression runs) |

---

## 13. Composants frontend

### Conservés (réutilisés tels quels ou avec adaptations mineures)
Viewer3D, MoleculeCard, ParetoFront, OptimizationView, InteractionDiagram, InteractionView, ScaffoldAnalyzer, ADMETRadar, SafetyReport, ConfidenceBreakdown, SynthesisTree, ScoringWeightsEditor, InfoTip, PedagogicalTip, ColumnSelector, FilterBar, Badge, Toast, ProgressBar, MethodologyPage, ReferencesPage

### Modifiés
- ResultsDashboard → **PhaseDashboard** (tableau cumulatif par phase)
- InputForm → **RunCreator** (modal 3 étapes)
- Sidebar → **SidebarLayout** (arbre Project > Target > Campaign > Phase + Experimental Results)
- AgentChat/CampaignAgentPanel → adapté au scope campagne multi-phases
- ExportModal → enrichi (PDF rapport + raccourcis)
- MoleculeCard → **MoleculeDrawer** (drawer latéral 40%)

### Nouveaux
- **TargetPreviewCard** : résumé target dans le projet (structure, poches, ChEMBL)
- **PhaseCreator** : création de phase (type + config)
- **SelectionToolbar** : barre d'actions sur sélection (New Run, Export, Bookmark All)
- **RunProgress** : bandeau progression avec étape courante
- **RunHistory** : liste des runs passés + statut + durée
- **FreezeDialog** : confirmation freeze avec warnings
- **QueueIndicator** : message file d'attente si run déjà en cours
- **SafetyPopup** : popup détaillée profil sécurité
- **SynthesisPopup** : popup arbre rétrosynthétique
- **ConfidencePopup** : popup breakdown confiance
- **GenerationInterface** : interface dédiée génération (mode batch + mode molécule)

---

## 14. API Endpoints

### V9
```
-- Projects
POST   /api/v9/projects
GET    /api/v9/projects
GET    /api/v9/projects/:id
PUT    /api/v9/projects/:id
DELETE /api/v9/projects/:id
POST   /api/v9/projects/:id/preview-target  (lance le Target Preview)

-- Campaigns
POST   /api/v9/projects/:id/campaigns
GET    /api/v9/projects/:id/campaigns
GET    /api/v9/campaigns/:id
PUT    /api/v9/campaigns/:id

-- Phases
POST   /api/v9/campaigns/:id/phases
GET    /api/v9/campaigns/:id/phases
GET    /api/v9/phases/:id
PUT    /api/v9/phases/:id
POST   /api/v9/phases/:id/freeze
POST   /api/v9/phases/:id/unfreeze

-- Runs
POST   /api/v9/phases/:id/runs
GET    /api/v9/phases/:id/runs
GET    /api/v9/runs/:id
GET    /api/v9/runs/:id/progress

-- Molecules
GET    /api/v9/phases/:id/molecules
PUT    /api/v9/molecules/:id/bookmark
GET    /api/v9/molecules/:id/safety          (profil sécurité détaillé)
GET    /api/v9/molecules/:id/synthesis        (arbre rétrosynthétique)
GET    /api/v9/molecules/:id/confidence       (breakdown confiance)

-- Export
POST   /api/v9/phases/:id/export             (CSV/SDF/PDF)

-- Agent
POST   /api/v9/campaigns/:id/agent/chat
GET    /api/v9/campaigns/:id/agent/recommendations

-- Pareto
POST   /api/v9/phases/:id/pareto

-- Audit
GET    /api/v9/audit-log
```

### Legacy
```
GET /api/legacy/projects
GET /api/legacy/projects/:id
GET /api/legacy/projects/:id/jobs
GET /api/legacy/jobs/:id/results
```

---

## 15. Backlog (estimations en heures Claude Code)

**Epic 1 — Fondations Supabase** : setup, schéma, migration depuis SQLAlchemy, auth Supabase, RLS

**Epic 2 — Modèle Project/Campaign/Phase** : CRUD, Target Preview (structure.py + pockets.py + ChEMBL), freeze/unfreeze

**Epic 3 — Runs par type** : import (fichier + interne + connecté + pré-filtres), calcul multi-coché (docking/ADMET/scoring/enrichment/clustering/off-target/confidence/retrosynthesis/safety), génération (batch + molécule), cache, queue, progress tracking

**Epic 4 — Dashboard Phase** : table cumulative EAV, sélection, bookmark, tri, filtres, column selector, viewer 3D drawer, popups (safety/synthesis/confidence)

**Epic 5 — Navigation UX** : sidebar arbre, RunCreator modal 3 étapes, RunProgress, RunHistory, QueueIndicator, InfoTips partout, notifications email

**Epic 6 — Polish & intégration** : legacy view, agent IA campagne, Pareto, export (CSV/SDF/PDF), audit log, tests

---

## 16. Exemple — EGFR Kinase

**Projet** : EGFR Inhibitors  
**Target** : UniProt P00533 → auto-détecté PDB 1M17 (expérimental, 2.6Å, co-cristallisé erlotinib), 5 poches détectées, 3400 actifs ChEMBL  

**Campagne** : ATP Pocket (pocket #1, druggability 0.92)

**Phase A** : import 20000 mols SDF → run calcul [docking GNINA GPU + ADMET + confidence] ~30min → run calcul [scoring] → dashboard : trier par composite_score, filtrer confidence > 0.7, bookmarker top 50 → freeze

**Phase B** : run génération batch sur 50 hits (3 itérations × 50 mols) → 150 nouvelles mols taguées "IA generated" avec docking+ADMET auto → Pareto (affinity vs ADMET) → 15 leads → freeze

**Phase C** : run calcul [off-target + retrosynthesis + safety] sur 15 leads → popups safety : vérifier hERG/hepatotox → export PDF rapport complet

---

## 17. Risques & mitigations

| Risque | Mitigation |
|--------|------------|
| EAV lent > 5k mols | Vue matérialisée refresh async + pagination côté serveur + objectif scalabilité maximale |
| Migration Supabase | Schéma prêt, bonne doc, pas de migration auto (V9 = clean slate) |
| Navigation 4 niveaux | Sidebar compacte + header persistant avec breadcrumb |
| Freeze edge cases | Warnings explicites + confirmation dialog |
| RunPod indisponible / cold start | Erreur explicative claire avec message d'aide, PAS de fallback silencieux |
| Trop de colonnes | ColumnSelector + colonnes masquables |
| Run calcul multi-coché complexe | Sous-tâches Celery chaînées, progress par sous-tâche |
| DiffDock installation lourde | Container Docker dédié, fallback = option non disponible dans UI |

---

## 18. Hors scope V9 (V10+)

- Model Factory ML (prédiction activité sans docking)
- Multi-cibles par projet (polypharmacologie)
- Phases multiples du même type
- SSO / SAML
- Concurrence runs illimitée
- Edge Functions pour micro-calculs
- Knowledge Layer (MegaPharmaDB)
- Intégration experimental results dans les phases
- PharmacoDB (extrait en projet standalone, fonctionnalité future via Knowledge Layer)

---

## 19. Décisions consolidées

| # | Décision |
|---|----------|
| D1 | Cible = niveau projet (rubrique Target avec auto-detection) |
| D2 | 3 inputs target : UniProt ID, Séquence FASTA, PDB ID |
| D3 | Target Preview auto à la création (structure, poches, ChEMBL) |
| D4 | Campagne = pocket sélectionnée parmi celles détectées |
| D5 | 1 campagne par défaut, ajout possible |
| D6 | Phases création manuelle, 1 par type en V9 |
| D7 | Freeze réversible avec warning |
| D8 | 3 types de run : import, calcul (multi-coché), génération |
| D9 | Run calcul : checkboxes pour sélectionner les calculs |
| D10 | Docking : 4 moteurs (Vina, GNINA CPU, GNINA GPU, DiffDock) |
| D11 | Pas de multi-source par run |
| D12 | 1 molécule = 1 ligne dédupliquée par SMILES canonique |
| D13 | Colonnes affichées par défaut, masquables par l'utilisateur |
| D14 | Viewer 3D : drawer 40% + plein écran |
| D15 | Popups détaillées : safety, retrosynthesis, confidence |
| D16 | InfoTip pédagogique sur chaque concept |
| D17 | Scoring : poids configurés par run (pas de défauts campagne) |
| D18 | Agent IA : niveau campagne, cross-phases |
| D19 | Pareto : niveau phase |
| D20 | Export = sélection cochée (CSV/SDF/PDF) |
| D21 | Supabase Cloud (PostgreSQL + Auth + Storage) |
| D22 | Celery + Redis conservés |
| D23 | Legacy V8 lecture seule |
| D24 | Archive seulement (pas de suppression de runs) |
| D25 | RunPod erreur explicite si indisponible (pas de fallback) |
| D26 | Notifications email pour runs longs |
| D27 | Audit log conservé |
| D28 | Cache par hash (smiles + type + params) |
| D29 | Experimental Results = placeholder V9 |
| D30 | Génération : mode batch + mode molécule individuelle |

---

*Fin du CDC BindX V9 — Prêt pour implémentation.*
