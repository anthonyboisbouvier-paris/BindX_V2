# BindX V9 — AFVS Module Specification
## Adaptive Focused Virtual Screening via AdaptiveFlow + AWS Batch

**Version** : 2.2 — COMPLÈTE (presets stratégie, cap top_n, Celery long-running)
**Date** : 25 mars 2026
**Statut** : ✅ Tous les blocs validés
**Linear** : BDX-41

---

## BLOC 1 — Vision & Périmètre ✅

L'AFVS dans BindX est une **intégration légère** : BindX pilote AdaptiveFlow (outil open-source, VirtualFlow 2.0) qui tourne sur AWS Batch. BindX ne réimplémente aucun algorithme de docking — il orchestre l'outil existant, attend les résultats, et les importe dans le dashboard Phase A.

**Du point de vue utilisateur :**
1. Il configure son run AFVS dans BindX (cible, filtres, stratégie)
2. Il clique "Lancer"
3. Il attend quelques heures (durée à calibrer au premier run de test)
4. Les top hits apparaissent dans son dashboard Phase A avec leurs scores

**Ce que BindX code :**
- Un formulaire de config AFVS (RunCreator)
- Une tâche Celery qui pilote l'EC2 AFVS via SSH
- Un polling qui surveille la fin du job sur S3
- Un import des résultats dans le dashboard

**Modèle de coût :**
- BindX utilise un **compte AWS centralisé** (pas le compte de l'utilisateur)
- Les jobs tournent sur des **instances Spot AWS Batch** (~70–90% moins cher que On-Demand), possible car AdaptiveFlow est conçu pour les interruptions : chaque work unit est indépendant et redémarrable automatiquement
- BindX facture le coût AWS réel **× 2** à l'utilisateur
- Coût fixe mensuel de l'infra : ~$30/mois (EC2 controller) + ~$2/mois (S3)
- Coût variable par run : à calibrer lors du premier run de test

---

## BLOC 2 — Architecture ✅

```
Utilisateur BindX
      │ Lance un Run AFVS
      ▼
BindX Backend (FastAPI + Celery)
      │ SSH (paramiko)
      ▼
EC2 "AFVS Controller" (t3.medium, toujours allumée ~$30/mois)
  └── AdaptiveFlow installé + AWS CLI configuré
      │ Soumet des jobs AWS Batch
      ▼
AWS Batch (instances Spot CPU, scale automatique 0 → N)
  ├── Lit la librairie depuis S3 VirtualFlow (Enamine REAL Space, 69B mol)
  └── Écrit les résultats dans S3 BindX
      │
      ▼
S3 Bucket BindX — résultats (top_results.csv)
      │ Polling Celery toutes les 60 secondes
      ▼
BindX importe les top hits → Dashboard Phase A
```

**Détection de fin de job :** polling S3 toutes les 60 secondes. Sur un run de plusieurs heures, 1 minute de latence est négligeable.

---

## BLOC 3 — Prérequis infra ✅

Setup à faire **une seule fois** par Anthony.

**Étape 1 — Compte AWS**
Créer un compte AWS sur aws.amazon.com. C'est le compte centralisé BindX.

**Étape 2 — Accès à la librairie REAL Space**
Demander l'accès sur virtual-flow.org/real-space-2022q12. Gratuit, délai ~1 jour ouvré. ⚠️ Nécessite un AWS Account ID → faire étape 1 d'abord.

**Étape 3 — EC2 AFVS Controller**
Lancer une instance EC2 t3.medium (Ubuntu 22.04, région us-east-1), allumée en permanence (~$30/mois).

**Étape 4 — Installer AdaptiveFlow sur l'EC2**
```bash
git clone -b afvs-2 https://github.com/LigandUniverse/AFVS/
cd AFVS
python3 -m virtualenv ~/afvs_env && source ~/afvs_env/bin/activate
pip install boto3 pandas pyarrow jinja2
aws configure
./tools/afvs_build_docker.sh
```

**Étape 5 — Bucket S3 outputs**
Créer un bucket S3 `bindx-afvs-outputs` (région us-east-1).

**Étape 6 — Rôles IAM**
- EC2 Controller : `AWSBatchFullAccess` + `AmazonS3FullAccess` (outputs) + `AmazonS3ReadOnlyAccess` (librairie VirtualFlow)
- Jobs AWS Batch : S3 read (librairie) + S3 write (outputs)

**Étape 7 — Run de test manuel**
Lancer 1 work unit depuis l'EC2 (Tutorial 1 AdaptiveFlow) pour valider la chaîne et calibrer les coûts.

---

## BLOC 4 — Base de données ✅

**1 modification sur `runs` :**
```sql
ALTER TABLE runs DROP CONSTRAINT IF EXISTS runs_type_check;
ALTER TABLE runs ADD CONSTRAINT runs_type_check
  CHECK (type IN ('import', 'calculation', 'generation', 'afvs'));
```

**1 nouvelle table `afvs_jobs` :**
```sql
CREATE TABLE afvs_jobs (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  run_id UUID UNIQUE REFERENCES runs(id) ON DELETE CASCADE,
  receptor_s3_path TEXT NOT NULL,
  docking_box JSONB NOT NULL,
  docking_scenario TEXT DEFAULT 'smina_vinardo',
  top_n_results INTEGER DEFAULT 10000,
  budget_cap_usd FLOAT,               -- plafond fixé par l'utilisateur
  s3_output_prefix TEXT,
  afvs_work_units INTEGER,
  afvs_phase TEXT DEFAULT 'prescreen'
    CHECK (afvs_phase IN ('prescreen', 'primary_screen', 'completed', 'failed')),
  molecules_screened BIGINT,
  molecules_imported INTEGER,
  actual_cost_usd FLOAT,
  created_at TIMESTAMPTZ DEFAULT now(),
  updated_at TIMESTAMPTZ DEFAULT now()
);
CREATE INDEX idx_afvs_jobs_run_id ON afvs_jobs(run_id);
```

---

## BLOC 5 — Backend ✅

**Fichiers :**
```
backend/pipeline/afvs.py     ← AFVSController (SSH) + utilitaires
backend/routers/v9/afvs.py   ← 3 endpoints
backend/tasks_v9.py          ← afvs_run_task + afvs_cancel_task
```

**AFVSController** : classe paramiko SSH. Méthodes : `run()`, `write_file()`, `build_all_ctrl()`, `upload_receptor_to_s3()`, `import_top_results()`. Protéine préparée + docking box récupérées automatiquement depuis la Campaign.

**afvs_run_task** (Celery) :
```
1. Récupérer protéine préparée + docking box depuis Campaign
2. Upload PDB → S3
3. SSH → all.ctrl + afvs_prepare_folders + afvs_prepare_workunits
4. SSH → afvs_submit_jobs (prescreen)
5. Poll S3/60s → fin prescreen (vérification budget_cap à chaque poll)
6. SSH → postprocess + afvs_submit_primary_screen
7. Poll S3/60s → fin primary screen (vérification budget_cap à chaque poll)
8. SSH → afvs_get_top_results --top N → CSV dans S3
9. Import CSV → molécules dans Phase A
10. Écrire actual_cost_usd
```

**afvs_cancel_task** : SSH → annuler jobs Batch + marquer run failed + écrire coût partiel.

**Endpoints :**
```
POST   /api/v9/phases/{phase_id}/runs/afvs   → lancer
DELETE /api/v9/runs/{run_id}/afvs            → annuler
GET    /api/v9/runs/{run_id}/afvs/status     → état détaillé
```

---

## BLOC 6 — Frontend ✅

### RunCreator — Étape 1
Carte "AFVS — Ultra-Large Screening (69B mol)" visible uniquement en Phase A.

### RunCreator — Étape 2 : Configuration complète

**Infos en lecture seule (depuis Campaign) :**
- Cible + nom de la pocket sélectionnée
- Protéine préparée utilisée

**Paramètres exposés à l'utilisateur :**

| Paramètre | Défaut | Description |
|---|---|---|
| Stratégie | Balanced ⭐ | Lean / Balanced / Aggressive (voir presets ci-dessous) |
| Scénario docking | smina_vinardo | smina_vinardo / qvina2 |
| Exhaustiveness prescreen | 1 | Slider 1–8 |
| Exhaustiveness primary | 8 | Slider 1–16 |
| Nb reps/tranche | 1 | Slider 1–10 |
| % tranches sélectionnées | 0.5% | Slider 0.1–2% |
| Filtre MW max | 500 Da | |
| Filtre logP max | 5 | |
| Filtre HBD max | 5 | |
| Filtre HBA max | 10 | |
| Top molécules à importer | 10 000 | Champ libre, max 100 000 |
| **Plafond de coût (USD)** | **500** | **L'utilisateur fixe son budget max. Le run s'arrête proprement si dépassé.** |

**Presets Stratégie** (sélectionner une stratégie pré-remplit les sliders, modifiables ensuite) :

| Paramètre | Lean | Balanced ⭐ | Aggressive |
|---|---|---|---|
| Exhaustiveness prescreen | 1 | 1 | 2 |
| Exhaustiveness primary | 4 | 8 | 16 |
| Reps/tranche | 1 | 1 | 3 |
| % tranches sélectionnées | 0.1% | 0.5% | 2% |
| Top N import | 5 000 | 10 000 | 50 000 |
| Budget cap (USD) | 100 | 500 | 2 000 |

**Estimation de coût en temps réel** sous les paramètres (indicative, recalculée à chaque changement).
Formule indicative : `coût ≈ nb_work_units × coût_moyen_par_WU`. Le `coût_moyen_par_WU` sera calibré lors du premier run de test (Bloc 3, Étape 7). En attendant, afficher "Estimation indisponible — premier run de calibration requis".

### RunCreator — Étape 3 : Review
Résumé config + estimation coût/durée + avertissement durée + rappel facturation × 2.

### Dashboard Phase A
- Barre de progression avec `runs.current_step`
- Bouton "Annuler" tant que `run.status = 'running'`
- Badge `AFVS` sur les molécules importées
- Coût réel dans RunHistory après complétion

---

## BLOC 7 — Variables d'environnement ✅

```bash
# EC2 AFVS Controller
AFVS_EC2_HOST=ec2-XX-XX-XX-XX.compute-1.amazonaws.com
AFVS_EC2_USER=ubuntu
AFVS_SSH_KEY_PATH=/app/secrets/afvs_ec2_key.pem

# S3
AFVS_S3_JOB_BUCKET=bindx-afvs-outputs
AFVS_S3_DATA_BUCKET=vf-real-space-bucket

# AWS credentials
AWS_ACCESS_KEY_ID=...
AWS_SECRET_ACCESS_KEY=...
AWS_DEFAULT_REGION=us-east-1
```

La clé SSH `afvs_ec2_key.pem` est générée à la création de l'EC2 et stockée hors git.

---

## BLOC 8 — Gestion des erreurs ✅

| Erreur | Comportement |
|---|---|
| EC2 inaccessible (SSH timeout) | `run.status = 'failed'` — "EC2 AFVS Controller inaccessible." |
| AWS Batch quota dépassé | stderr capturé dans `run.error_message` |
| Timeout prescreen (>48h) | `run.status = 'failed'` — "Prescreen timeout. Relancez le run." |
| Timeout primary screen (>72h) | Idem |
| CSV résultats vide | `run.status = 'failed'` — "Aucun résultat retourné." |
| Accès S3 librairie refusé | "Accès REAL Space refusé. Vérifiez vos credentials AWS." |
| Plafond de coût atteint | Run arrêté proprement, molécules déjà dockées importées, message : "Budget cap atteint — résultats partiels importés." |
| Annulation utilisateur | `run.status = 'failed'`, coût partiel écrit dans `actual_cost_usd` |
| Script AdaptiveFlow erreur | stderr → `run.error_message`, visible dans RunHistory |

Dans tous les cas : `actual_cost_usd` est écrit pour facturation même partielle.

---

## BLOC 9 — Checklist d'implémentation ✅

### Infra (one-time, dans l'ordre)
- [ ] Créer compte AWS + noter Account ID (12 chiffres)
- [ ] Demander accès librairie REAL Space (virtual-flow.org/real-space-2022q12)
- [ ] Lancer EC2 t3.medium us-east-1 + Security Group port 22 depuis IP BindX uniquement
- [ ] Installer AdaptiveFlow + builder containers Docker
- [ ] Créer bucket S3 `bindx-afvs-outputs`
- [ ] Configurer rôles IAM
- [ ] Tester run manuel 1 work unit → calibrer coûts réels

### Backend BindX
- [ ] Migration Supabase : contrainte `runs.type` + table `afvs_jobs`
- [ ] `backend/pipeline/afvs.py` — AFVSController + utilitaires
- [ ] `backend/tasks_v9.py` — `afvs_run_task` + `afvs_cancel_task`
- [ ] `backend/routers/v9/afvs.py` — 3 endpoints
- [ ] Variables d'environnement + validation au démarrage
- [ ] Tests : mock SSH (pytest-mock) + mock S3 (moto)

### Frontend BindX
- [ ] Carte "AFVS" dans RunCreator étape 1 (Phase A uniquement)
- [ ] `AFVSConfig.jsx` — formulaire complet + estimation coût temps réel
- [ ] `AFVSSummary.jsx` — review étape 3
- [ ] Badge `AFVS` dans PhaseDashboard
- [ ] Bouton "Annuler" dans RunHistory
- [ ] `actual_cost_usd` dans RunHistory post-complétion
- [ ] InfoTip AFVS sur les molécules

---

## BLOC 10 — Notes importantes ✅

**Protéine + docking box** viennent automatiquement de la Campaign, l'utilisateur ne ressaisit rien.

**Pas de pause** entre prescreen et primary screen : AdaptiveFlow sélectionne automatiquement les tranches selon le `top_tranche_pct` configuré à l'avance.

**Plafond de coût** : paramètre utilisateur (défaut 500 USD). Si atteint pendant le run, le job s'arrête proprement et les molécules déjà dockées sont importées — résultats partiels valables.

**Multi-utilisateurs** : plusieurs runs AFVS peuvent tourner en parallèle sans problème — chaque run a son propre `run_id` et préfixe S3. Pas de limite artificielle. La règle "1 run simultané par user" de V9 ne s'applique PAS aux runs AFVS (run long-running, pas de blocage du workflow).

**Top N cap** : max 100 000 molécules importables. Au-delà, le dashboard Phase A devient difficile à naviguer même avec la virtualisation. L'utilisateur peut relancer un run AFVS avec un top N plus élevé si besoin.

**Celery long-running** : un run AFVS dure 6–72h. Précautions :
- `acks_late=True` sur la task → si le worker crash, la task est redistribuée
- Reconnexion SSH automatique (paramiko keepalive + retry) à chaque commande
- État persisté dans `afvs_jobs` en DB (pas en mémoire Celery) → reprise possible
- `soft_time_limit=259200` (72h) + `time_limit=262800` (73h) comme garde-fou

**Références :**
- AdaptiveFlow doc : https://docs.virtual-flow.org/documentation-af
- GitHub AFVS : https://github.com/LigandUniverse/AFVS
- Tutorial AWS : https://docs.virtual-flow.org/tutorials-af/tutorial-1-afvs-on-aws/introduction
- Gorgulla et al. (2023). VirtualFlow 2.0. bioRxiv 10.1101/2023.04.25.537981
