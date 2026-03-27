3# AFVS — Infrastructure AWS : Setup Complet

> Document de handoff pour la suite du développement (BLOC 4-7 de la spec AFVS)
> Date : 2026-03-26

## Contexte

Le BLOC 3 (Prérequis infra) de `docs/AFVS_SPEC_BINDX_V9.md` a été réalisé via la console AWS.
Ce document détaille toutes les ressources créées et leur état actuel.

---

## 1. Compte AWS

- **Account ID** : `890004513839`
- **Région** : `us-east-1` (Virginie du Nord)
- **Plan** : Upgradé du Free Tier vers plan payant (nécessaire pour t3.medium)
- **User** : anthony.boisbouvier@hotmail.fr

---

## 2. Ressources AWS créées

### 2.1 EC2 — AFVS-Controller

| Propriété | Valeur |
|---|---|
| Instance ID | `i-073123fd8b2d63010` |
| Nom | `AFVS-Controller` |
| Type | `t3.medium` (2 vCPU, 4 GiB RAM) |
| AMI | Ubuntu 22.04 LTS (Jammy) |
| IP publique | `54.174.165.245` |
| DNS public | `ec2-54-174-165-245.compute-1.amazonaws.com` |
| IP privée | `172.31.92.22` |
| Key pair | `afvs-ec2-key` (RSA, fichier .pem dans `C:\Users\antho\Downloads\afvs-ec2-key.pem`) |
| Security Group | `launch-wizard-1` (SSH port 22 ouvert sur 0.0.0.0/0) |
| IAM Role | `AFVS-Controller-Role` (attaché) |
| Stockage | 8 GiB gp3 (défaut) |

**Connexion SSH** :
```bash
ssh -i afvs-ec2-key.pem ubuntu@54.174.165.245
```
Ou via EC2 Instance Connect dans la console AWS (testé et fonctionnel).

### 2.2 S3 — Bucket de sorties AFVS

| Propriété | Valeur |
|---|---|
| Nom | `bindx-afvs-outputs` |
| Région | `us-east-1` |
| État | Vide, accessible depuis l'EC2 via le rôle IAM |

### 2.3 IAM Role — AFVS-Controller-Role

| Propriété | Valeur |
|---|---|
| Nom | `AFVS-Controller-Role` |
| ARN | `arn:aws:iam::890004513839:role/AFVS-Controller-Role` |
| Type | EC2 (trusted entity: ec2.amazonaws.com) |
| Policies attachées | `AWSBatchFullAccess`, `AmazonS3FullAccess` |
| Attaché à | EC2 `i-073123fd8b2d63010` |

L'EC2 utilise ce rôle automatiquement (pas besoin d'access keys). Vérifié avec `aws sts get-caller-identity`.

### 2.4 Key Pair

| Propriété | Valeur |
|---|---|
| Nom | `afvs-ec2-key` |
| Type | RSA |
| Fichier local | `C:\Users\antho\Downloads\afvs-ec2-key.pem` |

**IMPORTANT** : Ce fichier .pem est nécessaire pour la connexion SSH depuis le backend BindX (via paramiko). Il doit être copié sur le serveur backend et référencé dans la variable d'env `AFVS_EC2_SSH_KEY_PATH`.

---

## 3. Logiciels installés sur l'EC2

### 3.1 Paquets système
```
python3-pip, python3-virtualenv, git, docker.io, awscli
```

### 3.2 AdaptiveFlow (AFVS)
- **Repo** : `~/AFVS/` (cloné depuis `https://github.com/LigandUniverse/AFVS/`)
- **Branche** : `develop` (branche par défaut)
  - Note : la branche `afvs-2` mentionnée dans la spec n'existe PAS sur le repo. Branches disponibles : `develop`, `feature/slurm-aws`, `nidhin_afvs`
- **Scripts disponibles** dans `~/AFVS/tools/` :
  - `afvs_submit_jobs.py` — soumission de jobs AWS Batch
  - `afvs_build_docker.sh` — build + push image Docker vers ECR
  - `afvs_get_status.py` — statut des jobs
  - `afvs_get_top_results.py` — récupération des meilleurs résultats
  - `afvs_prepare_folders.py`, `afvs_prepare_workunits.py`
  - `afvs_killall_awsbatch.py`
  - et d'autres (postprocess, gen_sparse, verify_collections...)

### 3.3 Virtualenv Python
- **Chemin** : `~/afvs_env/`
- **Activation** : `source ~/afvs_env/bin/activate`
- **Packages installés** :
  - boto3 1.42.77
  - pandas 2.3.3
  - pyarrow 23.0.1
  - jinja2 3.1.6
  - numpy 2.2.6
  - botocore, s3transfer, etc.

### 3.4 Docker
- **Version** : Docker 28.2.2
- **User ubuntu** : ajouté au groupe `docker`

### 3.5 AWS CLI
- **Région par défaut** : `us-east-1` (configurée via `aws configure set default.region us-east-1`)
- **Credentials** : via le rôle IAM (pas d'access keys, tout passe par le metadata service EC2)

---

## 4. Variables d'environnement pour le backend BindX

Comme décrit dans BLOC 7 de la spec, voici les valeurs à configurer dans le `.env` du backend :

```env
# AFVS EC2 Controller
AFVS_EC2_HOST=54.174.165.245
AFVS_EC2_USER=ubuntu
AFVS_EC2_SSH_KEY_PATH=/path/to/afvs-ec2-key.pem  # À copier sur le serveur backend

# AFVS S3
AFVS_S3_BUCKET=bindx-afvs-outputs
AFVS_S3_REGION=us-east-1

# AWS (si le backend a besoin d'accéder directement à S3)
AWS_ACCESS_KEY_ID=<à créer si nécessaire>
AWS_SECRET_ACCESS_KEY=<à créer si nécessaire>
AWS_DEFAULT_REGION=us-east-1

# AFVS paths sur l'EC2
AFVS_INSTALL_DIR=/home/ubuntu/AFVS
AFVS_VENV_PATH=/home/ubuntu/afvs_env
```

---

## 5. Infrastructure déployée (2026-03-26)

### CloudFormation Stacks

| Stack | Statut | Ressources |
|-------|--------|------------|
| `af-vpc` | CREATE_COMPLETE | VPC, subnets, NAT gateway, Internet gateway |
| `af` | CREATE_COMPLETE | 3 CE (Spot), Job Queue, 2 Job Definitions, 2 ECR repos |

### AWS Batch

| Ressource | Nom | Détails |
|-----------|-----|---------|
| Compute Env 1 | `af-CE1` | Spot, c5 instances, max 16 vCPUs |
| Compute Env 2 | `af-CE2` | Spot, c5 instances, max 16 vCPUs |
| Compute Env 3 | `af-CE3` | Spot, c5 instances, max 16 vCPUs |
| Job Queue | `af-queue1` | Priority 1, linked to CE1/CE2/CE3 |
| Job Definition | `af-jobdef-afvs` | Image ECR af-afvs-ecr:latest, timeout 7200s |
| Job Definition | `af-jobdef-aflp` | Image ECR af-aflp-ecr:latest (unused for now) |

### ECR

| Repo | URI |
|------|-----|
| `af-afvs-ecr` | `890004513839.dkr.ecr.us-east-1.amazonaws.com/af-afvs-ecr:latest` |
| `af-aflp-ecr` | `890004513839.dkr.ecr.us-east-1.amazonaws.com/af-aflp-ecr` (empty) |

Image Docker buildée et pushée (sha256:cf60d0594fd0...). Contient : smina, qvina02, gwovina, vina, idock, etc.

### Développement BindX (DONE)

| Bloc | Statut | Détails |
|------|--------|---------|
| BLOC 4 — DB | DONE | Table `afvs_jobs`, ORM `AFVSJobORM_V9`, migration appliquée |
| BLOC 5 — Backend | DONE | `pipeline/afvs.py` (AFVSController SSH), `routers/v9/afvs.py` (3 endpoints), `tasks_v9.py` (_run_afvs) |
| BLOC 6 — Frontend | DONE | `AFVSConfigForm.jsx`, RunCreator integration, PhaseDashboard wiring |
| BLOC 7 — Config | DONE | Docker env vars, SSH key mount, `.env` variables |

## 5b. Ce qui reste

1. **Accès librairie REAL Space** — Demande à faire sur virtual-flow.org avec Account ID `890004513839`
   - Sans cet accès, le bucket S3 data (`AFVS_S3_DATA_BUCKET`) est vide → les jobs n'ont pas de ligands à docker
2. **Test E2E** — Lancer 1 workunit de test après obtention REAL Space

---

## 6. Architecture de connexion

```
BindX Backend (Celery worker)
    │
    │ SSH (paramiko, port 22, clé afvs-ec2-key.pem)
    ▼
EC2 AFVS-Controller (54.174.165.245)
    │
    │ AWS SDK (boto3, via IAM Role)
    ├──► AWS Batch (submit/monitor jobs)
    ├──► S3 bindx-afvs-outputs (read results)
    └──► ECR (pull Docker image)
         │
         ▼
    AWS Batch Spot Instances
    (exécutent le docking AFVS dans des containers Docker)
         │
         │ Résultats écrits dans S3
         ▼
    S3 bindx-afvs-outputs/
         │
         │ Récupérés par le backend via boto3 ou SSH+script
         ▼
    Import dans la table molecules (dashboard phase)
```

---

## 7. Sécurité — Points d'attention

- Le Security Group `launch-wizard-1` ouvre SSH sur `0.0.0.0/0` — à restreindre en production à l'IP du backend BindX uniquement
- La clé SSH `.pem` ne doit jamais être commitée dans le repo
- Le rôle IAM a des permissions larges (FullAccess) — à affiner en production avec des policies plus restrictives
- L'IP publique de l'EC2 peut changer au reboot — envisager une Elastic IP pour la production
