# BindX V9 — Architecture

> Plateforme de découverte de médicaments in silico.
> Ce document explique comment toutes les pièces du projet s'assemblent.

---

## Vue d'ensemble

```mermaid
graph TB
    subgraph "👤 Utilisateur"
        Browser["🌐 Navigateur web"]
    end

    subgraph "☁️ Internet (hébergé gratuitement)"
        Pages["📄 GitHub Pages<br/><i>Frontend React</i><br/>anthonyboisbouvier-paris.github.io/BindX_V2"]
        Supabase["🗄️ Supabase<br/><i>Base de données + Auth</i><br/>PostgreSQL dans le cloud"]
    end

    subgraph "🔗 Pont Internet ↔ PC"
        Tunnel["🚇 LocalTunnel<br/><i>bindx-api.loca.lt</i><br/>Rend le PC accessible depuis internet"]
    end

    subgraph "💻 Ton PC (WSL + Docker)"
        API["⚡ FastAPI<br/><i>Serveur API</i><br/>Reçoit et traite les requêtes"]
        Celery["🔧 Celery Worker<br/><i>Calculs en arrière-plan</i><br/>Docking, ADMET, Scoring..."]
        Redis["📬 Redis<br/><i>File d'attente</i><br/>Transmet les tâches"]
    end

    subgraph "🌍 APIs externes"
        RCSB["🧬 RCSB PDB<br/><i>Structures 3D de protéines</i>"]
        UniProt["📋 UniProt<br/><i>Infos sur les protéines</i>"]
        ChEMBL["💊 ChEMBL<br/><i>Données de médicaments connus</i>"]
    end

    Browser -->|"1. Affiche le site"| Pages
    Browser -->|"2. Login / Auth"| Supabase
    Pages -->|"3. Appels API"| Tunnel
    Tunnel -->|"4. Transfère au PC"| API
    API -->|"5. Lit/écrit données"| Supabase
    API -->|"6. Lance un calcul"| Redis
    Redis -->|"7. Distribue la tâche"| Celery
    Celery -->|"8. Sauvegarde résultats"| Supabase
    API -->|"9. Cherche structure"| RCSB
    API -->|"10. Infos protéine"| UniProt
    API -->|"11. Médicaments connus"| ChEMBL

    style Pages fill:#1a1a2e,stroke:#00e6a0,color:#fff
    style Supabase fill:#1a1a2e,stroke:#3b82f6,color:#fff
    style Tunnel fill:#1a1a2e,stroke:#f59e0b,color:#fff
    style API fill:#1a1a2e,stroke:#06b6d4,color:#fff
    style Celery fill:#1a1a2e,stroke:#8b5cf6,color:#fff
    style Redis fill:#1a1a2e,stroke:#ef4444,color:#fff
    style RCSB fill:#1a1a2e,stroke:#6b7280,color:#fff
    style UniProt fill:#1a1a2e,stroke:#6b7280,color:#fff
    style ChEMBL fill:#1a1a2e,stroke:#6b7280,color:#fff
    style Browser fill:#0f0f23,stroke:#00e6a0,color:#fff
```

---

## Parcours d'une requête utilisateur

Quand tu cliques "Créer un projet" dans le navigateur, voici ce qui se passe :

```mermaid
sequenceDiagram
    participant U as 👤 Navigateur
    participant F as 📄 GitHub Pages
    participant T as 🚇 LocalTunnel
    participant A as ⚡ FastAPI
    participant S as 🗄️ Supabase

    U->>F: 1. Ouvre le site
    F-->>U: 2. Affiche l'interface React

    U->>S: 3. Login (email + mot de passe)
    S-->>U: 4. Token d'authentification ✅

    U->>T: 5. "Créer un projet" (avec le token)
    T->>A: 6. Transfère la requête au PC
    A->>S: 7. INSERT dans la base de données
    S-->>A: 8. Projet créé ✅
    A-->>T: 9. Réponse JSON
    T-->>U: 10. Affiche le nouveau projet
```

---

## Parcours d'un calcul (ex: ADMET)

Les calculs longs (docking, ADMET, scoring...) sont traités en arrière-plan pour ne pas bloquer l'interface :

```mermaid
sequenceDiagram
    participant U as 👤 Navigateur
    participant A as ⚡ FastAPI
    participant R as 📬 Redis
    participant C as 🔧 Celery
    participant S as 🗄️ Supabase

    U->>A: 1. "Lance un calcul ADMET sur 4 molécules"
    A->>R: 2. Dépose la tâche dans la file
    A-->>U: 3. "OK, calcul lancé" (réponse immédiate)

    R->>C: 4. Celery récupère la tâche
    C->>C: 5. Calcule ADMET pour chaque molécule
    Note over C: LogP, solubilité, toxicité,<br/>biodisponibilité, etc.
    C->>S: 6. Sauvegarde les résultats
    C->>S: 7. Marque le run "completed"

    U->>A: 8. Vérifie le statut (polling)
    A->>S: 9. Lit le statut
    A-->>U: 10. "Terminé ✅ — voici les résultats"
```

---

## Résolution d'une cible protéique

Quand tu configures la cible d'un projet (ex: EGFR), le système cherche automatiquement la meilleure structure 3D :

```mermaid
sequenceDiagram
    participant U as 👤 Navigateur
    participant A as ⚡ FastAPI
    participant UP as 📋 UniProt
    participant PDB as 🧬 RCSB PDB
    participant P2R as 🔍 P2Rank

    U->>A: 1. "Ma cible est P00533 (EGFR)"

    A->>UP: 2. Cherche infos protéine
    UP-->>A: 3. Nom, séquence, références

    A->>PDB: 4. Cherche structures 3D disponibles
    PDB-->>A: 5. Liste de structures (ex: 8A27, 1.07Å)

    A->>PDB: 6. Télécharge le fichier PDB
    PDB-->>A: 7. Fichier 3D de la protéine

    A->>P2R: 8. Détecte les poches (sites actifs)
    P2R-->>A: 9. 2 poches trouvées avec coordonnées [x,y,z]

    A-->>U: 10. Affiche structure + poches
    Note over U: L'utilisateur choisit la poche<br/>où docker les molécules
```

---

## Organisation des données

Le modèle de données suit la logique d'un projet de drug discovery :

```mermaid
graph TD
    Project["🎯 Projet<br/><i>1 cible protéique</i><br/>Ex: EGFR"]
    Campaign["📋 Campagne<br/><i>1 stratégie = cible + poche + règles</i><br/>Ex: Poche 1 du EGFR"]
    PhaseA["🔬 Phase A — Hit Discovery<br/><i>Trouver des molécules actives</i><br/>dans une librairie"]
    PhaseB["⚗️ Phase B — Hit-to-Lead<br/><i>Optimiser les meilleurs hits</i><br/>générer des analogues"]
    PhaseC["💎 Phase C — Lead Optimization<br/><i>Affiner avec ADMET complet</i><br/>préparer pour le labo"]
    Run["▶️ Run<br/><i>1 unité de calcul</i><br/>Import, Docking, ADMET..."]
    Mol["💊 Molécule<br/><i>1 ligne dans le tableau</i><br/>SMILES + propriétés calculées"]

    Project -->|"contient"| Campaign
    Campaign -->|"contient 3 phases"| PhaseA
    Campaign -->|""| PhaseB
    Campaign -->|""| PhaseC
    PhaseA -->|"contient des"| Run
    PhaseA -->|"contient des"| Mol
    PhaseB -->|"contient des"| Run
    PhaseB -->|"contient des"| Mol
    PhaseC -->|"contient des"| Run
    PhaseC -->|"contient des"| Mol

    style Project fill:#1e3a5f,stroke:#00e6a0,color:#fff
    style Campaign fill:#1e3a5f,stroke:#06b6d4,color:#fff
    style PhaseA fill:#1e3a5f,stroke:#3b82f6,color:#fff
    style PhaseB fill:#1e3a5f,stroke:#8b5cf6,color:#fff
    style PhaseC fill:#1e3a5f,stroke:#22c55e,color:#fff
    style Run fill:#2d2d44,stroke:#f59e0b,color:#fff
    style Mol fill:#2d2d44,stroke:#ec4899,color:#fff
```

---

## Les 9 types de calcul

Chaque **Run** effectue un type de calcul. Les résultats apparaissent comme colonnes dans le tableau :

| Type | Ce qu'il fait | Résultats |
|------|--------------|-----------|
| **Import** | Charge des molécules (SMILES, fichier SDF) | Nom, SMILES |
| **Docking** | Simule la fixation molécule ↔ protéine | Score de docking, score CNN |
| **ADMET** | Prédit absorption, distribution, métabolisme, excrétion, toxicité | LogP, solubilité, BBB, biodisponibilité |
| **Scoring** | Calcule un score composite pondéré | Score composite (0-1) |
| **Safety** | Évalue la toxicité (hERG, AMES, hépatotoxicité) | Code couleur sécurité (vert/jaune/rouge) |
| **Confidence** | Vérifie la fiabilité des prédictions (PAINS, domaine d'applicabilité) | Score de confiance (0-1) |
| **Clustering** | Regroupe les molécules similaires par scaffold | Cluster ID, similarité Tanimoto |
| **Enrichment** | Analyse les interactions protéine-ligand (ProLIF) | Nombre de contacts, scaffold |
| **Retrosynthesis** | Évalue la faisabilité de synthèse | Étapes, coût, disponibilité réactifs |
| **Off-target** | Teste la sélectivité (autres protéines) | Score sélectivité, hits off-target |

---

## Déploiement actuel (staging)

```mermaid
graph LR
    subgraph "☁️ GitHub (gratuit)"
        GH["📄 GitHub Pages<br/>Frontend React<br/><i>Toujours en ligne</i>"]
    end

    subgraph "☁️ Supabase (gratuit)"
        DB["🗄️ Base de données<br/>PostgreSQL<br/><i>Toujours en ligne</i>"]
    end

    subgraph "💻 Ton PC"
        LT["🚇 LocalTunnel<br/><i>bindx-api.loca.lt</i>"]
        DK["🐳 Docker<br/>FastAPI + Celery + Redis"]
    end

    GH -->|"API calls"| LT
    LT --> DK
    DK -->|"données"| DB

    style GH fill:#0d1117,stroke:#00e6a0,color:#fff
    style DB fill:#0d1117,stroke:#3b82f6,color:#fff
    style LT fill:#0d1117,stroke:#f59e0b,color:#fff
    style DK fill:#0d1117,stroke:#06b6d4,color:#fff
```

**Comment ça marche :**
- Le **frontend** (interface) est hébergé gratuitement sur GitHub → toujours accessible
- La **base de données** est sur Supabase → toujours accessible
- Le **backend** (calculs) tourne sur ton PC → accessible uniquement quand ton PC est allumé + tunnel lancé

**Pour démarrer :** double-clic sur "BindX Start" sur le bureau.

---

## Stack technique

| Couche | Technologie | Rôle |
|--------|------------|------|
| **Frontend** | React + Vite + Tailwind | Interface utilisateur (tableau, graphiques, 3D) |
| **State** | Zustand | Gestion de l'état côté navigateur |
| **API** | FastAPI (Python) | Serveur qui traite les requêtes |
| **Tâches** | Celery + Redis | Exécution des calculs longs en arrière-plan |
| **Base de données** | Supabase (PostgreSQL) | Stockage des projets, molécules, résultats |
| **Auth** | Supabase Auth (JWT) | Login / inscription sécurisé |
| **Chimie** | RDKit, GNINA, P2Rank | Calculs moléculaires, docking, détection de poches |
| **IA/ML** | PyTorch, admet-ai | Prédictions ADMET par intelligence artificielle |
| **Conteneurs** | Docker Compose | Emballe tout dans des boîtes isolées |
| **Hébergement** | GitHub Pages + LocalTunnel | Frontend gratuit + tunnel vers le PC |

---

*Dernière mise à jour : mars 2026*
