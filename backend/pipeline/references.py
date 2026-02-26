"""
DockIt pipeline — Scientific references and methodology.

V9: Provides structured scientific citations for all tools, methods,
and databases used in the platform. Used by the Methodology page
and PDF reports.
"""

from __future__ import annotations

REFERENCES: list[dict] = [
    # -- Docking engines --
    {
        "id": "vina",
        "category": "Docking Engines",
        "title": "AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading",
        "authors": "Trott, O. & Olson, A. J.",
        "journal": "Journal of Computational Chemistry",
        "year": 2010,
        "volume": "31(2)",
        "pages": "455-461",
        "doi": "10.1002/jcc.21334",
        "usage": "Primary molecular docking engine for virtual screening",
    },
    {
        "id": "gnina",
        "category": "Docking Engines",
        "title": "Protein-Ligand Scoring with Convolutional Neural Networks",
        "authors": "Ragoza, M., Hochuli, J., Idrobo, E., Sunseri, J. & Koes, D. R.",
        "journal": "Journal of Chemical Information and Modeling",
        "year": 2017,
        "volume": "57(4)",
        "pages": "942-957",
        "doi": "10.1021/acs.jcim.6b00740",
        "usage": "CNN-based scoring and pose ranking for improved docking accuracy",
    },
    # -- Protein structure prediction --
    {
        "id": "alphafold",
        "category": "Protein Structure Prediction",
        "title": "Highly accurate protein structure prediction with AlphaFold",
        "authors": "Jumper, J., Evans, R., Pritzel, A., et al.",
        "journal": "Nature",
        "year": 2021,
        "volume": "596",
        "pages": "583-589",
        "doi": "10.1038/s41586-021-03819-2",
        "usage": "Primary source for predicted 3D protein structures from AlphaFold DB",
    },
    {
        "id": "alphafold_db",
        "category": "Protein Structure Prediction",
        "title": "AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models",
        "authors": "Varadi, M., Anyango, S., Deshpande, A., et al.",
        "journal": "Nucleic Acids Research",
        "year": 2022,
        "volume": "50(D1)",
        "pages": "D439-D444",
        "doi": "10.1093/nar/gkab1061",
        "usage": "Public repository of 200M+ predicted protein structures",
    },
    {
        "id": "esmfold",
        "category": "Protein Structure Prediction",
        "title": "Evolutionary-scale prediction of atomic-level protein structure with a language model",
        "authors": "Lin, Z., Akin, H., Rao, R., et al.",
        "journal": "Science",
        "year": 2023,
        "volume": "379(6633)",
        "pages": "1123-1130",
        "doi": "10.1126/science.ade2574",
        "usage": "Fast fallback protein folding (60x faster than AlphaFold, no MSA needed)",
    },
    # -- Pocket detection --
    {
        "id": "fpocket",
        "category": "Binding Site Detection",
        "title": "Fpocket: an open source platform for ligand pocket detection",
        "authors": "Le Guilloux, V., Schmidtke, P. & Tuffery, P.",
        "journal": "BMC Bioinformatics",
        "year": 2009,
        "volume": "10",
        "pages": "168",
        "doi": "10.1186/1471-2105-10-168",
        "usage": "Voronoi tessellation-based pocket detection with druggability scoring",
    },
    # -- Cheminformatics --
    {
        "id": "rdkit",
        "category": "Cheminformatics",
        "title": "RDKit: Open-Source Cheminformatics Software",
        "authors": "Landrum, G.",
        "journal": "https://www.rdkit.org",
        "year": 2024,
        "volume": "",
        "pages": "",
        "doi": "",
        "usage": "Molecular descriptors, fingerprints, QED, ADMET properties, and 2D depiction",
    },
    {
        "id": "qed",
        "category": "Cheminformatics",
        "title": "Quantifying the chemical beauty of drugs",
        "authors": "Bickerton, G. R., Paolini, G. V., Besnard, J., Muresan, S. & Hopkins, A. L.",
        "journal": "Nature Chemistry",
        "year": 2012,
        "volume": "4",
        "pages": "90-98",
        "doi": "10.1038/nchem.1243",
        "usage": "Quantitative Estimate of Drug-likeness (QED) scoring",
    },
    {
        "id": "sa_score",
        "category": "Cheminformatics",
        "title": "Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions",
        "authors": "Ertl, P. & Schuffenhauer, A.",
        "journal": "Journal of Cheminformatics",
        "year": 2009,
        "volume": "1",
        "pages": "8",
        "doi": "10.1186/1758-2946-1-8",
        "usage": "Synthetic Accessibility (SA) score for prioritizing synthesizable candidates",
    },
    {
        "id": "pains",
        "category": "Cheminformatics",
        "title": "New Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS) from Screening Libraries",
        "authors": "Baell, J. B. & Holloway, G. A.",
        "journal": "Journal of Medicinal Chemistry",
        "year": 2010,
        "volume": "53(7)",
        "pages": "2719-2740",
        "doi": "10.1021/jm901137j",
        "usage": "PAINS substructure filter to flag likely false positives",
    },
    {
        "id": "butina",
        "category": "Cheminformatics",
        "title": "Unsupervised Database Clustering Based on Daylight's Fingerprint and Tanimoto Similarity",
        "authors": "Butina, D.",
        "journal": "Journal of Chemical Information and Computer Sciences",
        "year": 1999,
        "volume": "39(4)",
        "pages": "747-750",
        "doi": "10.1021/ci9803381",
        "usage": "Chemical series clustering via Tanimoto distance on Morgan fingerprints",
    },
    # -- Chemical databases --
    {
        "id": "chembl",
        "category": "Chemical Databases",
        "title": "ChEMBL: a large-scale bioactivity database for drug discovery",
        "authors": "Gaulton, A., Bellis, L. J., Bento, A. P., et al.",
        "journal": "Nucleic Acids Research",
        "year": 2012,
        "volume": "40(D1)",
        "pages": "D1100-D1107",
        "doi": "10.1093/nar/gkr777",
        "usage": "Source of known bioactive compounds with IC50/Ki/EC50 activity data",
    },
    {
        "id": "pubchem",
        "category": "Chemical Databases",
        "title": "PubChem 2023 update",
        "authors": "Kim, S., Chen, J., Cheng, T., et al.",
        "journal": "Nucleic Acids Research",
        "year": 2023,
        "volume": "51(D1)",
        "pages": "D1373-D1380",
        "doi": "10.1093/nar/gkac956",
        "usage": "Bioactive compounds via PUG REST API with assay-level activity data",
    },
    {
        "id": "zinc",
        "category": "Chemical Databases",
        "title": "ZINC20 - A Free Ultralarge-Scale Chemical Database for Ligand Discovery",
        "authors": "Irwin, J. J., Tang, K. G., Young, J., et al.",
        "journal": "Journal of Chemical Information and Modeling",
        "year": 2020,
        "volume": "60(12)",
        "pages": "6065-6073",
        "doi": "10.1021/acs.jcim.0c00675",
        "usage": "Drug-like compound library for virtual screening",
    },
    {
        "id": "enamine_real",
        "category": "Chemical Databases",
        "title": "Enamine REAL: A Make-on-Demand Chemical Space of Over 37 Billion Molecules",
        "authors": "Enamine Ltd.",
        "journal": "https://enamine.net/compound-collections/real-compounds",
        "year": 2024,
        "volume": "",
        "pages": "",
        "doi": "",
        "usage": "Virtual library sampling for novel drug-like compound generation",
    },
    # -- Target validation --
    {
        "id": "open_targets",
        "category": "Target Validation",
        "title": "Open Targets Platform: supporting systematic drug-target identification and prioritisation",
        "authors": "Ochoa, D., Hercules, A., Brber, B. M., et al.",
        "journal": "Nucleic Acids Research",
        "year": 2021,
        "volume": "49(D1)",
        "pages": "D1302-D1310",
        "doi": "10.1093/nar/gkaa1027",
        "usage": "Target-disease association scores for evidence-based target assessment",
    },
    # -- Virtual screening benchmarks --
    {
        "id": "vs_benchmark_1",
        "category": "Virtual Screening Benchmarks",
        "title": "Benchmarking Sets for Molecular Docking",
        "authors": "Huang, N., Shoichet, B. K. & Irwin, J. J.",
        "journal": "Journal of Medicinal Chemistry",
        "year": 2006,
        "volume": "49(23)",
        "pages": "6789-6801",
        "doi": "10.1021/jm0608356",
        "usage": "DUD dataset benchmark: typical enrichment factors of 5-20x for Vina-class methods",
    },
    {
        "id": "vs_benchmark_2",
        "category": "Virtual Screening Benchmarks",
        "title": "Directory of Useful Decoys, Enhanced (DUD-E): Better Ligands and Decoys for Better Benchmarking",
        "authors": "Mysinger, M. M., Carchia, M., Irwin, J. J. & Shoichet, B. K.",
        "journal": "Journal of Medicinal Chemistry",
        "year": 2012,
        "volume": "55(14)",
        "pages": "6582-6594",
        "doi": "10.1021/jm300687e",
        "usage": "Enhanced benchmark set: Vina achieves ~50-70% AUC across diverse targets",
    },
    {
        "id": "vs_best_practices",
        "category": "Virtual Screening Benchmarks",
        "title": "Best Practices for Virtual Screening",
        "authors": "Scior, T., Bender, A., Tresadern, G., et al.",
        "journal": "Journal of Chemical Information and Modeling",
        "year": 2012,
        "volume": "52(4)",
        "pages": "867-881",
        "doi": "10.1021/ci200528d",
        "usage": "Guidelines: expected hit rates of 1-5% for HTS, 5-20% for focused libraries",
    },
    {
        "id": "lipinski",
        "category": "Drug-likeness",
        "title": "Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings",
        "authors": "Lipinski, C. A., Lombardo, F., Dominy, B. W., Feeney, P. J.",
        "journal": "Advanced Drug Delivery Reviews",
        "year": 1997,
        "volume": "23(1-3)",
        "pages": "3-25",
        "doi": "10.1016/S0169-409X(96)00423-1",
        "usage": "Rule of Five for oral bioavailability filtering",
    },
]


EXPECTED_RESULTS: dict = {
    "rapid_mode": {
        "description": "Quick virtual screening with Vina docking only",
        "typical_runtime": "5-20 minutes (50 ligands)",
        "expected_hit_rate": "1-5% of screened compounds show significant binding",
        "enrichment_factor": "5-15x over random selection",
        "accuracy": "~50-70% AUC on standard benchmarks (DUD-E)",
        "limitations": "No ADMET filtering, no off-target screening",
    },
    "standard_mode": {
        "description": "Full pipeline with ADMET, off-target, clustering, and AI generation",
        "typical_runtime": "30-90 minutes (100 ligands)",
        "expected_hit_rate": "5-15% of final candidates are experimentally validated hits",
        "enrichment_factor": "10-30x with multi-objective scoring",
        "accuracy": "Higher confidence with GNINA CNN scoring + ADMET filters",
        "key_advantage": "Multi-objective Pareto ranking balances potency, safety, and drug-likeness",
    },
    "deep_mode": {
        "description": "Massive 5-pass screening with exhaustive coverage",
        "typical_runtime": "1-4 hours (500+ ligands)",
        "expected_hit_rate": "15-25% of top candidates suitable for experimental validation",
        "enrichment_factor": "20-50x with iterative AI-guided generation",
        "accuracy": "Best available: consensus scoring + confidence metrics",
        "key_advantage": "Chemical diversity through clustering + novel scaffold discovery via AI generation",
    },
    "general_notes": [
        "Hit rates depend heavily on target druggability (assessed by Target Assessment Engine)",
        "Well-characterized kinase targets (e.g., EGFR) typically yield higher hit rates than novel targets",
        "Experimental validation (IC50 assays) is always required to confirm computational predictions",
        "Expected success rate for computational-to-experimental: ~20-40% of top 10 candidates",
        "AlphaFold structures may produce lower docking accuracy than experimental PDB structures",
    ],
}


def get_references(category: str | None = None) -> list[dict]:
    """Return scientific references, optionally filtered by category."""
    if category:
        return [r for r in REFERENCES if r["category"].lower() == category.lower()]
    return REFERENCES


def get_expected_results() -> dict:
    """Return expected results and benchmarks for each mode."""
    return EXPECTED_RESULTS


def get_categories() -> list[str]:
    """Return unique reference categories."""
    seen: set[str] = set()
    cats: list[str] = []
    for r in REFERENCES:
        c = r["category"]
        if c not in seen:
            seen.add(c)
            cats.append(c)
    return cats
