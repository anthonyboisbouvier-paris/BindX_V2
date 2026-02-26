"""
Comprehensive kinase inhibitor benchmark dataset for DockIt validation.

Contains well-characterized inhibitors with published IC50/Ki values for 5 kinase targets,
plus common decoy compounds (non-kinase drugs) for selectivity assessment.

All SMILES retrieved from PubChem (canonical form). IC50 values from published literature
and vendor datasheets (Selleck, MedChemExpress, primary publications).

IC50 values are from cell-free biochemical kinase assays unless otherwise noted.
For covalent KRAS G12C inhibitors, cellular IC50 values are reported (as biochemical
IC50 is not meaningful for covalent-irreversible inhibitors of a non-enzyme target).

Sources:
    - PubChem REST API (https://pubchem.ncbi.nlm.nih.gov)
    - ChEMBL (https://www.ebi.ac.uk/chembl)
    - Selleck Chemicals datasheets
    - MedChemExpress datasheets
    - Primary literature (see per-compound source fields)
"""

from typing import Any

# ---------------------------------------------------------------------------
# 1. EGFR (P00533) -- Epidermal Growth Factor Receptor
# ---------------------------------------------------------------------------
# IC50 values from cell-free kinase assays against wild-type EGFR unless noted.
# References: Selleck, MedChemExpress, Moyer et al. 1997, Li et al. 2008,
# Engelman et al. 2007, Cross et al. 2014 (PMC4770737)

EGFR_ACTIVES: list[dict[str, Any]] = [
    {
        "name": "Erlotinib",
        "smiles": "COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC",
        "ic50_nM": 2.0,
        "source": "CHEMBL553 / Moyer et al. Cancer Res 1997",
        "notes": "1st-gen reversible TKI, cell-free IC50",
    },
    {
        "name": "Gefitinib",
        "smiles": "COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4",
        "ic50_nM": 33.0,
        "source": "CHEMBL939 / Wakeling et al. Cancer Res 2002",
        "notes": "1st-gen reversible TKI, cell-free kinase IC50",
    },
    {
        "name": "Lapatinib",
        "smiles": "CS(=O)(=O)CCNCC1=CC=C(O1)C2=CC3=C(C=C2)N=CN=C3NC4=CC(=C(C=C4)OCC5=CC(=CC=C5)F)Cl",
        "ic50_nM": 10.8,
        "source": "CHEMBL554 / Rusnak et al. Mol Cancer Ther 2001",
        "notes": "Dual EGFR/HER2 reversible inhibitor",
    },
    {
        "name": "Afatinib",
        "smiles": "CN(C)C/C=C/C(=O)NC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC(=C(C=C3)F)Cl)O[C@H]4CCOC4",
        "ic50_nM": 0.5,
        "source": "CHEMBL1173655 / Li et al. Oncogene 2008",
        "notes": "2nd-gen irreversible TKI, IC50 EGFR wt",
    },
    {
        "name": "Osimertinib",
        "smiles": "CN1C=C(C2=CC=CC=C21)C3=NC(=NC=C3)NC4=C(C=C(C(=C4)NC(=O)C=C)N(C)CCN(C)C)OC",
        "ic50_nM": 12.0,
        "source": "CHEMBL3353410 / Cross et al. Cancer Discov 2014",
        "notes": "3rd-gen, mutant-selective (T790M), cell-free IC50 wt EGFR",
    },
    {
        "name": "Icotinib",
        "smiles": "C#CC1=CC(=CC=C1)NC2=NC=NC3=CC4=C(C=C32)OCCOCCOCCO4",
        "ic50_nM": 5.0,
        "source": "CHEMBL1289930 / Tan et al. Lung Cancer 2012",
        "notes": "1st-gen Chinese-developed TKI",
    },
    {
        "name": "Neratinib",
        "smiles": "CCOC1=C(C=C2C(=C1)N=CC(=C2NC3=CC(=C(C=C3)OCC4=CC=CC=N4)Cl)C#N)NC(=O)/C=C/CN(C)C",
        "ic50_nM": 92.0,
        "source": "CHEMBL1289601 / Rabindran et al. Cancer Res 2004",
        "notes": "Pan-HER irreversible inhibitor, EGFR IC50 cell-free",
    },
    {
        "name": "Dacomitinib",
        "smiles": "COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)NC(=O)/C=C/CN4CCCCC4",
        "ic50_nM": 6.0,
        "source": "CHEMBL2110732 / Engelman et al. Cancer Res 2007",
        "notes": "2nd-gen irreversible pan-HER TKI",
    },
    {
        "name": "Vandetanib",
        "smiles": "CN1CCC(CC1)COC2=C(C=C3C(=C2)N=CN=C3NC4=C(C=C(C=C4)Br)F)OC",
        "ic50_nM": 500.0,
        "source": "CHEMBL24828 / Wedge et al. Cancer Res 2002",
        "notes": "Multi-target: EGFR/VEGFR2/RET",
    },
    {
        "name": "Canertinib",
        "smiles": "C=CC(=O)NC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4",
        "ic50_nM": 1.5,
        "source": "CHEMBL136876 / Smaill et al. J Med Chem 2000",
        "notes": "Pan-ErbB irreversible inhibitor (CI-1033)",
    },
    {
        "name": "Pelitinib",
        "smiles": "CCOC1=C(C=C2C(=C1)N=CC(=C2NC3=CC(=C(C=C3)F)Cl)C#N)NC(=O)/C=C/CN(C)C",
        "ic50_nM": 39.0,
        "source": "CHEMBL197279 / Torrance et al. Nat Med 2001",
        "notes": "Irreversible EGFR inhibitor (EKB-569)",
    },
    {
        "name": "AEE788",
        "smiles": "CCN1CCN(CC1)CC2=CC=C(C=C2)C3=CC4=C(N3)N=CN=C4N[C@H](C)C5=CC=CC=C5",
        "ic50_nM": 2.0,
        "source": "CHEMBL477772 / Traxler et al. Cancer Res 2004",
        "notes": "EGFR/VEGFR dual inhibitor (NVP-AEE788)",
    },
    {
        "name": "WZ4002",
        "smiles": "CN1CCN(CC1)C2=CC(=C(C=C2)NC3=NC=C(C(=N3)OC4=CC=CC(=C4)NC(=O)C=C)Cl)OC",
        "ic50_nM": 2.0,
        "source": "CHEMBL2023995 / Zhou et al. Nature 2009",
        "notes": "Mutant-selective 3rd-gen tool compound (T790M)",
    },
]

# ---------------------------------------------------------------------------
# 2. CDK2 (P24941) -- Cyclin-Dependent Kinase 2
# ---------------------------------------------------------------------------
# IC50 values from cell-free kinase assays (CDK2/CycA or CDK2/CycE complexes)

CDK2_ACTIVES: list[dict[str, Any]] = [
    {
        "name": "Dinaciclib",
        "smiles": "CCC1=C2N=C(C=C(N2N=C1)NCC3=C[N+](=CC=C3)[O-])N4CCCC[C@H]4CCO",
        "ic50_nM": 1.0,
        "source": "CHEMBL515634 / Parry et al. Mol Cancer Ther 2010",
        "notes": "CDK2/CycE IC50=1 nM; also CDK1(3nM), CDK5(1nM), CDK9(4nM)",
    },
    {
        "name": "Palbociclib",
        "smiles": "CC1=C(C(=O)N(C2=NC(=NC=C12)NC3=NC=C(C=C3)N4CCNCC4)C5CCCC5)C(=O)C",
        "ic50_nM": 11.0,
        "source": "CHEMBL1229517 / Fry et al. Mol Cancer Ther 2004",
        "notes": "Primarily CDK4/6 selective; CDK2 IC50 much higher",
    },
    {
        "name": "Ribociclib",
        "smiles": "CN(C)C(=O)C1=CC2=CN=C(N=C2N1C3CCCC3)NC4=NC=C(C=C4)N5CCNCC5",
        "ic50_nM": 39.0,
        "source": "CHEMBL3545110 / Hortobagyi et al. NEJM 2016",
        "notes": "CDK4/6 selective; weak CDK2 activity",
    },
    {
        "name": "Flavopiridol",
        "smiles": "CN1CC[C@@H]([C@@H](C1)O)C2=C(C=C(C3=C2OC(=CC3=O)C4=CC=CC=C4Cl)O)O",
        "ic50_nM": 12.0,
        "source": "CHEMBL44551 / Senderowicz et al. J Natl Cancer Inst 1998",
        "notes": "Pan-CDK 1st-gen; CDK2/CycA IC50=12 nM",
    },
    {
        "name": "Roscovitine",
        "smiles": "CC[C@H](CO)NC1=NC(=C2C(=N1)N(C=N2)C(C)C)NCC3=CC=CC=C3",
        "ic50_nM": 700.0,
        "source": "CHEMBL32485 / Meijer et al. Eur J Biochem 1997",
        "notes": "CDK2/CycE IC50; also known as Seliciclib / CYC202",
    },
    {
        "name": "SNS-032",
        "smiles": "CC(C)(C)C1=CN=C(O1)CSC2=CN=C(S2)NC(=O)C3CCNCC3",
        "ic50_nM": 48.0,
        "source": "CHEMBL477781 / Conroy et al. Anticancer Drugs 2009",
        "notes": "CDK2/CycE IC50=48 nM; also known as BMS-387032",
    },
    {
        "name": "AT7519",
        "smiles": "C1CNCCC1NC(=O)C2=C(C=NN2)NC(=O)C3=C(C=CC=C3Cl)Cl",
        "ic50_nM": 47.0,
        "source": "CHEMBL510068 / Wyatt et al. Bioorg Med Chem Lett 2008",
        "notes": "Multi-CDK inhibitor; CDK2 IC50=47 nM",
    },
    {
        "name": "AZD5438",
        "smiles": "CC1=NC=C(N1C(C)C)C2=NC(=NC=C2)NC3=CC=C(C=C3)S(=O)(=O)C",
        "ic50_nM": 6.0,
        "source": "CHEMBL575417 / Byth et al. J Med Chem 2009",
        "notes": "CDK1/2/9 inhibitor; CDK2 IC50=6 nM (RCSB:6GUE)",
    },
    {
        "name": "Milciclib",
        "smiles": "CC1(CC2=CN=C(N=C2C3=C1C(=NN3C)C(=O)NC)NC4=CC=C(C=C4)N5CCN(CC5)C)C",
        "ic50_nM": 45.0,
        "source": "CHEMBL496880 / Schiewer et al. Cancer Res 2012",
        "notes": "CDK2/CycA IC50=45 nM; also PHA-848125",
    },
    {
        "name": "Roniciclib",
        "smiles": "C[C@H]([C@@H](C)OC1=NC(=NC=C1C(F)(F)F)NC2=CC=C(C=C2)[S@](=N)(=O)C3CC3)O",
        "ic50_nM": 7.0,
        "source": "CHEMBL3301607 / Seki et al. ACS Chem Biol 2016",
        "notes": "Pan-CDK; CDK2 IC50=7 nM; also BAY-1000394",
    },
    {
        "name": "PHA-793887",
        "smiles": "CC(C)CC(=O)NC1=NNC2=C1CN(C2(C)C)C(=O)C3CCN(CC3)C",
        "ic50_nM": 8.0,
        "source": "CHEMBL590744 / Brasca et al. J Med Chem 2009",
        "notes": "Pan-CDK; CDK2 IC50=8 nM",
    },
]

# ---------------------------------------------------------------------------
# 3. BRAF V600E (P15056) -- B-Raf Proto-oncogene (V600E mutant)
# ---------------------------------------------------------------------------
# IC50 values from cell-free kinase assays against BRAF V600E mutant

BRAF_ACTIVES: list[dict[str, Any]] = [
    {
        "name": "Vemurafenib",
        "smiles": "CCCS(=O)(=O)NC1=C(C(=C(C=C1)F)C(=O)C2=CNC3=C2C=C(C=N3)C4=CC=C(C=C4)Cl)F",
        "ic50_nM": 31.0,
        "source": "CHEMBL1229211 / Bollag et al. Nature 2010",
        "notes": "V600E-selective; PLX4032",
    },
    {
        "name": "Dabrafenib",
        "smiles": "CC(C)(C)C1=NC(=C(S1)C2=NC(=NC=C2)N)C3=C(C(=CC=C3)NS(=O)(=O)C4=C(C=CC=C4F)F)F",
        "ic50_nM": 0.65,
        "source": "CHEMBL2028663 / Laquerre et al. AACR 2009",
        "notes": "V600E IC50=0.65 nM, GSK2118436",
    },
    {
        "name": "Encorafenib",
        "smiles": "C[C@@H](CNC1=NC=CC(=N1)C2=CN(N=C2C3=C(C(=CC(=C3)Cl)NS(=O)(=O)C)F)C(C)C)NC(=O)OC",
        "ic50_nM": 0.35,
        "source": "CHEMBL3137309 / Stuart et al. Cancer Discov 2012",
        "notes": "V600E IC50=0.35 nM; longest diss. t1/2 >30h; LGX818",
    },
    {
        "name": "Sorafenib",
        "smiles": "CNC(=O)C1=NC=CC(=C1)OC2=CC=C(C=C2)NC(=O)NC3=CC(=C(C=C3)Cl)C(F)(F)F",
        "ic50_nM": 38.0,
        "source": "CHEMBL1336 / Wilhelm et al. Cancer Res 2004",
        "notes": "Multi-kinase; BRAF-V600E IC50=38 nM; CRAF IC50=6 nM",
    },
    {
        "name": "PLX4720",
        "smiles": "CCCS(=O)(=O)NC1=C(C(=C(C=C1)F)C(=O)C2=CNC3=C2C=C(C=N3)Cl)F",
        "ic50_nM": 13.0,
        "source": "CHEMBL1230011 / Tsai et al. PNAS 2008",
        "notes": "V600E-selective tool compound; precursor to vemurafenib",
    },
    {
        "name": "TAK-632",
        "smiles": "C1CC1C(=O)NC2=NC3=C(S2)C(=C(C=C3)OC4=CC(=C(C=C4)F)NC(=O)CC5=CC(=CC=C5)C(F)(F)F)C#N",
        "ic50_nM": 2.4,
        "source": "CHEMBL3310823 / Nakamura et al. Cancer Lett 2013",
        "notes": "Pan-RAF type II; V600E IC50=2.4 nM; CRAF IC50=1.4 nM",
    },
    {
        "name": "LY3009120",
        "smiles": "CC1=CC(=C(C=C1C2=C(N=C3C(=C2)C=NC(=N3)NC)C)NC(=O)NCCC(C)(C)C)F",
        "ic50_nM": 5.8,
        "source": "CHEMBL3590575 / Henry et al. J Med Chem 2015",
        "notes": "Pan-RAF dimer inhibitor; BRAF-V600E IC50=5.8 nM",
    },
    {
        "name": "AZ628",
        "smiles": "CC1=C(C=C(C=C1)NC(=O)C2=CC(=CC=C2)C(C)(C)C#N)NC3=CC4=C(C=C3)N=CN(C4=O)C",
        "ic50_nM": 34.0,
        "source": "CHEMBL2032433 / Hatzivassiliou et al. Nature 2010",
        "notes": "ATP-competitive RAF inhibitor; V600E IC50=34 nM",
    },
    {
        "name": "GDC-0879",
        "smiles": "C1C/C(=N\\O)/C2=C1C=C(C=C2)C3=CN(N=C3C4=CC=NC=C4)CCO",
        "ic50_nM": 0.13,
        "source": "CHEMBL1213492 / Hoeflich et al. Cancer Res 2009",
        "notes": "V600E enzyme IC50=0.13 nM; cell pERK IC50=63 nM",
    },
    {
        "name": "SB590885",
        "smiles": "CN(C)CCOC1=CC=C(C=C1)C2=NC(=C(N2)C3=CC=NC=C3)C4=CC5=C(C=C4)C(=NO)CC5",
        "ic50_nM": 0.16,
        "source": "CHEMBL372310 / King et al. Cancer Res 2006",
        "notes": "Ki=0.16 nM B-Raf; 11x selective over c-Raf",
    },
    {
        "name": "RAF265",
        "smiles": "CN1C2=C(C=C(C=C2)OC3=CC(=NC=C3)C4=NC=C(N4)C(F)(F)F)N=C1NC5=CC=C(C=C5)C(F)(F)F",
        "ic50_nM": 30.0,
        "source": "CHEMBL2001430 / Stuart et al. Cancer Res 2008",
        "notes": "Multi-kinase RAF/VEGFR; BRAF IC50=3-60 nM (CHIR-265)",
    },
]

# ---------------------------------------------------------------------------
# 4. JAK2 (P52333) -- Janus Kinase 2
# ---------------------------------------------------------------------------
# IC50 values from cell-free enzymatic kinase assays against JAK2

JAK2_ACTIVES: list[dict[str, Any]] = [
    {
        "name": "Ruxolitinib",
        "smiles": "C1CCC(C1)[C@@H](CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3",
        "ic50_nM": 2.8,
        "source": "CHEMBL1789941 / Quintas-Cardama et al. Blood 2010",
        "notes": "JAK1/2 selective; JAK2 IC50=2.8 nM (cell-free)",
    },
    {
        "name": "Tofacitinib",
        "smiles": "C[C@@H]1CCN(C[C@@H]1N(C)C2=NC=NC3=C2C=CN3)C(=O)CC#N",
        "ic50_nM": 20.0,
        "source": "CHEMBL221959 / Flanagan et al. J Med Chem 2010",
        "notes": "Pan-JAK; JAK2 IC50 ~20 nM (primary JAK3 IC50=1 nM)",
    },
    {
        "name": "Baricitinib",
        "smiles": "CCS(=O)(=O)N1CC(C1)(CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3",
        "ic50_nM": 5.7,
        "source": "CHEMBL2105735 / Fridman et al. J Immunol 2010",
        "notes": "JAK1/2 selective; JAK2 IC50=5.7 nM",
    },
    {
        "name": "Fedratinib",
        "smiles": "CC1=CN=C(N=C1NC2=CC(=CC=C2)S(=O)(=O)NC(C)(C)C)NC3=CC=C(C=C3)OCCN4CCCC4",
        "ic50_nM": 3.0,
        "source": "CHEMBL1229896 / Wernig et al. Cancer Cell 2008",
        "notes": "JAK2-selective; IC50=3 nM (also known as TG101348/SAR302503)",
    },
    {
        "name": "Pacritinib",
        "smiles": "C1CCN(C1)CCOC2=C3COC/C=C/COCC4=CC(=CC=C4)C5=NC(=NC=C5)NC(=C3)C=C2",
        "ic50_nM": 23.0,
        "source": "CHEMBL2170845 / Hart et al. Leukemia 2011",
        "notes": "JAK2/FLT3 dual inhibitor; JAK2 IC50=23 nM",
    },
    {
        "name": "Momelotinib",
        "smiles": "C1COCCN1C2=CC=C(C=C2)NC3=NC=CC(=N3)C4=CC=C(C=C4)C(=O)NCC#N",
        "ic50_nM": 18.0,
        "source": "CHEMBL2388354 / Pardanani et al. Leukemia 2009",
        "notes": "JAK1/2/ACVR1 inhibitor; also CYT387",
    },
    {
        "name": "Gandotinib",
        "smiles": "CC1=CC(=NN1)NC2=NN3C(=C(N=C3C(=C2)CN4CCOCC4)C)CC5=C(C=C(C=C5)Cl)F",
        "ic50_nM": 3.0,
        "source": "CHEMBL2007613 / Ma et al. Blood 2013",
        "notes": "JAK2-V617F selective; also LY2784544",
    },
    {
        "name": "Cerdulatinib",
        "smiles": "CCS(=O)(=O)N1CCN(CC1)C2=CC=C(C=C2)NC3=NC=C(C(=N3)NC4CC4)C(=O)N",
        "ic50_nM": 6.0,
        "source": "CHEMBL3331629 / Ma et al. J Pharmacol Exp Ther 2015",
        "notes": "Dual SYK/JAK; JAK2 IC50=6 nM (cell-free)",
    },
    {
        "name": "Filgotinib",
        "smiles": "C1CC1C(=O)NC2=NN3C(=N2)C=CC=C3C4=CC=C(C=C4)CN5CCS(=O)(=O)CC5",
        "ic50_nM": 28.0,
        "source": "CHEMBL3301607 / Van Rompaey et al. J Immunol 2013",
        "notes": "JAK1-selective; JAK2 IC50=28 nM (GLPG0634)",
    },
    {
        "name": "Itacitinib",
        "smiles": "C1CN(CCC1N2CC(C2)(CC#N)N3C=C(C=N3)C4=C5C=CNC5=NC=N4)C(=O)C6=C(C(=NC=C6)C(F)(F)F)F",
        "ic50_nM": 63.0,
        "source": "CHEMBL3545265 / Bissonnette et al. Br J Clin Pharmacol 2020",
        "notes": "JAK1-selective; JAK2 IC50=63 nM (INCB039110)",
    },
    {
        "name": "AZD1480",
        "smiles": "CC1=CC(=NN1)NC2=NC(=NC=C2Cl)N[C@@H](C)C3=NC=C(C=N3)F",
        "ic50_nM": 0.26,
        "source": "CHEMBL1229517 / Derenzini et al. Oncotarget 2011",
        "notes": "JAK2-selective; IC50=0.26 nM (enzymatic)",
    },
    {
        "name": "TG101348",
        "smiles": "CC1=CN=C(N=C1NC2=CC(=CC=C2)S(=O)(=O)NC(C)(C)C)NC3=CC=C(C=C3)OCCN4CCCC4",
        "ic50_nM": 3.0,
        "source": "CHEMBL1229896 / Wernig et al. Cancer Cell 2008",
        "notes": "Same as Fedratinib (SAR302503). Included for alternate name.",
    },
]

# ---------------------------------------------------------------------------
# 5. KRAS G12C (P01116) -- KRAS GTPase, G12C mutant
# ---------------------------------------------------------------------------
# Note: KRAS G12C is NOT a kinase. These are covalent inhibitors that bind
# the GDP-bound inactive form. IC50 values are from cellular assays (pERK
# inhibition or KRAS-GTP reduction) since biochemical IC50 for covalent
# inhibitors targeting a non-enzymatic pocket is assay-dependent.

KRAS_ACTIVES: list[dict[str, Any]] = [
    {
        "name": "Sotorasib",
        "smiles": "C[C@H]1CN(CCN1C2=NC(=O)N(C3=NC(=C(C=C32)F)C4=C(C=CC=C4F)O)C5=C(C=CN=C5C(C)C)C)C(=O)C=C",
        "ic50_nM": 35.0,
        "source": "CHEMBL4523824 / Lanman et al. J Med Chem 2020",
        "notes": "AMG-510; 1st FDA-approved KRAS G12C inhibitor; cellular KRAS-GTP IC50",
    },
    {
        "name": "Adagrasib",
        "smiles": "CN1CCC[C@H]1COC2=NC3=C(CCN(C3)C4=CC=CC5=C4C(=CC=C5)Cl)C(=N2)N6CCN([C@H](C6)CC#N)C(=O)C(=C)F",
        "ic50_nM": 78.0,
        "source": "CHEMBL4594389 / Fell et al. J Med Chem 2020",
        "notes": "MRTX849; 2nd approved; cellular KRAS-GTP IC50",
    },
    {
        "name": "Divarasib",
        "smiles": "C[C@H]1CN(CCN1C2=NC(=NC3=C(C(=C(C=C32)Cl)C4=C(C(=CC(=N4)N)C)C(F)(F)F)F)OC[C@@H]5CCCN5C)C(=O)C=C",
        "ic50_nM": 0.32,
        "source": "CHEMBL4802955 / Sacher et al. NEJM 2023",
        "notes": "GDC-6036; next-gen; KRAS-alkylation IC50=0.32 nM (NCI-HCC1171)",
    },
    {
        "name": "ARS-1620",
        "smiles": "C=CC(=O)N1CCN(CC1)C2=NC=NC3=C(C(=C(C=C32)Cl)C4=C(C=CC=C4F)O)F",
        "ic50_nM": 120.0,
        "source": "CHEMBL4078668 / Janes et al. Cell 2018",
        "notes": "1st in vivo proof-of-concept tool compound; H358 cell pERK IC50",
    },
    {
        "name": "ARS-853",
        "smiles": "CC1(CC1)C2=CC(=C(C=C2Cl)O)NCC(=O)N3CCN(CC3)C4CN(C4)C(=O)C=C",
        "ic50_nM": 1700.0,
        "source": "CHEMBL4072803 / Patricelli et al. Cancer Discov 2016",
        "notes": "Early tool compound; cell pERK IC50; poor bioavailability",
    },
    {
        "name": "MRTX1257",
        "smiles": "CC1=C2C(=CC=C1)C=CC=C2N3CCC4=C(C3)N=C(N=C4N5CCN([C@H](C5)CC#N)C(=O)C=C)OC[C@@H]6CCCN6C",
        "ic50_nM": 0.9,
        "source": "CHEMBL4594388 / Fell et al. MCR 2020",
        "notes": "Preclinical tool compound from Mirati; pERK IC50=0.9 nM (H358)",
    },
    {
        "name": "BI-2852",
        "smiles": "CN1C=C(N=C1)CN2C=CC3=C2C=C(C=C3)CNCC4=C(C5=CC=CC=C5N4)[C@@H]6C7=C(C=CC(=C7)O)C(=O)N6",
        "ic50_nM": 490.0,
        "source": "CHEMBL4563722 / Kessler et al. PNAS 2019",
        "notes": "Non-covalent switch I/II pocket binder; KRAS-G12D SOS1 displacement IC50",
    },
    {
        "name": "JNJ-74699157",
        "smiles": "C=CC(=O)N1CCN(CC1)C2=NC=NC3=C(C(=C(C=C32)Cl)C4=CC=CC(=C4)OC)F",
        "ic50_nM": 12.4,
        "source": "Lorthiois et al. J Med Chem 2022 (ARS-3248)",
        "notes": "ARS-3248; covalent G12C inhib; Ba/F3 cellular IC50",
    },
    {
        "name": "D3S-001",
        "smiles": "C=CC(=O)N1CCN(CC1)C2=NC=NC3=C(C(=C(C=C32)F)C4=CC(=CC=C4)NC(=O)C5=CC=NO5)F",
        "ic50_nM": 0.6,
        "source": "Zhu et al. Cancer Discov 2024",
        "notes": "Rapid target engagement; KRAS-GTP IC50=0.6 nM (cellular)",
    },
]

# ---------------------------------------------------------------------------
# 6. DECOYS -- Common non-kinase drugs (negative controls)
# ---------------------------------------------------------------------------
# These drugs have no known activity against any of the 5 targets above.
# SMILES from PubChem canonical representation.

DECOYS: list[dict[str, Any]] = [
    {
        "name": "Ibuprofen",
        "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "target": "COX-1/2 (NSAID)",
    },
    {
        "name": "Metformin",
        "smiles": "CN(C)C(=N)N=C(N)N",
        "target": "AMPK/complex I (anti-diabetic)",
    },
    {
        "name": "Atorvastatin",
        "smiles": "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
        "target": "HMG-CoA reductase (statin)",
    },
    {
        "name": "Omeprazole",
        "smiles": "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC",
        "target": "H+/K+ ATPase (PPI)",
    },
    {
        "name": "Amoxicillin",
        "smiles": "CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C",
        "target": "Penicillin-binding proteins (beta-lactam antibiotic)",
    },
    {
        "name": "Lisinopril",
        "smiles": "C1CC(N(C1)C(=O)C(CCCCN)NC(CCC2=CC=CC=C2)C(=O)O)C(=O)O",
        "target": "ACE (antihypertensive)",
    },
    {
        "name": "Loratadine",
        "smiles": "CCOC(=O)N1CCC(=C2C3=C(CCC4=C2N=CC=C4)C=C(C=C3)Cl)CC1",
        "target": "H1 histamine receptor (antihistamine)",
    },
    {
        "name": "Metoprolol",
        "smiles": "CC(C)NCC(COC1=CC=C(C=C1)CCOC)O",
        "target": "Beta-1 adrenergic receptor (beta-blocker)",
    },
    {
        "name": "Fluoxetine",
        "smiles": "CNCCC(C1=CC=CC=C1)OC2=CC=C(C=C2)C(F)(F)F",
        "target": "SERT (SSRI antidepressant)",
    },
    {
        "name": "Amlodipine",
        "smiles": "CCOC(=O)C1=C(NC(=C(C1C2=CC=CC=C2Cl)C(=O)OC)C)COCCN",
        "target": "L-type Ca2+ channel (CCB antihypertensive)",
    },
]

# ---------------------------------------------------------------------------
# Master dictionary for programmatic access
# ---------------------------------------------------------------------------

TARGETS: dict[str, dict[str, Any]] = {
    "EGFR": {
        "uniprot": "P00533",
        "gene": "EGFR",
        "pdb_reference": "1M17",
        "actives": EGFR_ACTIVES,
    },
    "CDK2": {
        "uniprot": "P24941",
        "gene": "CDK2",
        "pdb_reference": "1HCK",
        "actives": CDK2_ACTIVES,
    },
    "BRAF_V600E": {
        "uniprot": "P15056",
        "gene": "BRAF",
        "mutation": "V600E",
        "pdb_reference": "3OG7",
        "actives": BRAF_ACTIVES,
    },
    "JAK2": {
        "uniprot": "P52333",
        "gene": "JAK2",
        "pdb_reference": "3FUP",
        "actives": JAK2_ACTIVES,
    },
    "KRAS_G12C": {
        "uniprot": "P01116",
        "gene": "KRAS",
        "mutation": "G12C",
        "pdb_reference": "6OIM",
        "actives": KRAS_ACTIVES,
        "note": "Covalent-irreversible inhibitors; IC50 from cellular assays",
    },
}


def get_all_actives() -> list[dict[str, Any]]:
    """Return flat list of all active compounds across all targets."""
    result = []
    for target_name, target_data in TARGETS.items():
        for mol in target_data["actives"]:
            result.append({**mol, "target": target_name})
    return result


def get_benchmark_smiles(target: str) -> list[str]:
    """Return list of SMILES for a given target name."""
    if target not in TARGETS:
        raise ValueError(f"Unknown target '{target}'. Choose from: {list(TARGETS.keys())}")
    return [mol["smiles"] for mol in TARGETS[target]["actives"]]


def get_enrichment_set(target: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Return (actives, decoys) pair for enrichment factor calculation."""
    if target not in TARGETS:
        raise ValueError(f"Unknown target '{target}'. Choose from: {list(TARGETS.keys())}")
    return TARGETS[target]["actives"], DECOYS


def summary() -> None:
    """Print summary statistics for the benchmark dataset."""
    print("=" * 70)
    print("DockIt Kinase Inhibitor Benchmark Dataset")
    print("=" * 70)
    for name, data in TARGETS.items():
        actives = data["actives"]
        ic50s = [m["ic50_nM"] for m in actives]
        print(
            f"  {name:15s} ({data['uniprot']}) : "
            f"{len(actives):2d} actives, "
            f"IC50 range: {min(ic50s):.2f} - {max(ic50s):.1f} nM"
        )
    print(f"  {'DECOYS':15s}             : {len(DECOYS):2d} compounds")
    total = sum(len(d["actives"]) for d in TARGETS.values()) + len(DECOYS)
    print(f"\n  Total compounds: {total}")
    print("=" * 70)


if __name__ == "__main__":
    summary()
