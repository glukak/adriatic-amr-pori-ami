#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Adriatic coastal One Health AMR workflow - single-file pipeline
===============================================================

Input : one Excel workbook (default: adriaticlast.xlsx) with 36 coastal stations
        nested in 17 watersheds (sheets: Physical data, Chemical microbiological
        data, Animal Abundances, Antibiotics, args, args data, Bacterial community
        (16S), bacteria data).
Output: <OUT_DIR>/AMR_results.xlsx, <OUT_DIR>/figures/*.png (incl. Fig0 PORI scheme, Fig3c PORI table,
        Fig7_AMI_table, supplementary AMI tables), <OUT_DIR>/run_log.txt

Workflow (each step maps to one Methods subsection)
---------------------------------------------------
  Step 1  Exploratory correlations between antibiotic consumption and resistome
  Step 1b Gene-level test (BH-FDR, basin permutations): multidrug / cross-resistance evidence
  Step 2  Distribution and outlier screening (justifies rank-based, robust methods)
  Step 3  PORI  - Pathogen-and-resistome One-health Risk Index (descriptive)
  Step 4  Sentinel taxa and origin-based (W / E / U) source signals
  Step 5  AMI   - Antimicrobial Mitigation Index (main, step6s): PORI points avoided if the anthropogenic-taxa signal
                  falls to the level of the cleaner stations (max cut 50%), elasticities estimated from the data
                  (log-log OLS with salinity, basin-cluster bootstrap), where and how to act per station.
  Supplementary (conditional, assumed effect sizes; not in the main text): step6 = consumption-based AMI (total / critical,
                  composite dPORI); 6b = multi-class gene assignment sensitivity; 6c = multidrug AMI;

Design principles
-----------------
* PORI and AMI are DESCRIPTIVE / RELATIVE indices: scores are scaled 0-100 over the
  36 stations of this campaign and are NOT absolute risk or efficacy estimates.
* No predictive model is fitted. Basin-level predictors (17 independent units) are
  too few for a validated machine-learning model; all inference uses rank-based
  statistics with Benjamini-Hochberg false-discovery-rate control.
* Station and basin levels are kept apart: consumption and demographic data exist
  at basin level only, so basin-level tests (n = 17) are the primary inference.

Usage:  AMR_XLSX=adriaticlast.xlsx AMR_OUT=results python adriatic_amr_pipeline.py
"""

import os
import re
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata, spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

warnings.filterwarnings("ignore")

# ======================================================================
# CONFIGURATION
# ======================================================================
INPUT_XLSX = os.environ.get("AMR_XLSX", "adriaticlast.xlsx")
OUT_DIR = os.environ.get("AMR_OUT", "results")
SEED = 42
N_STATIONS = 36

EPS = 0.1                 # offset inside the geometric means of the PORI (avoids zeros)
N_MC = 3000               # Monte-Carlo draws for the PORI sensitivity analysis
MIN_PREVALENCE = 0.50     # sentinel candidates must occur in >= 50% of stations
FDR_Q_SENTINEL = 0.05     # FDR threshold for sentinel candidates
SENT_MIN_LOBO = 0.30      # sentinel must keep rho > 0.30 when any one basin is dropped
FDR_Q_GATE = 0.10         # evidence gate: class-matched consumption-ARG link (basin level)
GENE_MIN_BASINS = 5       # step 1b: a gene is tested if detected in >= 5 of 17 watersheds
GENE_MIN_BASINS_CRIT = 3  # step 1b: critical ARGs (rare by nature) are tested separately from >= 3 watersheds
GENE_K_SENS = (3, 5, 8)   # step 1b: sensitivity of the prevalence threshold
N_PERM_GENE = 2000        # step 1b: basin permutations (consumption profiles shuffled jointly across classes)
BETAS = (0.25, 0.5, 1.0)  # elasticity of resistome load to class consumption (scenario)
REDUCTION = 0.5           # relative cut in a group's use of a class (scenario)

# Relative tiers on the 0-100 scale (apply to both PORI and AMI)
TIER_EDGES = [20, 40, 60, 80]
SOURCE_CUT = 0.50         # source-based AMI: maximum relative reduction of the anthropogenic-taxa signal at a station
SOURCE_REF_Q = 0.25       # reference level of the signal = this quantile of the 36 stations (cleaner stations); no gain below it
SOURCE_BOOT = 2000        # cluster (basin) bootstrap draws for the elasticities
VALID_PERM = 1000         # permutations for the leave-one-watershed-out check of the source-based AMI
TIER_LABELS = ["Low", "Moderate", "High", "Very high", "Extreme"]

# Definition of "critical" ARGs: "curated" (49 genes listed with references in the
# workbook sheet 'High Risk Pathogens and ARGs') or "keyword" (broad name matching).
CRITICAL_ARG_MODE = os.environ.get("AMR_CRIT", "curated")   # "curated" or "keyword"

# --- Clinically relevant ("high-alert") pathogen species (substring match on taxon name).
# Rationale: ESKAPE pathogens (Rice 2008), WHO Bacterial Priority Pathogens List 2024,
# and waterborne / foodborne / zoonotic agents relevant to coastal waters.
# High-alert pathogens (list of the authors, sheet "High Risk Pathogens and ARGs"). Species are matched by name; the genera are matched
# by the first word of the taxon name, so that all their species (and unresolved genus-level reads) count, but not the family level
# (e.g. Enterobacteriaceae) nor other genera that merely contain the string (e.g. Pseudocitrobacter).
HIGH_ALERT_SPECIES = [
    "enterococcus faecium", "staphylococcus aureus", "klebsiella pneumoniae", "klebsiella aerogenes", "klebsiella quasipneumoniae",
    "acinetobacter baumannii", "pseudomonas aeruginosa",
    "stenotrophomonas maltophilia", "acinetobacter pittii", "acinetobacter junii",
    "acinetobacter lwoffii", "ochrobactrum anthropi", "providencia rettgeri",
    "morganella morganii", "serratia marcescens", "vibrio cholerae", "vibrio vulnificus",
    "vibrio parahaemolyticus", "aeromonas dhakensis", "aeromonas hydrophila",
    "aeromonas caviae", "aeromonas veronii", "aliarcobacter butzleri",
    "cronobacter sakazakii", "staphylococcus epidermidis", "staphylococcus haemolyticus",
    "escherichia coli",
]
HIGH_ALERT_GENERA = ["enterobacter", "citrobacter", "streptococcus"]


def is_high_alert(name):
    tl = str(name).lower()
    return any(sp in tl for sp in HIGH_ALERT_SPECIES) or (tl.split()[0] if tl.split() else "") in HIGH_ALERT_GENERA

# --- Curated critical ARGs (last-resort drugs, carbapenemases, ESBL/AmpC, VRE/MRSA,
# oxazolidinone, 16S-rRNA methyltransferases, plasmid-mediated quinolone, erm MLSB).
CURATED_CRITICAL_ARGS = [
    "mcr-1.10", "mcr-3.1", "mcr-9.1", "mcr-10.1", "tet(X)", "tet(X4)", "blaKPC-40",
    "blaNDM-20", "blaNDM-28", "blaNDM-43", "blaNDM-44", "blaVIM-1", "blaVIM-2",
    "blaVIM-79", "blaIMP-4", "blaIMP-15", "blaOXA-23", "blaOXA-48", "blaOXA-58",
    "blaOXA-72", "blaOXA-181", "blaOXA-232", "vanA", "vanB", "vanM", "vanX-vanA",
    "mecA", "mecC", "optrA", "cfrA", "cfrC", "cfrE", "armA", "rmtB", "rmtD2", "rmtF",
    "rmtH", "AAC(6')-Ib-cr3", "AAC(6')-Ib-cr7", "blaACT-5", "blaACT-40",
    "blaCTX-M-55", "blaCTX-M-199", "blaPER-1", "blaPER-16", "blaVEB-1", "blaSHV-12",
    "erm(B)", "erm(F)",
]
# Broad alternative: any gene whose name contains one of these fragments.
CRITICAL_ARG_KEYWORDS = [
    "qnr", "qep", "blaact", "blakpc", "blandm", "blavim", "blaimp", "blages", "blagim",
    "bladim", "blaoxa", "blabkc", "blafri", "blabic", "blacfi", "blal1", "blalmb",
    "blatmb", "blakhm", "vim-64", "vana", "vanb", "vanf", "vang", "vanm", "vann", "vanh",
    "vanx", "vany", "vanz", "meca", "mecc", "blaz-mecc", "mcr-", "tet(x", "optra", "cfra",
    "cfrc", "cfre", "arma", "rmtb", "rmtd2", "rmtf", "rmth", "-cr3", "-cr6", "-cr7", "erm(",
]

# --- Antibiotic classes: regex to find consumption columns in the 'Antibiotics' sheet,
# and the matching pre-aggregated resistome columns in the 'args' sheet.
AB_CLASSES = {
    "Tetracyclines":         (r"tetracycl",                                ["Tetracyclines"]),
    "Amphenicols":           (r"amphen",                                   ["Phenicols"]),
    "Penicillins":           (r"penicill",                                 ["Beta-lactams"]),
    "Cephalosporins I-II":   (r"cephalosporins.*\b(i-ii|i and ii)\b",      ["Cephalosporines 1-2"]),
    "Cephalosporins III-IV": (r"cephalosporins.*(iii-iv|iii and iv)",      ["Cephalosporines 3-4"]),
    "Sulfonamides+TMP":      (r"sul[fp]onamides",                          ["Sulfonamides", "Trimethoprim"]),
    "Macrolides":            (r"macrolid",                                 ["Macrolides"]),
    "Lincosamides":          (r"l[yi]ncosamid",                            ["Lincosamides"]),
    "Aminoglycosides":       (r"aminoglyc",                                ["Aminoglycosides"]),
    "Fluoroquinolones":      (r"fluoroquin",                               ["Fluoroquinolones"]),
    "Polymyxins":            (r"pol[iy]mixins",                            ["Polymyxins"]),
}
# Animal groups: key used to find the group's columns in the 'Antibiotics' sheet.
ANIMAL_KEYS = {"Cattle": "cattle", "Pigs": "pigs", "Poultry": "poultry",
               "Sheep&Goats": "sheep", "Finfish": "finfish", "Molluscs": "molluscs",
               "Horses": "horses", "Rabbits": "rabbits", "Dogs": "dogs"}
GROUPS = ["Humans"] + list(ANIMAL_KEYS)

# Palette (validated reference palette; sequential blue ramp for ordinal tiers)
PAL = dict(ink="#0b0b0b", ink2="#52514e", surface="#fcfcfb", grid="#e3e2de",
           blue="#2a78d6", orange="#eb6834", aqua="#1baf7a", yellow="#eda100",
           magenta="#e87ba4")
TIER_COL = {"Low": "#f6c5a1", "Moderate": "#ee9c63", "High": "#e0702a",
            "Very high": "#a8460f", "Extreme": "#3a1a0c"}
SEQ_RAMP = ["#fff3e8", "#f6c5a1", "#ee9c63", "#e0702a", "#a8460f", "#3a1a0c"]   # light peach -> dark brown
DIV_RAMP = ["#3f7f86", "#f2efe9", "#c24f12"]                                    # teal - neutral - orange-brown
ACCENT, TEAL = "#e0702a", "#3f7f86"
# Outlier-count bins used on the map (count -> colour); 0 flagged = white
OUT_BINS = [(1, 1, "#fbe3d0"), (2, 2, "#f6c5a1"), (3, 3, "#ee9c63"), (4, 4, "#e0702a"),
            (5, 5, "#c24f12"), (6, 7, "#8a3a10"), (8, 99, "#3a1a0c")]

LOG = []
QUIET = [False]   # set to True to silence logging during the sensitivity re-run


def log(*a):
    """Print and keep a copy for run_log.txt."""
    if QUIET[0]:
        return
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.append(s)


# ======================================================================
# HELPERS
# ======================================================================
def clean_cols(df):
    """Strip and collapse whitespace in column names."""
    df.columns = [re.sub(r"\s+", " ", str(c)).strip() for c in df.columns]
    return df


def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values (q-values); NaNs are ignored."""
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    if pv.size == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order] * pv.size / np.arange(1, pv.size + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(pv)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def minmax(v):
    """Rescale to 0-1 (0 if the vector is constant)."""
    v = np.asarray(v, float)
    r = np.nanmax(v) - np.nanmin(v)
    return (v - np.nanmin(v)) / r if r > 0 else np.zeros_like(v)


def tier(score):
    """Map a 0-100 score to its relative tier label."""
    idx = np.digitize(np.asarray(score, float), TIER_EDGES)
    return np.array(TIER_LABELS)[idx]


def spear(x, y):
    """Spearman rho and p-value; (nan, nan) when a vector is constant or n < 4."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4 or np.ptp(x[m]) == 0 or np.ptp(y[m]) == 0:
        return np.nan, np.nan
    r, p = spearmanr(x[m], y[m])
    return float(r), float(p)


def basin_mean(x, wsh):
    """Average a station-level vector within watersheds (sorted by watershed name)."""
    return pd.Series(np.asarray(x, float)).groupby(wsh).mean()


def find_col(df, prefix):
    """First column whose name starts with `prefix` (robust to unit symbols)."""
    hit = [c for c in df.columns if c.lower().startswith(prefix.lower())]
    if not hit:
        raise KeyError(f"No column starting with '{prefix}'")
    return hit[0]


# ======================================================================
# 0. DATA LOADING AND HARMONISATION
# ======================================================================
STATION_NAMES = {"Delta Po Sud": "Delta Po South", "Delta Po centro": "Delta Po Centre", "Delta Po Nord": "Delta Po North"}


def read_sheet(name, **kw):
    df = clean_cols(pd.read_excel(INPUT_XLSX, sheet_name=name, **kw))
    for c in df.columns:                                   # English station names in every table and figure
        if str(c).strip().lower() == "location":
            df[c] = df[c].replace(STATION_NAMES)
    return df


def load_data():
    """Read all sheets and align every station-level table on the Location key."""
    D = {}
    phys = read_sheet("Physical data").iloc[:N_STATIONS].copy()
    for c in ("Location", "Substrate", "Watershed"):
        phys[c] = phys[c].astype(str).str.strip()
    D["meta"] = phys[["Location", "Substrate", "Watershed"]].reset_index(drop=True)
    loc = D["meta"]["Location"].values

    def aligned(name):
        df = read_sheet(name).iloc[:N_STATIONS].copy()
        df["Location"] = df["Location"].astype(str).str.strip()
        return df.set_index("Location").loc[loc].reset_index()  # KeyError if a station is missing

    D["phys"] = aligned("Physical data")
    D["chem"] = aligned("Chemical microbiological data")
    D["animals"] = aligned("Animal Abundances")
    D["ab"] = aligned("Antibiotics")
    D["args"] = aligned("args")
    D["b16"] = aligned("Bacterial community (16S)")

    # Taxon annotation (Pathogenicity P/PP/N, Origin W/E/U, habitat flags Y/N)
    tax = clean_cols(pd.read_excel(INPUT_XLSX, sheet_name="bacteria data", index_col=0).T)
    tax.index = tax.index.astype(str).str.strip()
    D["tax"] = tax

    # Gene -> resistance-class text (row 3 of the 'args data' sheet)
    ad = pd.read_excel(INPUT_XLSX, sheet_name="args data", header=None, nrows=3)
    names = ad.iloc[0, 1:].astype(str).str.strip()
    classes = ad.iloc[2, 1:].fillna("").astype(str)
    D["gene_class_text"] = dict(zip(names, classes))
    return D


# ======================================================================
# 0b. FEATURE CONSTRUCTION
# ======================================================================
def build_outcomes(D, mode=None):
    """
    Five resistome / pathogen components used by the PORI:
      Pathogens_Total  : relative abundance (reads per million) of taxa annotated P or PP
      Alert_Pathogens  : relative abundance of clinically relevant (high-alert) species
      ARGs_Total       : summed normalised abundance of all ARGs
      Critical_ARGs    : summed normalised abundance of critical ARGs
      Mobility         : mobile / (mobile + chromosomal) ARG fraction
    16S counts are converted to reads per million so that sequencing depth does not
    enter the signal (microbiome data are compositional; Gloor et al. 2017).
    """
    b16, tax, args = D["b16"], D["tax"], D["args"]
    taxa = [c for c in b16.columns if c not in ("Location", "Substrate", "Watershed")]
    cnt = b16[taxa].apply(pd.to_numeric, errors="coerce").fillna(0)
    depth = cnt.sum(axis=1).values

    patho = tax["Pathogenicity"].astype(str).str.strip().str.upper()
    path_all = [t for t in taxa if patho.get(t) in ("P", "PP")]
    path_alert = [t for t in path_all if is_high_alert(t)]

    # ARG columns = all gene columns after the 'Mobile (mobilizable)' summary column
    cols = list(args.columns)
    first_gene = cols.index("Mobile (mobilizable)") + 1
    genes = cols[first_gene:]
    GM = args[genes].apply(pd.to_numeric, errors="coerce").fillna(0).values

    mode = mode or CRITICAL_ARG_MODE
    low = {g.lower(): i for i, g in enumerate(genes)}
    cur_idx = [low[g.lower()] for g in CURATED_CRITICAL_ARGS if g.lower() in low]
    kw_idx = [i for i, g in enumerate(genes) if any(k in g.lower() for k in CRITICAL_ARG_KEYWORDS)]
    crit_idx = cur_idx if mode == "curated" else kw_idx

    chrom = pd.to_numeric(args["Chromosomal"], errors="coerce").fillna(0).values
    mobile = pd.to_numeric(args["Mobile (mobilizable)"], errors="coerce").fillna(0).values

    rel = lambda cs: cnt[cs].sum(axis=1).values / depth * 1e6
    raw = pd.DataFrame({
        "Pathogens_Total": rel(path_all),
        "Alert_Pathogens": rel(path_alert),
        "ARGs_Total": GM.sum(axis=1),
        "Critical_ARGs": GM[:, crit_idx].sum(axis=1),
        "Mobility": mobile / (mobile + chrom + 1e-9),
    })
    C = dict(taxa=taxa, cnt=cnt, depth=depth, path_all=path_all, path_alert=path_alert,
             genes=genes, GM=GM, crit_idx=crit_idx, cur_idx=cur_idx, kw_idx=kw_idx, raw=raw, mode=mode)
    log(f"[data] pathogens={len(path_all)} (high-alert={len(path_alert)}), ARGs={len(genes)}, "
        f"critical ARGs ({mode})={len(crit_idx)} (curated={len(cur_idx)}, keyword={len(kw_idx)})")
    return C


def build_ab_structure(ab):
    """
    Consumption by antibiotic class x host group (humans + 9 animal groups) per station.
    Consumption is a basin-level quantity (total weight / annual water flow), so stations
    of the same basin share the same values. The sum over groups must equal the class
    total ('TOT ...' columns); this is checked.
    """
    cols = list(ab.columns)
    num = lambda cs: (ab[cs].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1).values
                      if cs else np.zeros(len(ab)))
    cls = {}
    for c, (rx, _) in AB_CLASSES.items():
        hs = [x for x in cols if x.startswith("AB HS ") and re.search(rx, x, re.I)]
        tot = [x for x in cols if x.startswith("TOT ") and re.search(rx, x, re.I)]
        d = {"Humans": num(hs)}
        for g, k in ANIMAL_KEYS.items():
            d[g] = num([x for x in cols
                        if not x.startswith(("AB HS", "AB VET", "TOT", "Overall"))
                        and re.search(rx, x, re.I) and k in x.lower()])
        d["TOT"] = num(tot)
        cls[c] = d
    worst = max(np.nanmax(np.abs(sum(v for g, v in d.items() if g != "TOT") - d["TOT"])
                          / (d["TOT"] + 1e-30)) for d in cls.values())
    log(f"[data] consumption check: max relative gap between sum of groups and class total = {worst:.2e}")
    return cls


# ======================================================================
# STEP 1 - EXPLORATORY CORRELATIONS (consumption vs resistome)
# ======================================================================
def step1_correlations(D, C, AB):
    """
    Spearman correlations between antibiotic consumption and resistome.
    Basin level (n = 17 watersheds) is the primary analysis because consumption is
    a basin attribute; station-level rho (n = 36) is shown for comparison only
    (stations of one basin are not independent).
      Family A: class-matched tests (consumption of class X vs ARG load of class X)
      Family B: consumption / animal abundance by host group vs the five components
    Each family is corrected with Benjamini-Hochberg FDR.
    """
    wsh = D["meta"]["Watershed"].values
    args = D["args"]

    rows = []
    for c, (_, argcols) in AB_CLASSES.items():
        arg = args[argcols].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1).values
        x = AB[c]["TOT"]
        rb, pb = spear(basin_mean(x, wsh), basin_mean(arg, wsh))
        rs, ps = spear(x, arg)
        rows.append(dict(Class=c, rho_basin=rb, p_basin=pb, rho_station=rs, p_station=ps))
    A = pd.DataFrame(rows)
    A["q_basin"] = bh_fdr(A["p_basin"])
    A["Supported"] = (A["q_basin"] < FDR_Q_GATE) & (A["rho_basin"] > 0)

    comp = {"Pathogens_Total": C["raw"]["Pathogens_Total"].values,
            "Alert_Pathogens": C["raw"]["Alert_Pathogens"].values,
            "ARGs_Total": C["raw"]["ARGs_Total"].values,
            "Critical_ARGs": C["raw"]["Critical_ARGs"].values,
            "Mobility": C["raw"]["Mobility"].values}
    ab = D["ab"]
    preds = {"Overall AB": ab["Overall AB"], "Overall AB (veterinary)": ab["Overall AB VET"],
             "Overall AB (human)": ab["Overall AB HS"]}
    for g in GROUPS:
        preds[f"Consumption - {g}"] = sum(AB[c][g] for c in AB_CLASSES)
    an = D["animals"]
    for c in [c for c in an.columns if "Abund" in c and "Rel Prop" not in c]:
        preds[f"Abundance - {c.split(' ')[0]}"] = an[c]
    rows = []
    for pn, pv in preds.items():
        pv = pd.to_numeric(pd.Series(np.asarray(pv)), errors="coerce").values
        for cn, cv in comp.items():
            rb, pb = spear(basin_mean(pv, wsh), basin_mean(cv, wsh))
            rs, ps = spear(pv, cv)
            rows.append(dict(Predictor=pn, Component=cn, rho_basin=rb, p_basin=pb,
                             rho_station=rs, p_station=ps))
    B = pd.DataFrame(rows)
    B["q_basin"] = bh_fdr(B["p_basin"])

    # Heatmap matrices (watershed level, n = 17): host-group population vs class consumption
    an = D["animals"]
    pref = {"Humans": "Humans", "Cattle": "Cattle", "Pigs": "Pigs", "Poultry": "Poultry", "Sheep&Goats": "Sheep",
            "Finfish": "Finfish", "Molluscs": "Molluscs", "Horses": "Horses", "Rabbits": "Rabbits", "Dogs": "Dogs"}
    pop = {g: pd.to_numeric(an[[c for c in an.columns if "Abund" in c and c.startswith(pf)][0]]).values
           for g, pf in pref.items()}
    cons = {"Total": {c: AB[c]["TOT"] for c in AB_CLASSES},
            "Human use": {c: AB[c]["Humans"] for c in AB_CLASSES},
            "Veterinary use": {c: sum(AB[c][g] for g in ANIMAL_KEYS) for c in AB_CLASSES}}
    rows = []
    for panel, rowgroups in (("Total", GROUPS), ("Human use", ["Humans"]), ("Veterinary use", list(ANIMAL_KEYS))):
        for g in rowgroups:
            for c in AB_CLASSES:
                r, pv = spear(basin_mean(pop[g], wsh), basin_mean(cons[panel][c], wsh))
                rows.append(dict(Panel=panel, Population=g, Class=c, rho_basin=r, p_basin=pv))
    H = pd.DataFrame(rows)
    H["q_basin_within_panel"] = H.groupby("Panel")["p_basin"].transform(lambda x: bh_fdr(x.values))

    # Full matrix: every total-consumption class (13 'TOT' columns) vs every ARG class column (18)
    tot_cols = [c for c in D["ab"].columns if c.startswith("TOT ")]
    arg_cls = list(args.columns)[list(args.columns).index("Watershed") + 1: list(args.columns).index("Chromosomal")]
    rows = []
    for tc in tot_cols:
        x = pd.to_numeric(D["ab"][tc], errors="coerce").values
        for ac in arg_cls:
            y = pd.to_numeric(args[ac], errors="coerce").values
            rb, pb = spear(basin_mean(x, wsh), basin_mean(y, wsh))
            rs, ps = spear(x, y)
            rows.append(dict(Consumption=tc.replace("TOT ", ""), ARG_class=ac, rho_basin=rb, p_basin=pb, rho_station=rs, p_station=ps))
    M = pd.DataFrame(rows)
    M["q_basin"] = bh_fdr(M["p_basin"])

    log(f"[step1] family A (class-matched, n=17 basins): {int(A['Supported'].sum())}/{len(A)} classes "
        f"supported (q<{FDR_Q_GATE}, rho>0); nominal p<0.05: {int((A.p_basin < 0.05).sum())}; "
        f"median rho = {A.rho_basin.median():.2f}")
    log(f"[step1] family B ({len(B)} tests): nominal p<0.05 = {int((B.p_basin < 0.05).sum())} "
        f"(expected by chance ~{0.05 * B.p_basin.notna().sum():.1f}); q<0.10 = {int((B.q_basin < 0.10).sum())}")
    return dict(A=A, B=B, H=H, M=M, gate=dict(zip(A["Class"], A["Supported"])))



# ======================================================================
# STEP 1b - GENE-LEVEL CONSUMPTION-RESISTOME TEST (captures cross- and multi-resistance)
# ======================================================================
def _rank_z(M):
    """Column-wise ranks (ties averaged) standardised to mean 0, unit SD; constant columns -> NaN."""
    R = np.apply_along_axis(rankdata, 0, M)
    sd = R.std(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        Z = (R - R.mean(axis=0)) / sd
    return Z


def _rho_p(Zg, Zx, n):
    """Spearman rho matrix (genes x classes) and two-sided p from the t approximation (df = n - 2)."""
    rho = Zg.T @ Zx / n
    r = np.clip(rho, -0.999999, 0.999999)
    t = r * np.sqrt((n - 2) / (1 - r ** 2))
    return rho, 2 * stats.t.sf(np.abs(t), n - 2)


def step1b_gene_level(D, C, AB):
    """
    Gene-level test of the consumption-resistome link at watershed level (n = 17).
    Class-level tests assume that a gene responds only to its own antibiotic class; in the
    environment multidrug genes, cross-resistance and co-selection on mobile elements break
    that assumption. Here every gene is tested against the consumption of EVERY class and the
    genes that pass are then labelled "matched" (the class tested is one the gene confers
    resistance to) or "cross" (it is not).
      * genes detected in >= GENE_MIN_BASINS watersheds (family 1) and the curated critical
        ARGs detected in >= GENE_MIN_BASINS_CRIT watersheds (family 2, exploratory);
      * Spearman rho per gene x class, BH-FDR within each family;
      * permutation of watersheds (consumption profiles shuffled jointly, so the correlation
        structure among genes and among classes is kept) for: (i) a global max-statistic test,
        (ii) an empirical false-discovery estimate, (iii) the class-specificity test, i.e. is
        the mean rho of matched pairs higher than that of non-matched pairs?
    Descriptive and exploratory: with 17 watersheds, consumption classes covary strongly.
    """
    rng = np.random.default_rng(SEED)
    wsh = D["meta"]["Watershed"].values
    genes, GM = C["genes"], C["GM"]
    classes = list(AB_CLASSES)
    Bm = pd.DataFrame(GM, columns=genes).groupby(wsh).mean()
    n = len(Bm)
    X = pd.DataFrame({c: basin_mean(AB[c]["TOT"], wsh) for c in classes}).loc[Bm.index]
    nb = (Bm > 0).sum(axis=0).values
    Zx = _rank_z(X.values)
    gcl = {g: classes_of(D["gene_class_text"].get(g, "")) for g in genes}
    crit_set = {genes[i] for i in C["cur_idx"]}
    isc = np.array([g in crit_set for g in genes])
    match = np.array([[c in gcl[g] for c in classes] for g in genes])
    assigned = match.any(axis=1)
    nclass = np.array([len(gcl[g]) for g in genes])

    def analyse(sel, label, n_perm=N_PERM_GENE):
        idx = np.where(sel)[0]
        Zg = _rank_z(Bm.values[:, idx])
        ok = np.isfinite(Zg).all(axis=0)
        idx, Zg = idx[ok], Zg[:, ok]
        rho, p = _rho_p(Zg, Zx, n)
        q = bh_fdr(p.ravel()).reshape(p.shape)
        m = match[idx]
        a = assigned[idx]
        def spec(rh, msk):                       # class-specificity statistic
            u = a[:, None] & ~msk
            v = a[:, None] & msk
            return rh[v].mean() - rh[u].mean() if v.any() and u.any() else np.nan
        stat = {"all": spec(rho, m)}
        multi = (nclass[idx] > 1)
        sing = (nclass[idx] == 1)
        def spec_sub(rh, msk, rows):
            if rows.sum() == 0:
                return np.nan
            sub_a = a & rows
            u = sub_a[:, None] & ~msk
            v = sub_a[:, None] & msk
            return rh[v].mean() - rh[u].mean() if v.any() and u.any() else np.nan
        stat["single-class genes"] = spec_sub(rho, m, sing)
        stat["multi-class genes"] = spec_sub(rho, m, multi)
        # permutation of watersheds (joint across classes)
        mx_obs = np.nanmax(rho)
        n_nom_obs = int(((p <= 0.01) & (rho > 0)).sum())
        pm_max, pm_nom = [], []
        pm_stat = {k: [] for k in stat}
        for _ in range(n_perm):
            perm = rng.permutation(n)
            rp, pp = _rho_p(Zg, Zx[perm], n)
            pm_max.append(np.nanmax(rp))
            pm_nom.append(int(((pp <= 0.01) & (rp > 0)).sum()))
            pm_stat["all"].append(spec(rp, m))
            pm_stat["single-class genes"].append(spec_sub(rp, m, sing))
            pm_stat["multi-class genes"].append(spec_sub(rp, m, multi))
        pm_max, pm_nom = np.array(pm_max), np.array(pm_nom)
        out = dict(idx=idx, rho=rho, p=p, q=q, match=m, assigned=a, label=label,
                   max_rho=mx_obs, p_global=(1 + (pm_max >= mx_obs).sum()) / (n_perm + 1),
                   n_nom_obs=n_nom_obs, n_nom_null=float(pm_nom.mean()),
                   fdr_emp=(min(1.0, pm_nom.mean() / n_nom_obs) if n_nom_obs else np.nan),
                   spec={k: (stat[k], (1 + (np.array(pm_stat[k]) >= stat[k]).sum()) / (n_perm + 1)
                             if np.isfinite(stat[k]) else np.nan) for k in stat},
                   null_stat=np.array(pm_stat["all"]))
        return out

    fam1 = analyse(nb >= GENE_MIN_BASINS, f"Family 1: genes in >= {GENE_MIN_BASINS} watersheds")
    fam2 = analyse(isc & (nb >= GENE_MIN_BASINS_CRIT),
                   f"Family 2 (exploratory): critical ARGs in >= {GENE_MIN_BASINS_CRIT} watersheds")

    def table(F):
        rows = []
        for gi, g in enumerate(F["idx"]):
            for ci, c in enumerate(classes):
                rows.append(dict(Gene=genes[g], Gene_classes="; ".join(sorted(gcl[genes[g]])) or "(unassigned)",
                                 Critical=bool(isc[g]), Basins_detected=int(nb[g]), Class_tested=c,
                                 Pair_type=("matched" if F["match"][gi, ci] else ("cross" if F["assigned"][gi] else "unassigned gene")),
                                 rho=F["rho"][gi, ci], p=F["p"][gi, ci], q_BH=F["q"][gi, ci]))
        return pd.DataFrame(rows).sort_values("p").reset_index(drop=True)
    T1, T2 = table(fam1), table(fam2)
    sig = lambda T: T[(T.q_BH < FDR_Q_GATE) & (T.rho > 0)]
    S1_, S2_ = sig(T1), sig(T2)

    # class-level summary and gate: a class is supported if >= 1 gene that confers resistance to it
    # (matched pair) passes q < 0.10 with rho > 0 in either family
    rows = []
    gate_n = {}
    for c in classes:
        mt1 = T1[(T1.Class_tested == c) & (T1.Pair_type == "matched")]
        mt2 = T2[(T2.Class_tested == c) & (T2.Pair_type == "matched")]
        npass = int(((mt1.q_BH < FDR_Q_GATE) & (mt1.rho > 0)).sum() + ((mt2.q_BH < FDR_Q_GATE) & (mt2.rho > 0)).sum())
        allp = T1[T1.Class_tested == c]
        rows.append(dict(Class=c, Genes_matched_tested=len(mt1), Matched_with_rho_pos_p_lt_0_05=int(((mt1.p < 0.05) & (mt1.rho > 0)).sum()),
                         Median_rho_matched=mt1.rho.median(), Median_rho_non_matched=allp[allp.Pair_type != "matched"].rho.median(),
                         Matched_passing_FDR_fam1=int(((mt1.q_BH < FDR_Q_GATE) & (mt1.rho > 0)).sum()),
                         Matched_passing_FDR_fam2=int(((mt2.q_BH < FDR_Q_GATE) & (mt2.rho > 0)).sum()),
                         Cross_genes_passing_FDR=int(((allp.Pair_type == "cross") & (allp.q_BH < FDR_Q_GATE) & (allp.rho > 0)).sum())))
        gate_n[c] = npass
    CL = pd.DataFrame(rows)

    # global / specificity summary
    srows = []
    for F in (fam1, fam2):
        srows.append(dict(Family=F["label"], Genes_tested=len(F["idx"]), Pairs_tested=F["rho"].size,
                          Pairs_q_lt_0_10_rho_pos=int(((F["q"] < FDR_Q_GATE) & (F["rho"] > 0)).sum()),
                          Max_rho=F["max_rho"], P_global_max_rho_perm=F["p_global"],
                          Pairs_p_le_0_01_rho_pos_observed=F["n_nom_obs"], Same_expected_under_permutation=F["n_nom_null"],
                          Empirical_FDR_at_p_0_01=F["fdr_emp"],
                          Specificity_all_genes_mean_rho_diff=F["spec"]["all"][0], P_specificity_all=F["spec"]["all"][1],
                          Specificity_single_class_diff=F["spec"]["single-class genes"][0], P_single=F["spec"]["single-class genes"][1],
                          Specificity_multi_class_diff=F["spec"]["multi-class genes"][0], P_multi=F["spec"]["multi-class genes"][1]))
    SUM = pd.DataFrame(srows)

    # sensitivity of the prevalence threshold (family 1 definition)
    krows = []
    for k in GENE_K_SENS:
        F = analyse(nb >= k, f"k={k}", n_perm=500)
        krows.append(dict(Min_basins=k, Genes_tested=len(F["idx"]), Pairs_tested=F["rho"].size,
                          Pairs_q_lt_0_10_rho_pos=int(((F["q"] < FDR_Q_GATE) & (F["rho"] > 0)).sum()),
                          P_global_max_rho_perm=F["p_global"], Specificity_diff=F["spec"]["all"][0], P_specificity=F["spec"]["all"][1]))
    SENS = pd.DataFrame(krows)

    log(f"[step1b] genes tested: family 1 = {len(fam1['idx'])} (>= {GENE_MIN_BASINS} watersheds; {fam1['rho'].size} pairs), "
        f"family 2 (critical) = {len(fam2['idx'])} (>= {GENE_MIN_BASINS_CRIT} watersheds; {fam2['rho'].size} pairs)")
    log(f"[step1b] pairs with q<{FDR_Q_GATE} & rho>0: family 1 = {len(S1_)}, family 2 = {len(S2_)}; "
        f"max rho (fam 1) = {fam1['max_rho']:.2f}, global permutation p = {fam1['p_global']:.3f}")
    log(f"[step1b] class specificity (matched minus non-matched mean rho): all genes = {fam1['spec']['all'][0]:+.3f} "
        f"(perm p = {fam1['spec']['all'][1]:.3f}); multi-class genes = {fam1['spec']['multi-class genes'][0]:+.3f} "
        f"(perm p = {fam1['spec']['multi-class genes'][1]:.3f})")
    log(f"[step1b] classes with >= 1 matched gene passing the gene-level FDR: {sum(v > 0 for v in gate_n.values())}/{len(gate_n)}")
    return dict(T1=T1, T2=T2, sig=pd.concat([S1_.assign(Family=1), S2_.assign(Family=2)]), CL=CL, SUM=SUM, SENS=SENS,
                gate_n=gate_n, fam1=fam1, fam2=fam2, classes=classes)


# ======================================================================
# STEP 2 - DISTRIBUTION AND OUTLIER SCREENING
# ======================================================================
def describe_and_flag(frame, level):
    """
    Per-variable descriptors: skewness, Shapiro-Wilk (raw and log1p), and outliers by
    Tukey fences (1.5 x IQR) and by the modified z-score (|0.6745 (x - median)/MAD| > 3.5;
    Iglewicz & Hoaglin 1993; Leys et al. 2013). Returns a summary table and a boolean
    matrix of Tukey flags (rows = units, columns = variables).
    """
    rows, flags = [], {}
    for c in frame.columns:
        s = frame[c].astype(float)
        x = s.dropna().values
        q1, q3 = np.percentile(x, [25, 75])
        iqr = q3 - q1
        tuk = (s < q1 - 1.5 * iqr) | (s > q3 + 1.5 * iqr)
        med = np.median(x)
        mad = np.median(np.abs(x - med))
        mz = (0.6745 * (s - med) / mad).abs() > 3.5 if mad > 0 else pd.Series(False, index=s.index)
        sw = stats.shapiro(x)[1] if len(x) >= 3 and np.ptp(x) > 0 else np.nan
        swl = (stats.shapiro(np.log1p(x))[1]
               if (x >= 0).all() and np.ptp(x) > 0 and len(x) >= 3 else np.nan)
        rows.append(dict(Level=level, Variable=c, n=len(x), Median=med,
                         Skewness=stats.skew(x), Shapiro_p=sw, Shapiro_p_log1p=swl,
                         Zero_share=float((x == 0).mean()),
                         Tukey_outliers=int(tuk.sum()), ModZ_outliers=int(mz.sum())))
        flags[c] = tuk.values
    return pd.DataFrame(rows), pd.DataFrame(flags, index=frame.index)


def step2_outliers(D, C, AB):
    """Screen basin-level and station-level variables; count outliers per station."""
    meta, wsh = D["meta"], D["meta"]["Watershed"].values
    an = D["animals"]
    bvars = pd.DataFrame({"Inhabitants (x1000)": pd.to_numeric(D["phys"][find_col(D["phys"], "Inhabitants")])})
    for c in [c for c in an.columns if "Abund" in c and "Rel Prop" not in c]:
        bvars[f"Abundance {c.split(' ')[0]}"] = pd.to_numeric(an[c])
    for k, col in (("Overall AB", "Overall AB"), ("AB veterinary", "Overall AB VET"), ("AB human", "Overall AB HS")):
        bvars[k] = pd.to_numeric(D["ab"][col])
    for c in AB_CLASSES:
        bvars[f"Consumption {c}"] = AB[c]["TOT"]
    bvars = bvars.groupby(wsh).mean()                      # one row per basin (n = 17)

    chem = D["chem"]
    svars = pd.DataFrame({
        "Salinity": pd.to_numeric(D["phys"]["Salinity"]),
        "NO2": pd.to_numeric(chem[find_col(chem, "NO2")]), "NO3": pd.to_numeric(chem[find_col(chem, "NO3")]),
        "NH4": pd.to_numeric(chem[find_col(chem, "NH4")]), "PO4": pd.to_numeric(chem[find_col(chem, "PO4")]),
        "E. coli (ARPA)": pd.to_numeric(chem[find_col(chem, "ecoli mean")]),
        "Enterococci (ARPA)": pd.to_numeric(chem[find_col(chem, "enterococci mean")]),
    })
    for k, v in C["raw"].items():
        svars[k] = v.values

    sb, fb = describe_and_flag(bvars, "Basin (n=17)")
    ss, fs = describe_and_flag(svars, "Station (n=36)")
    S = pd.concat([sb, ss], ignore_index=True)

    fb_st = fb.loc[wsh].reset_index(drop=True)             # basin flags inherited by stations
    out = meta.copy()
    out["Flags_station_vars"] = fs.sum(axis=1).values
    out["Flags_basin_vars"] = fb_st.sum(axis=1).values
    out["Flag_share_%"] = (100 * (out.Flags_station_vars + out.Flags_basin_vars)
                           / (fs.shape[1] + fb.shape[1])).round(1)
    out["Station_variables_flagged"] = [", ".join(fs.columns[fs.iloc[i].values]) for i in range(len(out))]
    out["Basin_variables_flagged"] = [", ".join(fb.columns[fb_st.iloc[i].values]) for i in range(len(out))]

    nonnorm = (S.Shapiro_p < 0.05).mean() * 100
    log(f"[step2] {len(S)} variables screened: {nonnorm:.0f}% depart from normality (Shapiro p<0.05); "
        f"median |skewness| = {S.Skewness.abs().median():.2f}; "
        f"{(S.Tukey_outliers > 0).mean() * 100:.0f}% have >=1 Tukey outlier; "
        f"after log1p, {(S.Shapiro_p_log1p < 0.05).mean() * 100:.0f}% still non-normal")
    return dict(summary=S, stations=out)


def step2_zflags(D, C, AB, z_thr=3.5):
    """
    Station-level 'unusually high' flags used for the outlier map (Fig. 2).
    For each variable the robust z-score is computed on log(1 + x):
        z = (x - median) / (1.4826 * MAD)
    and a value is flagged when z > 3.5 (the modified-z cut-off of Iglewicz & Hoaglin 1993,
    high side only: we look for unusually HIGH loads or pressures; p = two-sided normal tail). Station-level variables use the 36 stations;
    basin-level variables (animal abundances) use the 17 watersheds and the flag is
    inherited by every station of the watershed. Variables with MAD = 0 are skipped.
    Flags are descriptive, not corrected for multiplicity.
    """
    args, chem, an, meta = D["args"], D["chem"], D["animals"], D["meta"]
    wsh = meta["Watershed"].values
    st_vars, ba_vars = {}, {}
    cols = list(args.columns)
    for c in cols[cols.index("Watershed") + 1: cols.index("Chromosomal")]:
        st_vars[c] = pd.to_numeric(args[c], errors="coerce").values
    st_vars["Mobile ARGs"] = C["raw"]["Mobility"].values
    st_vars["Salinity"] = pd.to_numeric(D["phys"]["Salinity"], errors="coerce").values
    for k in ("NO2", "NO3", "NH4", "PO4"):
        st_vars[k] = pd.to_numeric(chem[find_col(chem, k)], errors="coerce").values
    st_vars["E. coli (ARPA)"] = pd.to_numeric(chem[find_col(chem, "ecoli mean")], errors="coerce").values
    st_vars["Enterococci (ARPA)"] = pd.to_numeric(chem[find_col(chem, "enterococci mean")], errors="coerce").values
    for c in [c for c in an.columns if "Abund" in c and "Rel Prop" not in c]:
        ba_vars[c.split(" ")[0] + " (nr)"] = pd.to_numeric(an[c], errors="coerce").values

    def robust_z(x):
        x = np.log1p(np.maximum(np.asarray(x, float), 0))
        med = np.nanmedian(x)
        mad = np.nanmedian(np.abs(x - med)) * 1.4826
        return (x - med) / mad if mad > 0 else np.full_like(x, np.nan)

    rows = []
    for name, v in st_vars.items():
        z = robust_z(v)
        for i in np.where(z > z_thr)[0]:
            rows.append((meta["Location"][i], name, z[i], 2 * stats.norm.sf(z[i])))
    for name, v in ba_vars.items():
        bm = basin_mean(v, wsh)
        zb = robust_z(bm.values)
        zmap = dict(zip(bm.index, zb))
        for i, w in enumerate(wsh):
            if zmap[w] > z_thr:
                rows.append((meta["Location"][i], name, zmap[w], 2 * stats.norm.sf(zmap[w])))
    Z = pd.DataFrame(rows, columns=["Location", "Variable", "Z", "p"]).sort_values(["Location", "Z"], ascending=[True, False])
    cnt = Z.groupby("Location").size().reindex(meta["Location"]).fillna(0).astype(int)
    log(f"[step2] z-score flags (z>{z_thr}): {len(Z)} flags over {int((cnt > 0).sum())} stations; "
        f"stations with >=3 flags: {cnt[cnt >= 3].to_dict()}")
    return dict(flags=Z, counts=cnt)


# ======================================================================
# STEP 3 - PORI (descriptive, relative to the 36 stations)
# ======================================================================
def step3_pori(D, C, S1):
    """
    PORI combines five components rescaled 0-1 over the stations:
      pathogens (total, high-alert), ARGs (total, critical) and ARG mobility.
    PORI_raw = sqrt((a+e)(b+e)) * sqrt((p+e)(q+e)) * (1 + m), where a, b = ARG pillars,
    p, q = pathogen pillars and m = mobile fraction (all 0-1, e = EPS): the geometric mean
    of the two ARG pillars times that of the two pathogen pillars, amplified by mobility. A geometric structure means that a high score
    needs BOTH resistome and pathogens to be elevated. PORI is then rescaled to 0-100
    and classed into tiers. Sensitivity: random weights / offset / aggregation form.
    """
    raw = C["raw"]
    wsh = D["meta"]["Watershed"].values
    loc = D["meta"]["Location"].values
    lt, lc = np.log1p(raw["ARGs_Total"].values), np.log1p(raw["Critical_ARGs"].values)
    lp, lq = np.log1p(raw["Pathogens_Total"].values), np.log1p(raw["Alert_Pathogens"].values)
    N = {"Pathogens_Total": minmax(lp), "Alert_Pathogens": minmax(lq),
         "ARGs_Total": minmax(lt), "Critical_ARGs": minmax(lc), "Mobility": minmax(raw["Mobility"].values)}
    sq = lambda a, b: np.sqrt((a + EPS) * (b + EPS))

    def pori_raw(a, b, i=slice(None)):
        return sq(a, b) * sq(N["Pathogens_Total"][i], N["Alert_Pathogens"][i]) * (1 + N["Mobility"][i])

    r0 = pori_raw(N["ARGs_Total"], N["Critical_ARGs"])
    score = minmax(r0) * 100
    rank = rankdata(-score, method="min").astype(int)

    # --- Monte-Carlo sensitivity: Dirichlet weights, epsilon, geometric vs arithmetic ---
    rng = np.random.default_rng(SEED)
    keys = list(N)
    base_tier = tier(score)
    base_idx = np.digitize(score, TIER_EDGES)
    rc, ranks, top5 = [], [], np.zeros(len(score))
    same_tier, within_one = np.zeros(len(score)), np.zeros(len(score))
    for m in range(N_MC):
        w, e = rng.dirichlet(np.ones(5)), rng.uniform(0.01, 0.3)
        if m % 2 == 0:
            s = np.exp(sum(w[j] * np.log(N[k] + e) for j, k in enumerate(keys)))
        else:
            s = sum(w[j] * N[k] for j, k in enumerate(keys))
        sc = minmax(s) * 100
        rc.append(spearmanr(s, score)[0])
        ranks.append(rankdata(-s, method="min"))
        top5[np.argsort(-s)[:5]] += 1
        sim_idx = np.digitize(sc, TIER_EDGES)
        same_tier += (sim_idx == base_idx)
        within_one += (np.abs(sim_idx - base_idx) <= 1)
    ranks = np.array(ranks)

    # --- Alternative critical-ARG definition (keyword vs curated) ---
    alt_idx = C["kw_idx"] if C["mode"] == "curated" else C["cur_idx"]
    alt = np.log1p(C["GM"][:, alt_idx].sum(axis=1))
    score_alt = minmax(pori_raw(N["ARGs_Total"], minmax(alt))) * 100

    # --- External validation with ARPA indicators (not used in the index) ---
    chem = D["chem"]
    ecoli = pd.to_numeric(chem[find_col(chem, "ecoli mean")]).values
    ent = pd.to_numeric(chem[find_col(chem, "enterococci mean")]).values
    val = []
    for nm, v in (("E. coli (ARPA 2024 mean)", ecoli), ("Enterococci (ARPA 2024 mean)", ent)):
        rs, ps = spear(v, score)
        rb, pb = spear(basin_mean(v, wsh), basin_mean(score, wsh))
        val.append(dict(Check=f"PORI vs {nm}", rho_station=rs, p_station=ps, rho_basin=rb, p_basin=pb))

    # --- Structure of the index ---
    sub = D["meta"]["Substrate"].values
    groups = [score[sub == s] for s in np.unique(sub)]
    H, pkw = stats.kruskal(*groups)
    sst = ((score - score.mean()) ** 2).sum()
    ssw = ((score - pd.Series(score).groupby(wsh).transform("mean").values) ** 2).sum()

    out = D["meta"].copy()
    for k in keys:
        out[k + " (0-100)"] = (N[k] * 100).round(1)
    out["PORI_raw"] = r0
    out["PORI (0-100)"] = score.round(1)
    out["PORI tier"] = base_tier
    out["Rank"] = rank
    out["Rank_p5"] = np.percentile(ranks, 5, axis=0).astype(int)
    out["Rank_p95"] = np.percentile(ranks, 95, axis=0).astype(int)
    out["Top5_frequency"] = (top5 / N_MC).round(3)
    out["Tier_agreement_exact"] = (same_tier / N_MC).round(3)
    out["Tier_agreement_within_1"] = (within_one / N_MC).round(3)
    out = out.sort_values("Rank").reset_index(drop=True)

    comp_corr = pd.DataFrame(N).corr(method="spearman").round(2)
    kv = pd.DataFrame([
        ("Sensitivity: mean Spearman vs base ranking", np.mean(rc)),
        ("Sensitivity: 5th percentile Spearman", np.percentile(rc, 5)),
        ("Sensitivity: median exact tier agreement across draws", np.median(same_tier / N_MC)),
        ("Sensitivity: median agreement within one tier", np.median(within_one / N_MC)),
        (f"PORI ({C['mode']}) vs alternative critical-ARG definition (rho)", spear(score, score_alt)[0]),
        ("Kruskal-Wallis H across substrates", H), ("Kruskal-Wallis p", pkw),
        ("Share of PORI variance between basins", 1 - ssw / sst),
    ], columns=["Metric", "Value"])
    log(f"[step3] PORI tiers: {pd.Series(base_tier).value_counts().to_dict()}")
    log(f"[step3] sensitivity: mean rho={np.mean(rc):.2f} (p5={np.percentile(rc, 5):.2f}); tier agreement exact/within-1 = "
        f"{np.median(same_tier / N_MC):.2f}/{np.median(within_one / N_MC):.2f}; "
        f"alt. critical-ARG definition rho={spear(score, score_alt)[0]:.2f}")
    for v in val:
        log(f"[step3] validation {v['Check']}: station rho={v['rho_station']:.2f} (p={v['p_station']:.2f}); "
            f"basin rho={v['rho_basin']:.2f} (p={v['p_basin']:.2f})")
    log(f"[step3] substrates KW p={pkw:.2f}; between-basin share of variance = {1 - ssw / sst:.2f}")
    return dict(table=out, score=score, N=N, validation=pd.DataFrame(val), kv=kv, corr=comp_corr,
                pori_raw=pori_raw, r0=r0, lt=lt, lc=lc)


def step3c_substrate(D, P):
    """
    PORI by substrate (river, estuary, beach, touristic port, industrial harbour):
    n, mean, median, interquartile range, range, tier of the mean, and the share of
    stations in each tier. Substrate groups are small (4-10 stations), so the summary is
    descriptive; the Kruskal-Wallis test across substrates is reported in S3_sensitivity_structure.
    """
    df = D["meta"][["Location", "Substrate"]].copy()
    df["PORI"] = P["score"]
    df["tier"] = tier(P["score"])
    g = df.groupby("Substrate")["PORI"]
    T = pd.DataFrame({"n_stations": g.size(), "Mean_PORI": g.mean(), "Median_PORI": g.median(),
                      "Q1": g.quantile(0.25), "Q3": g.quantile(0.75), "Min": g.min(), "Max": g.max(),
                      "Highest_station": df.loc[g.idxmax(), "Location"].values})
    T["Tier_of_mean"] = tier(T["Mean_PORI"].values)
    for lab in TIER_LABELS:
        T[f"Stations_{lab}"] = df[df.tier == lab].groupby("Substrate").size().reindex(T.index).fillna(0).astype(int)
    T = T.sort_values("Mean_PORI", ascending=False).reset_index().round(1)
    log("[step3] PORI by substrate (mean): " + ", ".join(f"{r.Substrate} {r.Mean_PORI:.1f} (n={r.n_stations})" for r in T.itertuples()))
    return dict(table=T, points=df)


def step3b_spatial(D, C, P, n_perm=999):
    """
    Spatial structure (does proximity predict similarity?).
    (a) Mantel test (Spearman) between great-circle distance among stations and the
        Bray-Curtis dissimilarity of (i) all OTUs, (ii) pathogenic OTUs, (iii) ARG profiles,
        and the absolute difference in PORI; permutation p-values (999 station relabelings).
    (b) Moran's I of the PORI with inverse-distance weights (permutation p-value).
    A weak distance-decay and no autocorrelation for the PORI support assessing
    risk station by station rather than inferring it from neighbours.
    """
    from scipy.spatial.distance import pdist, squareform
    lat = np.radians(pd.to_numeric(D["phys"]["Lat (N)"]).values.astype(float))
    lon = np.radians(pd.to_numeric(D["phys"]["Long (E)"]).values.astype(float))
    n = len(lat)
    a = np.sin((lat[:, None] - lat[None]) / 2) ** 2 + np.cos(lat)[:, None] * np.cos(lat)[None] * np.sin((lon[:, None] - lon[None]) / 2) ** 2
    G = 2 * 6371 * np.arcsin(np.sqrt(a))
    iu = np.triu_indices(n, 1)
    rng = np.random.default_rng(SEED)
    rel = C["cnt"].div(C["depth"], axis=0)
    mats = {"All OTUs (Bray-Curtis)": squareform(pdist(rel.values, "braycurtis")),
            "Pathogenic OTUs (Bray-Curtis)": squareform(pdist(rel[C["path_all"]].values, "braycurtis")),
            "ARG profiles (Bray-Curtis)": squareform(pdist(C["GM"], "braycurtis")),
            "|delta PORI|": np.abs(P["score"][:, None] - P["score"][None])}
    rows = []
    for k, Dm in mats.items():
        r = spearmanr(G[iu], Dm[iu])[0]
        ge = sum(abs(spearmanr(G[iu], Dm[np.ix_(pm, pm)][iu])[0]) >= abs(r) for pm in (rng.permutation(n) for _ in range(n_perm)))
        rows.append(dict(Test="Mantel (Spearman) vs geographic distance", Variable=k, Statistic=r, p_perm=(ge + 1) / (n_perm + 1)))
    W = 1 / np.where(G == 0, np.nan, G); W = np.nan_to_num(W); W = W / W.sum(1, keepdims=True)
    z = P["score"] - P["score"].mean()
    mi = lambda v: (n / W.sum()) * (v @ W @ v) / (v @ v)
    I = mi(z)
    ge = sum(mi(rng.permutation(z)) >= I for _ in range(n_perm))
    rows.append(dict(Test="Moran's I (inverse-distance weights)", Variable="PORI", Statistic=I, p_perm=(ge + 1) / (n_perm + 1)))
    SP = pd.DataFrame(rows)
    log("[step3] spatial structure: " + "; ".join(f"{r.Variable} {r.Statistic:.2f} (p={r.p_perm:.3f})" for r in SP.itertuples()))
    return SP


# ======================================================================
# STEP 4 - SENTINEL TAXA AND ORIGIN-BASED SIGNALS
# ======================================================================
def step4_sentinels(D, C, P):
    """
    (a) Origin fractions: share of reads from taxa of anthropogenic (W), environmental (E)
        and ubiquitous/unknown (U) origin; W is also computed excluding pathogens to
        avoid circularity with the PORI.
    (b) Source signals: share of reads from non-pathogenic taxa whose habitat annotation
        is human-only (HUMAN=Y, ANIMAL=N), animal-only, gut or waste.
    (c) Sentinel search among environmental (E) species-level taxa present in >= 50% of
        stations: Spearman of log1p(reads per million) with the PORI, BH-FDR, then
        robustness checks (basin-level rho, leave-one-basin-out minimum rho, correlation
        with sequencing depth, PORI recomputed without the critical-ARG pillar).
    Sensitivity: the screen is repeated with U taxa included.
    """
    tax, cnt, depth, taxa = D["tax"], C["cnt"], C["depth"], C["taxa"]
    wsh, loc = D["meta"]["Watershed"].values, D["meta"]["Location"].values
    score = P["score"]
    org = tax["Origin"].astype(str).str.strip().str.upper()
    flag = lambda col: tax[col].astype(str).str.strip().str.upper()
    path = set(C["path_all"])
    frac = lambda cols: cnt[cols].sum(axis=1).values / depth * 100

    O = D["meta"][["Location", "Watershed"]].copy()
    for o in "WEU":
        O[f"Frac_{o}_%"] = frac([t for t in taxa if org.get(t) == o])
    O["Frac_W_nonpathogen_%"] = frac([t for t in taxa if org.get(t) == "W" and t not in path])
    H, A_, G_, W_ = flag("HUMAN"), flag("ANIMAL"), flag("GUT"), flag("WASTE")
    nonp = [t for t in taxa if t not in path]
    O["Human_only_taxa_%"] = frac([t for t in nonp if H.get(t) == "Y" and A_.get(t) == "N"])
    O["Animal_only_taxa_%"] = frac([t for t in nonp if A_.get(t) == "Y" and H.get(t) == "N"])
    O["Gut_taxa_%"] = frac([t for t in nonp if G_.get(t) == "Y"])
    O["Waste_taxa_%"] = frac([t for t in nonp if W_.get(t) == "Y"])
    O["Human_share_of_host_associated_%"] = 100 * O["Human_only_taxa_%"] / (
        O["Human_only_taxa_%"] + O["Animal_only_taxa_%"] + 1e-9)
    O["PORI (0-100)"] = score.round(1)

    chem = D["chem"]
    ecoli = pd.to_numeric(chem[find_col(chem, "ecoli mean")]).values
    ent = pd.to_numeric(chem[find_col(chem, "enterococci mean")]).values
    rows = []
    for c in [c for c in O.columns if c.endswith("_%") and c != "Human_share_of_host_associated_%"] + ["Human_share_of_host_associated_%"]:
        rb, _ = spear(basin_mean(O[c], wsh), basin_mean(score, wsh))
        rows.append(dict(Variable=c, rho_PORI_station=spear(O[c], score)[0], rho_PORI_basin=rb,
                         rho_ARPA_Ecoli=spear(O[c], ecoli)[0], rho_ARPA_Enterococci=spear(O[c], ent)[0]))
    OS = pd.DataFrame(rows)

    # ---- sentinel screen ----
    relm = cnt.div(C["depth"], axis=0) * 1e6
    lg = np.log1p(relm)
    prev = (cnt > 0).mean()
    # PORI without the critical-ARG pillar (the ARG total is used twice) as an independence check
    score_noc = minmax(P["pori_raw"](P["N"]["ARGs_Total"], P["N"]["ARGs_Total"])) * 100

    def screen(origins, exclude=(), min_prev=MIN_PREVALENCE):
        cand = [t for t in taxa if org.get(t) in origins and t not in exclude and " " in t and prev[t] >= min_prev]
        rows = []
        for t in cand:
            r, p = spear(lg[t], score)
            lobo = [spear(lg[t][wsh != w], score[wsh != w])[0] for w in np.unique(wsh)]
            rb, _ = spear(basin_mean(lg[t], wsh), basin_mean(score, wsh))
            rows.append(dict(Species=t, Origin=org.get(t), Prevalence_pct=100 * prev[t], rho=r, p=p,
                             rho_basin=rb, rho_LOBO_min=np.nanmin(lobo),
                             rho_depth=spear(lg[t], C["depth"])[0], rho_PORI_noCriticalARG=spear(lg[t], score_noc)[0],
                             Mean_reads_per_M=relm[t].mean()))
        d = pd.DataFrame(rows)
        d["q_FDR"] = bh_fdr(d["p"])
        d["Robust_sentinel"] = (d.q_FDR < FDR_Q_SENTINEL) & (d.rho > 0) & (d.rho_LOBO_min > SENT_MIN_LOBO)
        return d.sort_values("p").reset_index(drop=True)

    SE, SEU = screen({"E"}), screen({"E", "U"})
    SW = screen({"W"}, exclude=path, min_prev=0.25)      # exploratory: anthropogenic non-pathogenic taxa (pathogens excluded to avoid circularity); prevalence >= 25% because W taxa are sparse (only 1 reaches 50%)
    log(f"[step4] origin: W/E/U mean share = {O['Frac_W_%'].mean():.1f}/{O['Frac_E_%'].mean():.1f}/{O['Frac_U_%'].mean():.1f}%; "
        f"rho(W non-pathogen, PORI)={spear(O['Frac_W_nonpathogen_%'], score)[0]:.2f}; "
        f"rho(human-only taxa, PORI)={spear(O['Human_only_taxa_%'], score)[0]:.2f}; "
        f"rho(human-only taxa, ARPA E. coli)={spear(O['Human_only_taxa_%'], ecoli)[0]:.2f}")
    log(f"[step4] sentinels: {len(SE)} E-taxa tested; q<{FDR_Q_SENTINEL} & rho>0: {int(((SE.q_FDR < FDR_Q_SENTINEL) & (SE.rho > 0)).sum())}; "
        f"robust: {SE.loc[SE.Robust_sentinel, 'Species'].tolist()}; with U included: {len(SEU)} tested, "
        f"{int(((SEU.q_FDR < FDR_Q_SENTINEL) & (SEU.rho > 0)).sum())} significant")
    SW["rho_ARPA_Ecoli"] = [spear(lg[t], ecoli)[0] for t in SW.Species]
    SW["rho_ARPA_Enterococci"] = [spear(lg[t], ent)[0] for t in SW.Species]
    log(f"[step4] exploratory W non-pathogen screen: {len(SW)} taxa tested; q<{FDR_Q_SENTINEL} & rho>0: {int(((SW.q_FDR < FDR_Q_SENTINEL) & (SW.rho > 0)).sum())}; "
        f"robust: {SW.loc[SW.Robust_sentinel, 'Species'].tolist()}")
    return dict(origin=O, origin_summary=OS, sent=SE, sent_EU=SEU, sent_W=SW, lg=lg, ecoli=ecoli, ent=ent)


# ======================================================================
# STEP 5 - RELATIVE IMPACT OF HOST GROUPS PER STATION
# ======================================================================
# ======================================================================
# STEP 6 - AMI (what to act on) + source context
# ======================================================================
def classes_of(text):
    """Antibiotic classes affected by a gene, from the 'Resistance' text of 'args data'."""
    t = str(text).lower()
    K = set()
    if "tetracycl" in t: K.add("Tetracyclines")
    if "phenicol" in t: K.add("Amphenicols")
    if re.search(r"penicillin|carbenicillin|cloxacillin|oxacillin", t): K.add("Penicillins")
    if "cephalosporin" in t:
        a, b = bool(re.search(r"1st|2nd", t)), bool(re.search(r"3rd|4th", t))
        if a: K.add("Cephalosporins I-II")
        if b: K.add("Cephalosporins III-IV")
        if not a and not b: K |= {"Cephalosporins I-II", "Cephalosporins III-IV"}
    if re.search(r"sulfonamide|trimethoprim", t): K.add("Sulfonamides+TMP")
    if re.search(r"macrolide|mls", t): K.add("Macrolides")
    if re.search(r"lincosamide|mls", t): K.add("Lincosamides")
    if "aminoglycoside" in t: K.add("Aminoglycosides")
    if "quinolone" in t: K.add("Fluoroquinolones")
    if re.search(r"polymyxin|colistin", t): K.add("Polymyxins")
    return K


def step6_ami(D, C, AB, S1, P, S4):
    """
    AMI = relative priority of a mitigation target at each station.
    Scenario engine: cut the use of antibiotic class c in host group g by REDUCTION.
    The class's consumption falls by REDUCTION x share(g, c); the ARG load of class c
    is assumed to fall by beta x that amount (beta = elasticity of resistome to
    consumption, not identifiable from these data: beta=0 means no effect, beta=1 is
    the upper bound). Only genes acting on a single antibiotic class are used here; multidrug
    (multi-class) genes are treated in a supplementary analysis (step6c) that points to a host group, not to a class.
    The PORI is recomputed (pillars floored at the observed minimum) and
        dPORI = PORI_base - PORI_scenario     (in PORI points)
    AMI(station) = largest dPORI over all (group, class) targets, rescaled 0-100 over
    the stations and classed into the same relative tiers as the PORI. The ranking of
    targets does not depend on beta; beta only scales the size of the reduction.
    Each target also carries an evidence label from Step 1 (is class-matched
    consumption associated with class ARG load?), and each station an action
    suggestion from the source quadrant (PORI tier x human-source signal).
    """
    n = len(D["meta"])
    loc, wsh = D["meta"]["Location"].values, D["meta"]["Watershed"].values
    genes, GM = C["genes"], C["GM"]
    N, score = P["N"], P["score"]
    classes = list(AB_CLASSES)
    gcl = {g: classes_of(D["gene_class_text"].get(g, "")) for g in genes}
    crit_mask = np.zeros(len(genes), bool)
    crit_mask[C["crit_idx"]] = True
    A_tot, A_crit = np.zeros((n, len(classes))), np.zeros((n, len(classes)))
    for ci, c in enumerate(classes):
        w = np.array([1.0 if (len(gcl[g]) == 1 and c in gcl[g]) else 0.0 for g in genes])   # single-class genes only
        A_tot[:, ci] = GM @ w
        A_crit[:, ci] = GM[:, crit_mask] @ w[crit_mask]
    tot0, crit0 = GM.sum(axis=1), GM[:, crit_mask].sum(axis=1)
    log(f"[step6] genes mapped to >=1 antibiotic class: {sum(1 for g in genes if gcl[g])}/{len(genes)}; "
        f"single-class genes used for the AMI: {sum(1 for g in genes if len(gcl[g]) == 1)}; multi-class (supplementary): {sum(1 for g in genes if len(gcl[g]) > 1)}")

    lt, lc = np.log1p(tot0), np.log1p(crit0)
    pillar = lambda x, lo, hi: np.clip((np.log1p(np.maximum(x, 0)) - lo) / (hi - lo), 0, None)
    sq = lambda a, b: np.sqrt((a + EPS) * (b + EPS))
    base_raw = sq(pillar(tot0, lt.min(), lt.max()), pillar(crit0, lc.min(), lc.max())) * \
        sq(N["Pathogens_Total"], N["Alert_Pathogens"]) * (1 + N["Mobility"])
    lo, hi = base_raw.min(), base_raw.max()
    to_score = lambda r: (r - lo) / (hi - lo) * 100
    base = to_score(base_raw)
    assert np.allclose(base, score, atol=1e-6), "baseline must reproduce the observed PORI"
    floor_raw = sq(pillar(np.full(n, np.expm1(lt.min())), lt.min(), lt.max()),
                   pillar(np.full(n, np.expm1(lc.min())), lc.min(), lc.max())) * \
        sq(N["Pathogens_Total"], N["Alert_Pathogens"]) * (1 + N["Mobility"])
    ceiling = base - to_score(floor_raw)   # max PORI reduction achievable by ANY action on ARGs

    rows = []
    for i in range(n):
        for ci, c in enumerate(classes):
            tc = AB[c]["TOT"][i]
            if tc <= 0:
                continue
            for g in GROUPS:
                sh = AB[c][g][i] / tc
                if sh > 0:
                    rows.append((i, g, c, ci, sh))
    R = pd.DataFrame(rows, columns=["i", "Group", "Class", "ci", "share"])
    ii, ci_ = R.i.values, R.ci.values

    def dpori(beta, cut):
        f = beta * cut * R.share.values
        t = tot0[ii] - f * A_tot[ii, ci_]
        cr = crit0[ii] - f * A_crit[ii, ci_]
        r = sq(pillar(t, lt.min(), lt.max()), pillar(cr, lc.min(), lc.max())) * \
            sq(N["Pathogens_Total"][ii], N["Alert_Pathogens"][ii]) * (1 + N["Mobility"][ii])
        return base[ii] - to_score(r)

    for b in BETAS:
        R[f"dPORI_beta{b}"] = dpori(b, REDUCTION)
    R["dPORI_use_zeroed_beta1"] = dpori(1.0, 1.0)
    key = f"dPORI_beta{BETAS[-1]}"
    # Outcome-specific scenarios: cutting class c in group g by REDUCTION lowers the ARG load of class c by
    # beta x REDUCTION x share(g, c). Two outcomes are tracked separately: (1) the TOTAL resistome and
    # (2) the CRITICAL ("dangerous") resistome. Absolute reductions are in normalised-abundance units, so they
    # are comparable among stations; percentages are relative to the station's own load.
    f1 = REDUCTION * R.share.values                       # beta = 1 (upper bound); other betas scale linearly
    R["dTot_abs"] = f1 * A_tot[ii, ci_]
    R["dCrit_abs"] = f1 * A_crit[ii, ci_]
    R["dTot_pct"] = 100 * R.dTot_abs / tot0[ii]
    R["dCrit_pct"] = 100 * R.dCrit_abs / np.where(crit0[ii] > 0, crit0[ii], np.nan)
    R["Location"], R["Watershed"] = loc[ii], wsh[ii]
    R["PORI (0-100)"] = score[ii].round(1)
    R["Share_of_class_use_%"] = (100 * R.share).round(1)
    gate = S1["gate"]
    R["Evidence_class_consumption_vs_ARG"] = R.Class.map(lambda c: "Supported" if gate[c] else "Not supported")
    gg = S1.get("gene_gate", {})
    R["Evidence_gene_level"] = R.Class.map(lambda c: f"Supported ({gg[c]} gene(s))" if gg.get(c, 0) > 0 else "Not supported")

    best = R.sort_values(key, ascending=False).groupby("i").head(3).copy()
    best["Target_rank"] = best.groupby("i")[key].rank(ascending=False, method="first").astype(int)
    top1 = best[best.Target_rank == 1].set_index("i").reindex(range(n))
    ami_raw = top1[key].values
    ami = ami_raw / np.nanmax(ami_raw) * 100

    # source quadrant (PORI tier >= High  x  human-source signal above the median)
    O = S4["origin"]
    hs = O["Human_share_of_host_associated_%"].values
    hi_human = hs >= np.median(hs)
    hi_pori = score >= TIER_EDGES[1]
    quad, act = [], []
    for a, b in zip(hi_pori, hi_human):
        if a and b:
            quad.append("High PORI + human-source signal")
            act.append("Priority 1: reduce urban wastewater inputs (treatment, river discharges)")
        elif a:
            quad.append("High PORI + non-human signal")
            act.append("Priority 2: identify the source (livestock, other) before acting")
        elif b:
            quad.append("Low PORI + human-source signal")
            act.append("Surveillance: human footprint without elevated load")
        else:
            quad.append("Low PORI + non-human signal")
            act.append("No priority action")

    # best target per station for each outcome
    outc = {}
    for tag, absc, pctc in (("total", "dTot_abs", "dTot_pct"), ("critical", "dCrit_abs", "dCrit_pct")):
        bt = R.sort_values(absc, ascending=False).groupby("i").head(1).set_index("i").reindex(range(n))
        vals = bt[absc].fillna(0).values
        mx = vals.max() if vals.max() > 0 else 1.0
        outc[tag] = dict(bt=bt, abs=vals, ami=vals / mx * 100, pct=bt[pctc].values)

    T = D["meta"].copy()
    T["PORI (0-100)"] = score.round(1)
    T["PORI tier"] = tier(score)
    T["Target_group"] = top1["Group"].values
    T["Target_class"] = top1["Class"].values
    T["Share_of_class_use_%"] = top1["Share_of_class_use_%"].values
    for b in BETAS:
        T[f"dPORI_beta{b}_cut{int(REDUCTION * 100)}%"] = top1[f"dPORI_beta{b}"].round(2).values
    T["dPORI_if_use_zeroed_beta1"] = top1["dPORI_use_zeroed_beta1"].round(2).values
    T["Ceiling_dPORI_any_ARG_action"] = ceiling.round(1)
    T["AMI (0-100)"] = ami.round(1)
    T["AMI tier"] = tier(ami)
    for tag, o in outc.items():
        bt = o["bt"]
        none = o["abs"] <= 0
        T[f"Target_group_{tag}"] = np.where(none, "none", bt["Group"].values)
        T[f"Target_class_{tag}"] = np.where(none, "none", bt["Class"].values)
        T[f"Share_of_class_use_%_{tag}"] = bt["Share_of_class_use_%"].values
        for b in BETAS:
            T[f"Reduction_{tag}_ARG_% beta{b}"] = (b * o["pct"]).round(2)
        T[f"AMI_{tag} (0-100)"] = o["ami"].round(1)
        T[f"AMI_{tag} tier"] = tier(o["ami"])
        T[f"Evidence_{tag}_target_class"] = bt["Evidence_class_consumption_vs_ARG"].values
    T["Evidence_for_target_class"] = top1["Evidence_class_consumption_vs_ARG"].values
    T["Evidence_for_target_class_gene_level"] = top1["Evidence_gene_level"].values
    T["Source quadrant"] = quad
    T["Suggested action (source)"] = act
    T = T.sort_values("AMI (0-100)", ascending=False).reset_index(drop=True)
    T.insert(0, "AMI rank", np.arange(1, n + 1))

    top3 = best.sort_values(["Target_rank", key], ascending=[True, False])[
        ["Target_rank", "Location", "Watershed", "PORI (0-100)", "Group", "Class", "Share_of_class_use_%",
         f"dPORI_beta{BETAS[0]}", f"dPORI_beta{BETAS[1]}", key, "Evidence_class_consumption_vs_ARG",
         "Evidence_gene_level"]].round(2)
    by_group = R.groupby("Group")[key].sum().div(n).sort_values(ascending=False).round(2).rename("Mean dPORI (sum over classes)").reset_index()
    by_class = R.groupby("Class")[key].sum().div(n).sort_values(ascending=False).round(2).rename("Mean dPORI (sum over groups)").reset_index()
    stab = np.mean([g.sort_values(f"dPORI_beta{BETAS[0]}").iloc[-1].name == g.sort_values(key).iloc[-1].name for _, g in R.groupby("i")])
    log(f"[step6] AMI tiers: {pd.Series(tier(ami)).value_counts().to_dict()}; top target identical for beta={BETAS[0]} and {BETAS[-1]} "
        f"in {100 * stab:.0f}% of stations; median dPORI of the top target = {np.nanmedian(ami_raw):.2f} points "
        f"(median PORI = {np.median(score):.1f}); evidence-supported classes: {sum(gate.values())}/{len(gate)}")
    log(f"[step6] target groups: {T['Target_group'].value_counts().to_dict()}; target classes: {T['Target_class'].value_counts().to_dict()}")
    for tag in outc:
        log(f"[step6] outcome '{tag}': AMI tiers {pd.Series(tier(outc[tag]['ami'])).value_counts().to_dict()}; "
            f"target groups {T[f'Target_group_{tag}'].value_counts().to_dict()}; "
            f"target classes {T[f'Target_class_{tag}'].value_counts().to_dict()}")
    same = ((T.Target_group_total == T.Target_group_critical) & (T.Target_class_total == T.Target_class_critical)).mean()
    log(f"[step6] same best target for total and critical resistome in {100 * same:.0f}% of stations")
    by_out = pd.concat([R.groupby("Class")[["dTot_pct", "dCrit_pct"]].mean().add_prefix("mean_").round(3)
                        .assign(Level="class").reset_index().rename(columns={"Class": "Name"}),
                        R.groupby("Group")[["dTot_pct", "dCrit_pct"]].mean().add_prefix("mean_").round(3)
                        .assign(Level="group").reset_index().rename(columns={"Group": "Name"})])
    log(f"[step6] source quadrants: {pd.Series(quad).value_counts().to_dict()}")
    return dict(table=T, top3=top3, by_group=by_group, by_class=by_class, by_outcome=by_out)


def step6s_ami_source(D, C, AB, P, S4):
    """
    SOURCE-BASED AMI (main AMI). The only driver that the data support is the anthropogenic signal (share of reads from taxa found only
    in humans/their waste; Step 4): it tracks total ARGs, high-alert pathogens and PORI across stations, robustly to salinity, substrate,
    basin aggregation and leave-one-out. Consumption, population and livestock abundance are not associated with the resistome (Step 1).
      1. Elasticity of each PORI component to the anthropogenic signal: OLS of log1p(component) on log1p(signal, % of reads) + salinity (z), with a
         basin-cluster bootstrap. A component responds in the scenario ONLY if its bootstrap 95% CI excludes zero (data-supported).
      2. Scenario: the anthropogenic signal is brought down to the reference level (SOURCE_REF_Q quantile of the stations) with a maximum cut of
         SOURCE_CUT (e.g. wastewater treatment upgrade / source control); stations already below the reference do not change. Supported components
         fall by slope x change in log1p(signal); the other components are unchanged. The PORI is recomputed.
      3. AMI = PORI points avoided, rescaled 0-100 across stations; 95% interval from the bootstrap slopes.
    Context (not part of the AMI): the two groups, humans included, with the largest share of antibiotic use in the basin (surveillance only).
    The relationship is cross-sectional; the scenario holds only if the signal is causally linked to the resistome.
    """
    n = len(D["meta"]); loc, wsh = D["meta"]["Location"].values, D["meta"]["Watershed"].values
    raw, score = C["raw"], P["score"]; o = S4["origin"]
    sal = pd.to_numeric(D["phys"]["Salinity"], errors="coerce").values; sal = (sal - np.nanmean(sal)) / np.nanstd(sal)
    sig = o["Human_only_taxa_%"].values
    comps = ["ARGs_Total", "Critical_ARGs", "Pathogens_Total", "Alert_Pathogens", "Mobility"]
    y = {c: raw[c].values.astype(float) for c in comps}
    lx = np.log1p(sig); X = np.column_stack([np.ones(n), lx, sal])
    ref = float(np.quantile(sig, SOURCE_REF_Q))
    def fit(idx):
        out = {}
        for c in comps:
            b = np.linalg.lstsq(X[idx], np.log1p(y[c][idx]), rcond=None)[0]; out[c] = b[1]
        return out
    b0 = fit(np.arange(n))
    rng = np.random.default_rng(SEED); ws = np.unique(wsh); boots = []
    for _ in range(SOURCE_BOOT):
        pick = rng.choice(ws, len(ws), replace=True)
        idx = np.concatenate([np.where(wsh == w)[0] for w in pick])
        if len(np.unique(lx[idx])) > 3: boots.append(fit(idx))
    B = pd.DataFrame(boots)
    el = pd.DataFrame({"Component": comps, "Slope (log1p component per log1p signal)": [b0[c] for c in comps],
                       "CI low (2.5%)": [B[c].quantile(.025) for c in comps], "CI high (97.5%)": [B[c].quantile(.975) for c in comps]})
    el["Data-supported (CI > 0)"] = np.where(el["CI low (2.5%)"] > 0, "Yes", "No")
    el["Median change at the scenario cut (%)"] = [round(float(np.nanmedian(100 * (np.expm1(np.log1p(y[c]) + max(b0[c], 0) * (np.log1p(np.where(sig > ref, np.maximum(ref, sig * (1 - SOURCE_CUT)), sig)) - lx)) / np.where(y[c] > 0, y[c], np.nan) - 1))), 1) if b0[c] > 0 else 0.0 for c in comps]
    sup = {c: bool(el.loc[el.Component == c, "CI low (2.5%)"].iloc[0] > 0) for c in comps}
    # PORI recomputation
    lt, lc, lp, lq = [np.log1p(raw[c].values) for c in ("ARGs_Total", "Critical_ARGs", "Pathogens_Total", "Alert_Pathogens")]
    mob0 = raw["Mobility"].values; mmin, mmax = mob0.min(), mob0.max()
    pil = lambda x, lo, hi: np.clip((x - lo) / (hi - lo), 0, None)
    sq = lambda a, b: np.sqrt((a + EPS) * (b + EPS))
    def pori_raw(t, cr, pt_, pa_, mb):
        return (sq(pil(np.log1p(t), lt.min(), lt.max()), pil(np.log1p(cr), lc.min(), lc.max())) *
                sq(pil(np.log1p(pt_), lp.min(), lp.max()), pil(np.log1p(pa_), lq.min(), lq.max())) * (1 + np.clip((mb - mmin) / (mmax - mmin), 0, None)))
    r0 = pori_raw(y["ARGs_Total"], y["Critical_ARGs"], y["Pathogens_Total"], y["Alert_Pathogens"], mob0)
    lo, hi = r0.min(), r0.max(); to = lambda r: (r - lo) / (hi - lo) * 100
    assert np.allclose(to(r0), score, atol=1e-6), "baseline must reproduce the observed PORI"
    def scenario(slopes, cut, ref_=None):
        r_ = ref if ref_ is None else ref_
        dls = np.log1p(np.maximum(np.minimum(sig, r_) if False else np.maximum(r_, sig * (1 - cut)), 0) * (sig > r_) + sig * (sig <= r_)) - lx
        new = {c: np.expm1(np.log1p(y[c]) + slopes[c] * dls) if slopes.get(c, 0) > 0 else y[c] for c in comps}
        return score - to(pori_raw(new["ARGs_Total"], new["Critical_ARGs"], new["Pathogens_Total"], new["Alert_Pathogens"], new["Mobility"]))
    slopes = {c: (b0[c] if sup[c] else 0.0) for c in comps}
    d = scenario(slopes, SOURCE_CUT)
    dd = np.array([scenario({c: max(bt[c], 0.0) for c in comps}, SOURCE_CUT) for bt in boots[:500]])
    sens = {f"dPORI_maxcut_{int(k * 100)}%": scenario(slopes, k) for k in (0.3, 0.5, 0.7)}
    sens["dPORI_ref_median"] = scenario(slopes, SOURCE_CUT, float(np.median(sig)))
    ami = d / d.max() * 100
    # indicative action focus from the composition of the anthropogenic signal
    waste, gut = o["Waste_taxa_%"].values, o["Gut_taxa_%"].values
    wm, gm_ = np.median(waste), np.median(gut)
    act = []
    for i in range(n):
        if d[i] < 0.5: act.append("Monitor")
        elif waste[i] >= wm and gut[i] >= gm_: act.append("WWTP upgrade + sewer control")
        elif waste[i] >= wm: act.append("WWTP upgrade")
        elif gut[i] >= gm_: act.append("Faecal input control")
        else: act.append("Source tracking first")
    # antibiotic-use context: the two groups (humans included) with the largest share of total antibiotic use in the basin
    classes = list(AB_CLASSES); gu = {g: sum(AB[c][g] for c in classes) for g in GROUPS}; gtot = sum(gu.values())
    ctx, ctx_live = [], []
    animal_groups = [g for g in GROUPS if g != "Humans"]
    for i in range(n):
        order = sorted(GROUPS, key=lambda g: -gu[g][i]); sh = lambda g: 100 * gu[g][i] / gtot[i] if gtot[i] > 0 else 0
        ctx.append(f"{order[0]} ({sh(order[0]):.0f}%), {order[1]} ({sh(order[1]):.0f}%)")
        lv = max(animal_groups, key=lambda g: gu[g][i]); ctx_live.append(f"{lv} ({sh(lv):.0f}% of total use)")
    T = D["meta"].copy(); T["PORI (0-100)"] = score.round(1); T["PORI tier"] = tier(score)
    T["Anthropogenic_taxa_%"] = sig.round(1)
    T["PORI_points_avoided"] = d.round(2); T["PORI_points_avoided_CI_low"] = np.quantile(dd, .025, axis=0).round(2); T["PORI_points_avoided_CI_high"] = np.quantile(dd, .975, axis=0).round(2)
    T["AMI (0-100)"] = ami.round(1); T["AMI tier"] = tier(ami); T["Action focus"] = act; T["Antibiotic-use context (surveillance)"] = ctx; T["Top livestock group"] = ctx_live
    for k, v in sens.items(): T[k] = np.round(v, 2)
    T["PORI_after_%"] = (100 * (1 - d / np.where(score > 0, score, np.nan))).round(0)
    T = T.sort_values("AMI (0-100)", ascending=False).reset_index(drop=True); T.insert(0, "AMI rank", np.arange(1, n + 1))
    log("[step6s] source-based AMI: supported components = " + ", ".join(c for c in comps if sup[c]) + f"; reference signal = {ref:.1f}% of reads; stations with a gain = {int((d > 0.005).sum())}; median PORI points avoided (all) = {np.median(d):.2f}, max = {d.max():.2f}")
    return dict(table=T, elasticity=el)


def step6v_ami_validation(D, C, P, S4):
    """
    OUT-OF-SAMPLE CHECK of the source-based AMI (no new data; it does not test an intervention).
      A. Leave-one-watershed-out (LOWO): elasticities fitted on 16 watersheds predict the stations of the held-out one; skill = out-of-sample R2
         against a salinity-only model, with a permutation test (signal shuffled among stations, VALID_PERM draws).
      B. Within-watershed slopes (watershed fixed effects): do differences in the signal between stations of the same watershed track the components?
      C. Pseudo-intervention: for pairs of stations of one watershed, the PORI drop predicted by setting the signal of the higher station to that
         of the lower one (slopes fitted without that watershed) is compared with the observed PORI gap.
    """
    n = len(D["meta"]); wsh = D["meta"]["Watershed"].values; raw, score = C["raw"], P["score"]; o = S4["origin"]
    loc = D["meta"]["Location"].values
    sal = pd.to_numeric(D["phys"]["Salinity"], errors="coerce").values; sal = (sal - np.nanmean(sal)) / np.nanstd(sal)
    sig = o["Human_only_taxa_%"].values; lx = np.log1p(sig)
    comps = ["ARGs_Total", "Critical_ARGs", "Pathogens_Total", "Alert_Pathogens", "Mobility"]
    y = {c: raw[c].values.astype(float) for c in comps}; ly = {c: np.log1p(y[c]) for c in comps}
    X = np.column_stack([np.ones(n), lx, sal]); X0 = np.column_stack([np.ones(n), sal])
    ws = np.unique(wsh); rng = np.random.default_rng(SEED)
    lt, lc, lp, lq = [np.log1p(raw[c].values) for c in ("ARGs_Total", "Critical_ARGs", "Pathogens_Total", "Alert_Pathogens")]
    mob0 = raw["Mobility"].values; mmin, mmax = mob0.min(), mob0.max()
    pil = lambda x, lo, hi: np.clip((x - lo) / (hi - lo), 0, None)
    sq = lambda a, b: np.sqrt((a + EPS) * (b + EPS))
    def pori_raw(t, cr, pt_, pa_, mb):
        return (sq(pil(np.log1p(t), lt.min(), lt.max()), pil(np.log1p(cr), lc.min(), lc.max())) *
                sq(pil(np.log1p(pt_), lp.min(), lp.max()), pil(np.log1p(pa_), lq.min(), lq.max())) * (1 + np.clip((mb - mmin) / (mmax - mmin), 0, None)))
    r0 = pori_raw(y["ARGs_Total"], y["Critical_ARGs"], y["Pathogens_Total"], y["Alert_Pathogens"], mob0)
    lo, hi = r0.min(), r0.max(); to = lambda r: (r - lo) / (hi - lo) * 100
    assert np.allclose(to(r0), score, atol=1e-6), "baseline must reproduce the observed PORI"
    def lowo(lxv):
        Xv = np.column_stack([np.ones(n), lxv, sal]); pred = {c: np.zeros(n) for c in comps}; null = {c: np.zeros(n) for c in comps}
        for w in ws:
            te = wsh == w; tr = ~te
            for c in comps:
                pred[c][te] = Xv[te] @ np.linalg.lstsq(Xv[tr], ly[c][tr], rcond=None)[0]
                null[c][te] = X0[te] @ np.linalg.lstsq(X0[tr], ly[c][tr], rcond=None)[0]
        return pred, null
    oos = lambda obs, pr, nl: 1 - np.sum((obs - pr) ** 2) / np.sum((obs - nl) ** 2)
    pred, null = lowo(lx); perm = {c: [] for c in comps}
    for _ in range(VALID_PERM):
        pp, nn_ = lowo(rng.permutation(lx))
        for c in comps: perm[c].append(oos(ly[c], pp[c], nn_[c]))
    A = pd.DataFrame([(c, oos(ly[c], pred[c], null[c]), spearmanr(ly[c], pred[c])[0],
                       (np.sum(np.array(perm[c]) >= oos(ly[c], pred[c], null[c])) + 1) / (VALID_PERM + 1)) for c in comps],
                     columns=["Component", "Out-of-sample R2 vs salinity-only", "Spearman (predicted, observed)", "Permutation p"])
    pori_from = lambda pr: to(pori_raw(*[np.expm1(pr[c]) for c in comps]))
    pm, pn = pori_from(pred), pori_from(null)
    PORI_ROW = pd.DataFrame([dict(Component="PORI (from predicted components)", **{"Out-of-sample R2 vs salinity-only": 1 - np.sum((score - pm) ** 2) / np.sum((score - pn) ** 2),
                                  "Spearman (predicted, observed)": spearmanr(score, pm)[0], "Permutation p": np.nan, "Spearman salinity-only": spearmanr(score, pn)[0]})])
    A = pd.concat([A, PORI_ROW], ignore_index=True)
    def within(idx_ws):
        keep = np.concatenate([np.where(wsh == w)[0] for w in idx_ws])
        d = pd.DataFrame({"w": wsh[keep], "lx": lx[keep], "sal": sal[keep], **{c: ly[c][keep] for c in comps}})
        d = d[d.groupby("w")["lx"].transform("count") > 1]
        for col in ["lx", "sal"] + comps: d[col] = d[col] - d.groupby("w")[col].transform("mean")
        return {c: np.linalg.lstsq(d[["lx", "sal"]].values, d[c].values, rcond=None)[0][0] for c in comps}
    w0 = within(ws); bw = pd.DataFrame([within(rng.choice(ws, len(ws), replace=True)) for _ in range(SOURCE_BOOT)])
    pooled = {c: np.linalg.lstsq(X, ly[c], rcond=None)[0][1] for c in comps}
    B = pd.DataFrame({"Component": comps, "Pooled slope": [pooled[c] for c in comps], "Within-watershed slope": [w0[c] for c in comps],
                      "CI low (2.5%)": [bw[c].quantile(.025) for c in comps], "CI high (97.5%)": [bw[c].quantile(.975) for c in comps]})
    cnt = pd.Series(wsh).value_counts(); rows = []
    for w in cnt[cnt > 1].index:
        idx = np.where(wsh == w)[0]; tr = wsh != w
        sl = {c: np.linalg.lstsq(X[tr], ly[c][tr], rcond=None)[0][1] for c in comps}
        for a in idx:
            for b in idx:
                if sig[a] > sig[b]:
                    dls = np.log1p(sig[b]) - lx[a]
                    new = {c: (np.expm1(ly[c][a] + sl[c] * dls) if sl[c] > 0 else y[c][a]) for c in comps}
                    pn_ = to(pori_raw(*[new[c] for c in comps]))
                    rows.append((w, loc[a], loc[b], score[a] - pn_, score[a] - score[b], sig[a] - sig[b]))
    Cp = pd.DataFrame(rows, columns=["Watershed", "Higher-signal station", "Lower-signal station", "Predicted PORI drop", "Observed PORI gap", "Signal gap (% of reads)"])
    Cw = Cp.loc[Cp.groupby("Watershed")["Signal gap (% of reads)"].idxmax()]
    rho_all = spearmanr(Cp["Predicted PORI drop"], Cp["Observed PORI gap"])[0]; rho_w = spearmanr(Cw["Predicted PORI drop"], Cw["Observed PORI gap"])[0]
    agree = int((Cw["Observed PORI gap"] > 0).sum())
    log("[step6v] LOWO out-of-sample R2 vs salinity-only: " + "; ".join(f"{a_} {b_:.2f}" for a_, b_ in zip(A["Component"], A["Out-of-sample R2 vs salinity-only"])) +
        f"; pseudo-intervention pairs = {len(Cp)} in {Cp.Watershed.nunique()} watersheds, Spearman all pairs = {rho_all:.2f}, one pair per watershed = {rho_w:.2f}, "
        f"higher-signal station has the higher PORI in {agree}/{len(Cw)} watersheds")
    S = pd.DataFrame([dict(Metric="Watersheds with >= 2 stations", Value=int(len(cnt[cnt > 1]))), dict(Metric="Pseudo-intervention pairs", Value=int(len(Cp))),
                      dict(Metric="Spearman, all pairs", Value=round(rho_all, 2)), dict(Metric="Spearman, one pair per watershed", Value=round(rho_w, 2)),
                      dict(Metric="Watersheds where the higher-signal station has the higher PORI", Value=f"{agree}/{len(Cw)}")])
    return dict(lowo=A, within=B, pairs=Cp, summary=S)


def step6c_multidrug_ami(D, C, AB, S6):
    """
    Supplementary AMI for MULTIDRUG (multi-class) genes. These genes respond to the use of any of
    several classes, so naming one class to cut would be arbitrary. The target is therefore a host
    group only: the group whose antibiotic use, across the classes affecting each gene, contributes most.
    For gene j (classes S_j) and group g, a 50% cut of g's use lowers the gene's load by
    beta x REDUCTION x mean_{c in S_j} share(g, c)   (beta = 1, upper bound; each class assumed equally
    likely to be the selecting one). Reductions are summed over genes, per station, in normalised-abundance
    units; AMI_multidrug = best group's reduction rescaled 0-100 across stations.
    """
    genes, GM = C["genes"], C["GM"]
    classes, n = list(AB_CLASSES), len(D["meta"])
    gcl = {g: classes_of(D["gene_class_text"].get(g, "")) for g in genes}
    multi = np.array([len(gcl[g]) > 1 for g in genes])
    shares = {g: np.zeros((n, len(classes))) for g in GROUPS}
    for ci, c in enumerate(classes):
        tc = np.asarray(AB[c]["TOT"], float)
        for g in GROUPS:
            with np.errstate(invalid="ignore", divide="ignore"):
                shares[g][:, ci] = np.where(tc > 0, np.asarray(AB[c][g], float) / tc, 0.0)
    cidx = {c: k for k, c in enumerate(classes)}
    # mean share of group g over the classes of each multidrug gene: matrix genes x stations per group
    red = {}
    for g in GROUPS:
        W = np.zeros((n, len(genes)))
        for j, gene in enumerate(genes):
            if multi[j]:
                W[:, j] = np.mean([shares[g][:, cidx[c]] for c in gcl[gene]], axis=0)
        red[g] = REDUCTION * (GM * W).sum(axis=1)          # absolute reduction of the multidrug load
    Rm = pd.DataFrame(red)
    best_g = Rm.idxmax(axis=1).values
    best_v = Rm.max(axis=1).values
    ami = best_v / best_v.max() * 100
    tot = GM.sum(axis=1)
    md_load = GM[:, multi].sum(axis=1)
    T = D["meta"][["Location", "Watershed"]].copy()
    T["PORI (0-100)"] = S6["table"].set_index("Location").loc[T.Location, "PORI (0-100)"].values
    T["PORI tier"] = S6["table"].set_index("Location").loc[T.Location, "PORI tier"].values
    T["Multidrug_share_of_ARG_load_%"] = (100 * md_load / tot).round(1)
    T["Target_group_multidrug"] = best_g
    T["Reduction_multidrug_load_%_beta1"] = (100 * best_v / np.where(md_load > 0, md_load, np.nan)).round(1)
    T["AMI_multidrug (0-100)"] = ami.round(1)
    T["AMI_multidrug tier"] = tier(ami)
    T = T.sort_values("AMI_multidrug (0-100)", ascending=False).reset_index(drop=True)
    T.insert(0, "Rank", np.arange(1, n + 1))
    log(f"[step6c] multidrug genes = {int(multi.sum())}, median share of ARG load = {np.median(100 * md_load / tot):.0f}%; "
        f"target groups: {T['Target_group_multidrug'].value_counts().to_dict()}; tiers: {T['AMI_multidrug tier'].value_counts().to_dict()}")
    return dict(table=T)


def step6b_assignment_sensitivity(D, C, AB):
    """
    How much do the AMI targets depend on how multi-class genes are assigned to antibiotic classes?
      single  : only genes acting on exactly one class are used (default of step 6)
      split   : a gene acting on k classes counts 1/k for each
      full    : a gene counts fully for each of its classes (no split)
    For each station and outcome (total / critical resistome) the best (group, class) target is
    recomputed and compared with the default.
    """
    genes, GM = C["genes"], C["GM"]
    classes, n = list(AB_CLASSES), len(D["meta"])
    gcl = {g: classes_of(D["gene_class_text"].get(g, "")) for g in genes}
    crit = np.zeros(len(genes), bool)
    crit[C["crit_idx"]] = True
    variants = {"split": lambda S, c: (1 / len(S) if c in S else 0.0),
                "single": lambda S, c: (1.0 if (len(S) == 1 and c in S) else 0.0),
                "full": lambda S, c: (1.0 if c in S else 0.0)}
    loc = D["meta"]["Location"].values
    rows, res = [], {}
    for vk, f in variants.items():
        Am = {"total": np.zeros((n, len(classes))), "critical": np.zeros((n, len(classes)))}
        for ci, c in enumerate(classes):
            w = np.array([f(gcl[g], c) for g in genes])
            Am["total"][:, ci] = GM @ w
            Am["critical"][:, ci] = GM[:, crit] @ w[crit]
        for tag in Am:
            for i in range(n):
                bv, bg, bc = 0.0, "none", "none"
                for ci, c in enumerate(classes):
                    tc = AB[c]["TOT"][i]
                    if tc <= 0:
                        continue
                    for g in GROUPS:
                        v = REDUCTION * AB[c][g][i] / tc * Am[tag][i, ci]
                        if v > bv:
                            bv, bg, bc = v, g, c
                res[(vk, tag, i)] = (bg, bc)
    T = pd.DataFrame({"Location": loc})
    for tag in ("total", "critical"):
        for vk in variants:
            T[f"{tag}_{vk}_group"] = [res[(vk, tag, i)][0] for i in range(n)]
            T[f"{tag}_{vk}_class"] = [res[(vk, tag, i)][1] for i in range(n)]
    srows = []
    for tag in ("total", "critical"):
        for vk in ("split", "full"):
            same_gc = np.mean([res[("single", tag, i)] == res[(vk, tag, i)] for i in range(n)])
            same_c = np.mean([res[("single", tag, i)][1] == res[(vk, tag, i)][1] for i in range(n)])
            top = pd.Series([res[(vk, tag, i)][1] for i in range(n)]).value_counts().to_dict()
            srows.append(dict(Outcome=tag, Assignment=vk, Same_group_and_class_as_default_pct=round(100 * same_gc),
                              Same_class_pct=round(100 * same_c), Target_classes=str(top)))
    SM = pd.DataFrame(srows)
    # composition of the fluoroquinolone class load
    fq = [g for g in genes if "Fluoroquinolones" in gcl[g]]
    ab_fq = np.array([GM[:, genes.index(g)].sum() for g in fq])
    multi_share = 100 * sum(a for a, g in zip(ab_fq, fq) if len(gcl[g]) > 1) / ab_fq.sum()
    log("[step6b] target stability under alternative multi-class assignment: " +
        "; ".join(f"{r.Outcome}/{r.Assignment}: {r.Same_group_and_class_as_default_pct}% same target" for r in SM.itertuples()))
    log(f"[step6b] fluoroquinolone class: {len(fq)} genes, {multi_share:.0f}% of their abundance from multi-class genes")
    return dict(stations=T, summary=SM, fq_multi_share=multi_share)


# ======================================================================
# OUTPUT: EXCEL
# ======================================================================
def write_excel(path, sheets, readme):
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter
    with pd.ExcelWriter(path, engine="openpyxl") as w:
        pd.DataFrame({"README": readme}).to_excel(w, sheet_name="README", index=False)
        for name, df in sheets.items():
            d = df.copy()
            for c in d.columns:
                if pd.api.types.is_float_dtype(d[c]):
                    d[c] = d[c].round(4)
            d.to_excel(w, sheet_name=name[:31], index=False)
        for ws in w.book.worksheets:
            for c in ws[1]:
                c.font = Font(name="Arial", bold=True, color="FFFFFF")
                c.fill = PatternFill("solid", fgColor="1F3864")
                c.alignment = Alignment(wrap_text=True, vertical="center")
            for row in ws.iter_rows(min_row=2):
                for c in row:
                    c.font = Font(name="Arial", size=10)
            for col in ws.columns:
                L = max(len(str(x.value or "")) for x in col[:60])
                ws.column_dimensions[get_column_letter(col[0].column)].width = min(max(L + 2, 11), 120 if ws.title == "README" else 34)
            ws.freeze_panes = "A2"


# ======================================================================
# OUTPUT: FIGURES
# ======================================================================
def set_style():
    plt.rcParams.update({
        "figure.facecolor": PAL["surface"], "axes.facecolor": PAL["surface"], "savefig.facecolor": PAL["surface"],
        "text.color": PAL["ink"], "axes.labelcolor": PAL["ink2"], "xtick.color": PAL["ink2"], "ytick.color": PAL["ink2"],
        "axes.edgecolor": PAL["grid"], "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": PAL["grid"], "grid.linewidth": 0.6, "axes.axisbelow": True,
        "font.size": 8.5, "axes.titlesize": 9.5, "axes.titleweight": "bold", "legend.frameon": False,
    })


def tier_legend(ax, loc="lower right", title=None):
    ax.legend(handles=[Patch(facecolor=TIER_COL[t], label=t) for t in TIER_LABELS[::-1]],
              loc=loc, title=title, fontsize=7.5)


def make_figures(fdir, D, C, S1, S2, ZF, P, S4, S6):
    os.makedirs(fdir, exist_ok=True)
    set_style()
    n = len(D["meta"])

    # ---- Fig. 1: correlation heatmaps (host-group population vs antibiotic-class use) ----
    Hm = S1["H"]
    cmap_d = LinearSegmentedColormap.from_list("div", DIV_RAMP)
    panels = [("Total", "a  Total population (incl. humans) vs total use"),
              ("Human use", "b  Human population vs human use"),
              ("Veterinary use", "c  Livestock and pets vs veterinary use")]
    fig, axs = plt.subplots(3, 1, figsize=(9, 11), gridspec_kw={"height_ratios": [10, 1.6, 9]}, constrained_layout=True)
    for ax, (pn, ttl) in zip(axs, panels):
        d = Hm[Hm.Panel == pn]
        rho = d.pivot(index="Population", columns="Class", values="rho_basin").loc[d.Population.unique(), list(AB_CLASSES)]
        pv = d.pivot(index="Population", columns="Class", values="p_basin").loc[rho.index, rho.columns]
        im = ax.imshow(rho.values, cmap=cmap_d, vmin=-1, vmax=1, aspect="auto")
        ax.grid(False)
        ax.set_xticks(range(rho.shape[1])); ax.set_yticks(range(rho.shape[0])); ax.set_yticklabels(rho.index)
        ax.set_xticklabels(rho.columns if pn == "Veterinary use" else [""] * rho.shape[1], rotation=40, ha="right")
        ax.set_title(ttl, loc="left")
        for (i, j), v in np.ndenumerate(rho.values):
            star = "*" if pv.values[i, j] < 0.05 else ""
            ax.text(j, i, f"{v:.2f}{star}", ha="center", va="center", fontsize=7, color="white" if abs(v) > 0.72 else PAL["ink"])
    fig.colorbar(im, ax=axs, shrink=0.5, label="Spearman rho (watershed level, n = 17)")
    fig.savefig(os.path.join(fdir, "Fig1_correlations_heatmap.png"), dpi=200); plt.close(fig)

    # ---- Fig. S1: all consumption classes vs all ARG classes (watershed level) ----
    Mx = S1["M"]
    rho = Mx.pivot(index="Consumption", columns="ARG_class", values="rho_basin").loc[Mx.Consumption.unique(), Mx.ARG_class.unique()]
    pv = Mx.pivot(index="Consumption", columns="ARG_class", values="p_basin").loc[rho.index, rho.columns]
    qv = Mx.pivot(index="Consumption", columns="ARG_class", values="q_basin").loc[rho.index, rho.columns]
    fig, ax = plt.subplots(figsize=(11, 6.5), constrained_layout=True)
    im = ax.imshow(rho.values, cmap=cmap_d, vmin=-1, vmax=1, aspect="auto"); ax.grid(False)
    ax.set_xticks(range(rho.shape[1])); ax.set_xticklabels(rho.columns, rotation=40, ha="right")
    ax.set_yticks(range(rho.shape[0])); ax.set_yticklabels(rho.index)
    for (i, j), v in np.ndenumerate(rho.values):
        ax.text(j, i, f"{v:.2f}" + ("*" if pv.values[i, j] < 0.05 else "") + ("\u2020" if qv.values[i, j] < 0.10 else ""),
                ha="center", va="center", fontsize=6.5, color="white" if abs(v) > 0.72 else PAL["ink"])
    ax.set_xlabel("ARG resistance class (summed normalised abundance)"); ax.set_ylabel("Total antibiotic consumption class")
    ax.set_title("Consumption classes vs ARG classes (Spearman, watershed level n = 17; * p<0.05, \u2020 FDR q<0.10)", loc="left")
    fig.colorbar(im, ax=ax, shrink=0.8, label="Spearman rho")
    fig.savefig(os.path.join(fdir, "FigS1_consumption_vs_ARG_classes.png"), dpi=200); plt.close(fig)

    # ---- Fig. 2: distribution and outliers --------------------------------
    S = S2["summary"].sort_values("Skewness")
    st = S2["stations"].sort_values("Flag_share_%")
    fig, ax = plt.subplots(1, 2, figsize=(11, 8), gridspec_kw={"width_ratios": [1.15, 1]})
    colors = [TEAL if l.startswith("Basin") else ACCENT for l in S.Level]
    ax[0].barh(S.Variable, S.Skewness, color=colors, height=0.7, edgecolor=PAL["surface"], linewidth=1)
    ax[0].set_xlabel("Skewness"); ax[0].set_title("a  Distribution shape of screened variables")
    ax[0].tick_params(axis="y", labelsize=6.5)
    ax[0].legend(handles=[Patch(facecolor=TEAL, label="Basin level (n=17)"),
                          Patch(facecolor=ACCENT, label="Station level (n=36)")], loc="lower right")
    ax[1].barh(st.Location, st["Flag_share_%"], color=ACCENT, height=0.7, edgecolor=PAL["surface"], linewidth=1)
    ax[1].set_xlabel("Variables flagged as Tukey outliers (% of screened)"); ax[1].set_title("b  Outlier burden per station")
    ax[1].tick_params(axis="y", labelsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(fdir, "FigS2_distribution_outliers.png"), dpi=200); plt.close(fig)

    # ---- Fig. 3: PORI bars + component heatmap ----------------------------
    T = P["table"]
    comp_cols = ["Pathogens_Total (0-100)", "Alert_Pathogens (0-100)", "ARGs_Total (0-100)",
                 "Critical_ARGs (0-100)", "Mobility (0-100)"]
    lab = ["Pathogens", "High-alert\npathogens", "ARGs", "Critical\nARGs", "Mobility"]
    order = T.index[::-1]
    fig, ax = plt.subplots(1, 2, figsize=(11, 9), gridspec_kw={"width_ratios": [1.1, 1]}, sharey=True)
    ax[0].barh(range(n), T.loc[order, "PORI (0-100)"], color=[TIER_COL[t] for t in T.loc[order, "PORI tier"]],
               height=0.72, edgecolor=PAL["surface"], linewidth=1)
    ax[0].set_yticks(range(n)); ax[0].set_yticklabels(T.loc[order, "Location"], fontsize=7.5)
    ax[0].set_xlabel("PORI (0-100, relative to the 36 stations)"); ax[0].set_title("a  PORI by station")
    for v in TIER_EDGES:
        ax[0].axvline(v, color=PAL["ink2"], lw=0.6, ls=":")
    tier_legend(ax[0], "center right")
    cmap = LinearSegmentedColormap.from_list("seq", SEQ_RAMP)
    M = T.loc[order, comp_cols].values.astype(float)
    ax[1].imshow(M, aspect="auto", cmap=cmap, vmin=0, vmax=100)
    ax[1].set_xticks(range(5)); ax[1].set_xticklabels(lab, fontsize=7.5); ax[1].grid(False)
    for (i, j), v in np.ndenumerate(M):
        ax[1].text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6, color="white" if v > 62 else PAL["ink"])
    ax[1].set_title("b  Components (0-100)")
    fig.tight_layout(); fig.savefig(os.path.join(fdir, "Fig3_PORI.png"), dpi=200); plt.close(fig)

    # ---- Fig. 4: sensitivity (rank intervals) -----------------------------
    fig, ax = plt.subplots(figsize=(7.5, 8.5))
    for k, i in enumerate(T.index[::-1]):
        ax.plot([T.Rank_p5[i], T.Rank_p95[i]], [k, k], color=PAL["ink2"], lw=1.4, solid_capstyle="round")
        ax.scatter(T.Rank[i], k, s=34, color=TIER_COL[T["PORI tier"][i]], edgecolor=PAL["surface"], linewidth=1.2, zorder=3)
    ax.set_yticks(range(n)); ax.set_yticklabels(T.loc[T.index[::-1], "Location"], fontsize=7.5)
    ax.set_xlabel("Rank (1 = highest PORI); bar = 5th-95th percentile across weight / offset / aggregation draws")
    ax.set_title("PORI sensitivity: rank stability"); tier_legend(ax, "lower right")
    fig.tight_layout(); fig.savefig(os.path.join(fdir, "Fig4_PORI_sensitivity.png"), dpi=200); plt.close(fig)

    # ---- Fig. 5: sentinels --------------------------------------------------
    SE = S4["sent"]; lg = S4["lg"]; score = P["score"]
    top = SE.head(2)
    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    for a, (_, r) in zip(ax[:2], top.iterrows()):
        a.scatter(score, lg[r.Species], s=30, color=ACCENT, edgecolor=PAL["surface"], linewidth=1)
        a.set_xlabel("PORI (0-100)"); a.set_ylabel("log(1 + reads per million)")
        a.set_title(f"{r.Species}", fontstyle="italic"); a.text(0.97, 0.95, f"rho = {r.rho:.2f}, q = {r.q_FDR:.3f}", transform=a.transAxes, va="top", ha="right", fontsize=8)
    d = SE.sort_values("rho").reset_index(drop=True)
    cols = [ACCENT if (q < FDR_Q_SENTINEL and r > 0) else PAL["grid"] for q, r in zip(d.q_FDR, d.rho)]
    ax[2].scatter(d.rho, range(len(d)), s=22, color=cols, edgecolor=PAL["ink2"], linewidth=0.4)
    ax[2].axvline(0, color=PAL["ink2"], lw=0.8); ax[2].set_yticks([]); ax[2].set_xlabel("Spearman rho with PORI")
    ax[2].set_title(f"All {len(d)} environmental candidates"); ax[2].legend(handles=[Patch(facecolor=ACCENT, label=f"q < {FDR_Q_SENTINEL}")], loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(fdir, "Fig5_sentinels.png"), dpi=200); plt.close(fig)

    # ---- Fig. 7b: composite AMI (scenario on the whole PORI) ---------------------
    A = S6["table"].sort_values("AMI (0-100)")
    fig, ax = plt.subplots(figsize=(9.5, 9))
    ax.barh(range(n), A["AMI (0-100)"], color=[TIER_COL[t] for t in A["AMI tier"]], height=0.72, edgecolor=PAL["surface"], linewidth=1)
    for k, (_, r) in enumerate(A.iterrows()):
        ax.text(r["AMI (0-100)"] + 1.2, k, f"{r.Target_group} - {r.Target_class}", va="center", fontsize=6.8, color=PAL["ink"])
    ax.set_yticks(range(n)); ax.set_yticklabels(A["Location"], fontsize=7.5); ax.set_xlim(0, 150); ax.set_xticks(range(0, 101, 20))
    ax.set_xlabel("AMI (0-100, relative); label = best mitigation target (group - antibiotic class)")
    ax.set_title("Composite AMI (scenario on the whole PORI): priority and target by station"); tier_legend(ax, "lower right")
    fig.tight_layout(); fig.savefig(os.path.join(fdir, "Fig7b_AMI_composite_PORI.png"), dpi=200); plt.close(fig)


def fig_gene_level(path, G):
    """Gene-level test: (a) rho of matched vs non-matched pairs, (b) observed vs permutation-expected p-values,
    (c) per-class median rho of matched genes and nominal hits."""
    set_style()
    F, CL = G["fam1"], G["CL"]
    rho, p, m, a = F["rho"], F["p"], F["match"], F["assigned"]
    fig, ax = plt.subplots(1, 3, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1, 1, 1.25]})
    # a) distributions
    rm, rn = rho[a[:, None] & m], rho[a[:, None] & ~m]
    parts = ax[0].violinplot([rn, rm], showmedians=True, showextrema=False)
    for b, c in zip(parts["bodies"], [PAL["grid"], ACCENT]):
        b.set_facecolor(c); b.set_edgecolor(PAL["ink2"]); b.set_alpha(0.9)
    parts["cmedians"].set_color(PAL["ink"])
    ax[0].axhline(0, color=PAL["ink2"], lw=0.7)
    ax[0].set_xticks([1, 2]); ax[0].set_xticklabels([f"cross\n(n={len(rn)})", f"matched\n(n={len(rm)})"])
    ax[0].set_ylabel("Spearman rho (gene vs class use, 17 watersheds)")
    d, pp = F["spec"]["all"]
    ax[0].set_title("a  Is the link class-specific?", loc="left")
    ax[0].text(0.03, 0.97, f"matched - cross = {d:+.3f}\nperm. p = {pp:.2f}", transform=ax[0].transAxes, fontsize=8, va="top")
    ax[0].set_ylim(-0.95, 0.95)
    # b) QQ-like: observed sorted -log10 p of positive pairs vs null expectation from uniform
    pv = np.sort(p.ravel())
    exp = (np.arange(1, len(pv) + 1) - 0.5) / len(pv)
    ax[1].scatter(-np.log10(exp), -np.log10(pv), s=5, color=ACCENT, linewidths=0)
    lim = max(-np.log10(pv.min()), -np.log10(exp.min())) + 0.2
    ax[1].plot([0, lim], [0, lim], color=PAL["ink2"], lw=0.8, ls=":")
    ax[1].set_xlabel("Expected -log10 p"); ax[1].set_ylabel("Observed -log10 p (both signs)")
    ax[1].set_title("b  All gene x class tests", loc="left")
    ax[1].text(0.03, 0.95, f"{rho.size} tests; q<0.10 with rho>0: {int(((F['q'] < FDR_Q_GATE) & (rho > 0)).sum())}\n"
               f"global max-rho perm. p = {F['p_global']:.2f}", transform=ax[1].transAxes, va="top", fontsize=8)
    # c) per class
    y = np.arange(len(CL))[::-1]
    ax[2].barh(y, CL.Median_rho_matched, color=ACCENT, height=0.55, label="median rho, matched genes")
    ax[2].scatter(CL.Median_rho_non_matched, y, color=PAL["ink"], s=18, zorder=3, label="median rho, other genes")
    ax[2].axvline(0, color=PAL["ink2"], lw=0.7)
    ax[2].set_yticks(y); ax[2].set_yticklabels(CL.Class, fontsize=8)
    for yi, k, t in zip(y, CL.Matched_with_rho_pos_p_lt_0_05, CL.Genes_matched_tested):
        ax[2].text(1.0, yi, f"{k}/{t}", transform=ax[2].get_yaxis_transform(), ha="right", va="center", fontsize=7, color=PAL["ink2"])
    ax[2].set_xlim(-0.35, 0.12)
    ax[2].set_xlabel("Median Spearman rho\n(numbers: matched genes with rho>0, p<0.05 / tested)")
    ax[2].set_title("c  By antibiotic class", loc="left"); ax[2].legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=7, frameon=False)
    fig.tight_layout(); fig.savefig(path, dpi=200); plt.close(fig)


def fig_ami_table(path, T, title, subtitle, mode, foot, g_col, a_col, t_col, c_col=None):
    """
    Colour-coded summary table (one row per station, sorted by AMI): PORI score chip, host-group target,
    antibiotic-class target (omitted for the multidrug analysis), AMI score with an inline bar, priority chip.
    mode: 'class' (columns for group and class) or 'multi' (group only + multidrug share of the ARG load).
    Data units are isotropic (1 unit = same length in x and y) so chips and dots stay round.
    """
    from matplotlib.patches import FancyBboxPatch, Rectangle
    HDR, ZEB, TRACK = "#17525e", "#f6f4ef", "#e9e5dd"
    GCOL = {"Humans": PAL["blue"], "Cattle": PAL["orange"], "Pigs": PAL["aqua"], "Poultry": PAL["yellow"]}
    T = T.sort_values(a_col, ascending=False).reset_index(drop=True)
    n = len(T)
    if mode == "class":
        cols = [("#", 0.55), ("Location", 1.9), ("Watershed", 2.0), ("PORI\nscore", 1.15), ("Mitigation\ngroup", 1.55),
                ("Mitigation target\n(antibiotic class)", 2.1), ("AMI score", 2.1), ("Priority\nlevel", 1.35)]
    else:
        cols = [("#", 0.55), ("Location", 1.9), ("Watershed", 2.0), ("PORI\nscore", 1.15), ("Mitigation\ngroup", 1.55),
                ("Multidrug share\nof ARG load", 2.1), ("AMI score", 2.1), ("Priority\nlevel", 1.35)]
    xs = np.concatenate([[0], np.cumsum([w for _, w in cols])])
    W = xs[-1]
    rh, hh = 0.38, 0.95
    fw = 12.2
    unit = fw * 0.94 / W                                   # inches per data unit (x and y)
    ax_h = (n * rh + hh) * unit
    fig_h = ax_h + 1.35 + 0.55
    fig = plt.figure(figsize=(fw, fig_h))
    ax = fig.add_axes([0.03, 0.55 / fig_h, 0.94, ax_h / fig_h])
    ax.set_xlim(0, W); ax.set_ylim(-(n * rh), hh); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0, 0), W, hh, boxstyle="round,pad=0,rounding_size=0.18", fc=HDR, ec="none", zorder=1))
    for (h, w), x0 in zip(cols, xs[:-1]):
        ax.text(x0 + w / 2, hh / 2, h, color="white", ha="center", va="center", fontsize=8.5, fontweight="bold", zorder=3)

    def pill(x, y, w, label, tr):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - 0.13), w, 0.26, boxstyle="round,pad=0,rounding_size=0.13",
                                    fc=TIER_COL[tr], ec="none", zorder=3))
        ax.text(x, y, label, ha="center", va="center", fontsize=7.6, fontweight="bold", zorder=4,
                color="white" if tr in ("High", "Very high", "Extreme") else PAL["ink"])

    def bar(x, y, w, frac, col):
        ax.add_patch(FancyBboxPatch((x, y - 0.055), w, 0.11, boxstyle="round,pad=0,rounding_size=0.055", fc=TRACK, ec="none", zorder=2))
        ax.add_patch(FancyBboxPatch((x, y - 0.055), max(w * frac, 0.06), 0.11, boxstyle="round,pad=0,rounding_size=0.055", fc=col, ec="none", zorder=3))

    for k, (_, r) in enumerate(T.iterrows()):
        yc = -(k + 0.5) * rh
        if k % 2 == 0:
            ax.add_patch(Rectangle((0, yc - rh / 2), W, rh, fc=ZEB, ec="none", zorder=0))
        ax.add_patch(Rectangle((0, yc - rh / 2), 0.09, rh, fc=TIER_COL[r[t_col]], ec="none", zorder=2))
        ax.text(cols[0][1] / 2, yc, str(k + 1), ha="center", va="center", fontsize=7.6, color=PAL["ink2"])
        ax.text(xs[1] + 0.1, yc, r["Location"], va="center", fontsize=8.2, fontweight="bold", color=PAL["ink"])
        ax.text(xs[2] + 0.1, yc, r["Watershed"], va="center", fontsize=7.8, color=PAL["ink2"])
        pill(xs[3] + cols[3][1] / 2, yc, 0.62, f"{r['PORI (0-100)']:.0f}", r["PORI tier"])
        grp = r[g_col]
        if grp == "none":
            ax.text(xs[4] + 0.15, yc, "\u2014", va="center", fontsize=8, color=PAL["ink2"])
        else:
            ax.add_patch(plt.Circle((xs[4] + 0.2, yc), 0.065, fc=GCOL.get(grp, PAL["magenta"]), ec="none", zorder=3))
            ax.text(xs[4] + 0.38, yc, grp, va="center", fontsize=8, color=PAL["ink"])
        bw = cols[6][1] - 0.85
        if mode == "class":
            cl = r[c_col]
            ax.text(xs[5] + 0.12, yc, "\u2014" if cl == "none" else cl, va="center", fontsize=8,
                    color=PAL["ink"] if cl != "none" else PAL["ink2"])
        else:
            v = float(r["Multidrug_share_of_ARG_load_%"])
            bar(xs[5] + 0.12, yc, cols[5][1] - 0.85, min(v, 100) / 100, TEAL)
            ax.text(xs[6] - 0.12, yc, f"{v:.0f}%", va="center", ha="right", fontsize=7.6, color=PAL["ink"])
        a = float(r[a_col])
        bar(xs[6] + 0.12, yc, bw, a / 100, TIER_COL[r[t_col]])
        ax.text(xs[7] - 0.12, yc, f"{a:.1f}", va="center", ha="right", fontsize=7.6, fontweight="bold", color=PAL["ink"])
        pill(xs[7] + cols[7][1] / 2, yc, 0.95, r[t_col], r[t_col])
    fig.text(0.03, 1 - 0.38 / fig_h, title, fontsize=13, fontweight="bold", color=PAL["ink"], va="center")
    fig.text(0.03, 1 - 0.72 / fig_h, subtitle, fontsize=8.5, color=PAL["ink2"], va="center")
    fig.text(0.03, 0.28 / fig_h, foot, fontsize=7.2, color=PAL["ink2"], va="center")
    fig.savefig(path, dpi=200, facecolor="white"); plt.close(fig)


def fig_ami_source_table(path, T, title, subtitle, foot):
    """Source-based AMI table: PORI, anthropogenic signal, where/how to act, avoidable PORI points (95% CI), AMI bar and tier."""
    from matplotlib.patches import FancyBboxPatch, Rectangle
    HDR, ZEB, TRACK = "#17525e", "#f6f4ef", "#e9e5dd"
    GCOL = {"Humans": PAL["blue"], "Cattle": PAL["orange"], "Pigs": PAL["aqua"], "Poultry": PAL["yellow"], "Finfish": PAL["magenta"]}
    T = T.sort_values("AMI (0-100)", ascending=False).reset_index(drop=True); n = len(T)
    cols = [("#", 0.5), ("Location", 1.7), ("Watershed", 1.9), ("PORI\nscore", 1.0), ("Anthropogenic\ntaxa (% reads)", 1.45),
            ("How to act\n(action focus)", 2.5), ("PORI points avoided\n(95% interval)", 2.1), ("AMI score", 2.0), ("Priority\nlevel", 1.2),
            ("Antibiotic-use context\n(top 2 groups, % total use)", 3.1)]
    xs = np.concatenate([[0], np.cumsum([w for _, w in cols])]); W = xs[-1]
    rh, hh, fw = 0.38, 0.95, 14.5
    unit = fw * 0.94 / W; ax_h = (n * rh + hh) * unit; fig_h = ax_h + 1.35 + 0.65
    fig = plt.figure(figsize=(fw, fig_h)); ax = fig.add_axes([0.03, 0.65 / fig_h, 0.94, ax_h / fig_h])
    ax.set_xlim(0, W); ax.set_ylim(-(n * rh), hh); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0, 0), W, hh, boxstyle="round,pad=0,rounding_size=0.18", fc=HDR, ec="none", zorder=1))
    for (h, w), x0 in zip(cols, xs[:-1]):
        ax.text(x0 + w / 2, hh / 2, h, color="white", ha="center", va="center", fontsize=8.3, fontweight="bold", zorder=3)
    def pill(x, y, w, label, fc, tc):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - 0.13), w, 0.26, boxstyle="round,pad=0,rounding_size=0.13", fc=fc, ec="none", zorder=3))
        ax.text(x, y, label, ha="center", va="center", fontsize=7.6, fontweight="bold", zorder=4, color=tc)
    for k, (_, r) in enumerate(T.iterrows()):
        yc = -(k + 0.5) * rh
        if k % 2 == 0: ax.add_patch(Rectangle((0, yc - rh / 2), W, rh, fc=ZEB, ec="none", zorder=0))
        ax.add_patch(Rectangle((0, yc - rh / 2), 0.09, rh, fc=TIER_COL[r["AMI tier"]], ec="none", zorder=2))
        ax.text(cols[0][1] / 2, yc, str(k + 1), ha="center", va="center", fontsize=7.6, color=PAL["ink2"])
        ax.text(xs[1] + 0.1, yc, r["Location"], va="center", fontsize=8.2, fontweight="bold", color=PAL["ink"])
        ax.text(xs[2] + 0.1, yc, r["Watershed"], va="center", fontsize=7.8, color=PAL["ink2"])
        pill(xs[3] + cols[3][1] / 2, yc, 0.62, f"{r['PORI (0-100)']:.0f}", TIER_COL[r["PORI tier"]],
             "white" if r["PORI tier"] in ("High", "Very high", "Extreme") else PAL["ink"])
        ax.text(xs[4] + cols[4][1] / 2, yc, f"{r['Anthropogenic_taxa_%']:.1f}", ha="center", va="center", fontsize=8, color=PAL["ink"])
        a_txt = r["Action focus"]; mon = a_txt == "Monitor"
        ax.text(xs[5] + 0.12, yc, a_txt, va="center", fontsize=8, color=PAL["ink2"] if mon else PAL["ink"], fontweight="normal" if mon else "bold")
        ax.text(xs[6] + cols[6][1] / 2, yc, f"{r['PORI_points_avoided']:.1f}  ({r['PORI_points_avoided_CI_low']:.1f}\u2013{r['PORI_points_avoided_CI_high']:.1f})",
                ha="center", va="center", fontsize=7.8, color=PAL["ink"])
        a = float(r["AMI (0-100)"]); bw = cols[7][1] - 0.85
        ax.add_patch(FancyBboxPatch((xs[7] + 0.12, yc - 0.055), bw, 0.11, boxstyle="round,pad=0,rounding_size=0.055", fc=TRACK, ec="none", zorder=2))
        ax.add_patch(FancyBboxPatch((xs[7] + 0.12, yc - 0.055), max(bw * a / 100, 0.06), 0.11, boxstyle="round,pad=0,rounding_size=0.055", fc=TIER_COL[r["AMI tier"]], ec="none", zorder=3))
        ax.text(xs[8] - 0.12, yc, f"{a:.1f}", va="center", ha="right", fontsize=7.6, fontweight="bold", color=PAL["ink"])
        pill(xs[8] + cols[8][1] / 2, yc, 0.95, r["AMI tier"], TIER_COL[r["AMI tier"]], "white" if r["AMI tier"] in ("High", "Very high", "Extreme") else PAL["ink"])
        ctx = r["Antibiotic-use context (surveillance)"]; g0 = ctx.split(" ")[0]
        ax.add_patch(plt.Circle((xs[9] + 0.16, yc), 0.06, fc=GCOL.get(g0, PAL["magenta"]), ec="none", zorder=3))
        ax.text(xs[9] + 0.32, yc, ctx, va="center", fontsize=7.6, color=PAL["ink2"])
    fig.text(0.03, 1 - 0.38 / fig_h, title, fontsize=13, fontweight="bold", color=PAL["ink"], va="center")
    fig.text(0.03, 1 - 0.72 / fig_h, subtitle, fontsize=8.5, color=PAL["ink2"], va="center")
    fig.text(0.03, 0.33 / fig_h, foot, fontsize=7.2, color=PAL["ink2"], va="center")
    fig.savefig(path, dpi=200, facecolor="white"); plt.close(fig)


def fig_pori_table(path, T, subtitle, foot):
    """
    Colour-coded PORI table: one row per station sorted by PORI, the five components (0-100, coloured by the
    same tier scale as PORI) grouped as resistome / pathogens / mobility, then PORI with an inline bar and tier chip.
    """
    from matplotlib.patches import FancyBboxPatch, Rectangle
    HDR, ZEB, TRACK = "#17525e", "#f6f4ef", "#e9e5dd"
    T = T.sort_values("PORI (0-100)", ascending=False).reset_index(drop=True)
    n = len(T)
    comps = [("ARGs_Total (0-100)", "Total\nARGs"), ("Critical_ARGs (0-100)", "Critical\nARGs"),
             ("Pathogens_Total (0-100)", "Total\npathogens"), ("Alert_Pathogens (0-100)", "High-alert\npathogens"),
             ("Mobility (0-100)", "Genetic\nmobility")]
    cw = 1.2
    cols = [("#", 0.55), ("Location", 1.9), ("Watershed", 2.0)] + [(h, cw) for _, h in comps] + [("PORI score", 2.3), ("PORI\ntier", 1.35)]
    xs = np.concatenate([[0], np.cumsum([w for _, w in cols])])
    W = xs[-1]
    rh, hh, gh = 0.38, 0.95, 0.42                      # row, header and group-band heights
    fw = 12.2
    unit = fw * 0.94 / W
    ax_h = (n * rh + hh + gh) * unit
    fig_h = ax_h + 1.35 + 0.55
    fig = plt.figure(figsize=(fw, fig_h))
    ax = fig.add_axes([0.03, 0.55 / fig_h, 0.94, ax_h / fig_h])
    ax.set_xlim(0, W); ax.set_ylim(-(n * rh), hh + gh); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0, 0), W, hh, boxstyle="round,pad=0,rounding_size=0.18", fc=HDR, ec="none", zorder=1))
    for (h, w), x0 in zip(cols, xs[:-1]):
        ax.text(x0 + w / 2, hh / 2, h, color="white", ha="center", va="center", fontsize=8.5, fontweight="bold", zorder=3)
    # group bands above the component columns
    for (a, b, lab, col, tc) in ((3, 5, "Resistome", PAL["blue"], "white"), (5, 7, "Pathogens", PAL["aqua"], "white"),
                                 (7, 8, "Mobility", PAL["yellow"], PAL["ink"])):
        ax.add_patch(FancyBboxPatch((xs[a] + 0.05, hh + 0.06), xs[b] - xs[a] - 0.1, gh - 0.12,
                                    boxstyle="round,pad=0,rounding_size=0.12", fc=col, ec="none", zorder=2))
        ax.text((xs[a] + xs[b]) / 2, hh + gh / 2, lab, ha="center", va="center", fontsize=8, fontweight="bold", color=tc, zorder=3)

    def pill(x, y, w, label, tr, h=0.26, fs=7.6):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle=f"round,pad=0,rounding_size={h / 2}",
                                    fc=TIER_COL[tr], ec="none", zorder=3))
        ax.text(x, y, label, ha="center", va="center", fontsize=fs, fontweight="bold", zorder=4,
                color="white" if tr in ("High", "Very high", "Extreme") else PAL["ink"])

    def bar(x, y, w, frac, col):
        ax.add_patch(FancyBboxPatch((x, y - 0.055), w, 0.11, boxstyle="round,pad=0,rounding_size=0.055", fc=TRACK, ec="none", zorder=2))
        ax.add_patch(FancyBboxPatch((x, y - 0.055), max(w * frac, 0.06), 0.11, boxstyle="round,pad=0,rounding_size=0.055", fc=col, ec="none", zorder=3))

    for k, (_, r) in enumerate(T.iterrows()):
        yc = -(k + 0.5) * rh
        if k % 2 == 0:
            ax.add_patch(Rectangle((0, yc - rh / 2), W, rh, fc=ZEB, ec="none", zorder=0))
        ax.add_patch(Rectangle((0, yc - rh / 2), 0.09, rh, fc=TIER_COL[r["PORI tier"]], ec="none", zorder=2))
        ax.text(cols[0][1] / 2, yc, str(k + 1), ha="center", va="center", fontsize=7.6, color=PAL["ink2"])
        ax.text(xs[1] + 0.1, yc, r["Location"], va="center", fontsize=8.2, fontweight="bold", color=PAL["ink"])
        ax.text(xs[2] + 0.1, yc, r["Watershed"], va="center", fontsize=7.8, color=PAL["ink2"])
        for j, (c, _) in enumerate(comps):
            v = float(r[c])
            pill(xs[3 + j] + cw / 2, yc, 0.8, f"{v:.0f}", tier(np.array([v]))[0])
        s = float(r["PORI (0-100)"])
        bar(xs[8] + 0.12, yc, cols[8][1] - 0.85, s / 100, TIER_COL[r["PORI tier"]])
        ax.text(xs[9] - 0.12, yc, f"{s:.1f}", va="center", ha="right", fontsize=7.6, fontweight="bold", color=PAL["ink"])
        pill(xs[9] + cols[9][1] / 2, yc, 0.95, r["PORI tier"], r["PORI tier"])
    fig.text(0.03, 1 - 0.38 / fig_h, "PORI – Pathogen-and-resistome One-health Risk Index by station", fontsize=13, fontweight="bold", color=PAL["ink"], va="center")
    fig.text(0.03, 1 - 0.72 / fig_h, subtitle, fontsize=8.5, color=PAL["ink2"], va="center")
    fig.text(0.03, 0.28 / fig_h, foot, fontsize=7.2, color=PAL["ink2"], va="center")
    fig.savefig(path, dpi=200, facecolor="white"); plt.close(fig)


def fig_pori_scheme(path):
    """
    Conceptual scheme of the five PORI components: every component is read as LOW or HIGH, giving 2^5 = 32 profiles
    (binary tree, left to right); the leaf colour counts how many components are HIGH. Soft gradient fills.
    """
    from matplotlib.patches import FancyBboxPatch, Polygon, PathPatch
    from matplotlib.path import Path as MPath
    from matplotlib.colors import LinearSegmentedColormap, to_rgb

    def mix(c1, c2, t):
        a, b = np.array(to_rgb(c1)), np.array(to_rgb(c2))
        return tuple(a * (1 - t) + b * t)

    # soft palettes: (light, deep) pairs -> horizontal gradient inside each chip
    LOW_G, HIGH_G = ("#e3f0ef", "#9fcfcc"), ("#fbe0c8", "#f0a272")
    LOW_T, HIGH_T = "#27585e", "#8f3b0c"
    HEAD = [(("#d6e6fa", "#7fb0ec"), "#16457d"), (("#d6e6fa", "#7fb0ec"), "#16457d"),
            (("#d5f1e6", "#7fd0ae"), "#0f5c3e"), (("#d5f1e6", "#7fd0ae"), "#0f5c3e"),
            (("#fcecc4", "#f3c76a"), "#6b4700")]
    names = ["ARGs", "High-risk\nARGs", "Pathogens", "High-alert\npathogens", "Genetic\nmobility"]
    L, nrow, rh = 5, 32, 0.2
    cx = [0.0, 1.75, 3.5, 5.25, 7.0]
    cwid = 1.2
    W = 10.3
    top_y = 1.55
    fw = 11.0
    unit = fw * 0.94 / (W + 1.6)
    ax_h = (nrow * rh + top_y + 0.9) * unit
    fig_h = ax_h + 1.35 + 0.55
    fig = plt.figure(figsize=(fw, fig_h))
    ax = fig.add_axes([0.03, 0.55 / fig_h, 0.94, ax_h / fig_h])
    XL, YL = (-1.0, W + 0.6), (-(nrow * rh) - 0.9, top_y)
    ax.set_xlim(*XL); ax.set_ylim(*YL); ax.axis("off")

    def gchip(x, y, w, h, g, z=3, rad=None, shadow=True, direction="h", gloss=False):
        """rounded chip with a soft gradient fill and a faint drop shadow"""
        rad = h / 2 if rad is None else rad
        if shadow:
            ax.add_patch(FancyBboxPatch((x + 0.02, y - 0.035), w, h, boxstyle=f"round,pad=0,rounding_size={rad}",
                                        fc="#000000", ec="none", alpha=0.07, zorder=z - 0.5))
        p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={rad}", fc="none", ec="none", zorder=z)
        ax.add_patch(p)
        cm = LinearSegmentedColormap.from_list("g", [g[0], g[1]])
        arr = np.linspace(0, 1, 64).reshape(1, -1) if direction == "h" else np.linspace(0, 1, 64).reshape(-1, 1)
        im = ax.imshow(arr, cmap=cm, extent=[x, x + w, y, y + h], aspect="auto", zorder=z, origin="lower" if direction == "h" else "upper")
        im.set_clip_path(p)
        # thin glossy highlight on the upper edge
        if gloss:
            ax.plot([x + rad * 0.9, x + w - rad * 0.9], [y + h * 0.86] * 2, color="white", lw=0.6, alpha=0.5, zorder=z + 0.2, solid_capstyle="round")

    # block bands: faint tint behind each ARGs block, to show the two big branches
    for b, col in ((0, "#eef6f5"), (1, "#fdf1e6")):
        y0 = -(b * 16 * rh)
        ax.add_patch(FancyBboxPatch((cx[0] - 0.75, y0 - 16 * rh + 0.03), (7.0 + cwid / 2 + 1.0) - (cx[0] - 0.75) + 0.3 + 0.55, 16 * rh - 0.06,
                                    boxstyle="round,pad=0,rounding_size=0.2", fc=col, ec="none", zorder=0))
    # headers
    for j, (nm, (g, tc)) in enumerate(zip(names, HEAD)):
        gchip(cx[j] - 0.78, 0.35, 1.56, 0.8, g, rad=0.2, direction="h", gloss=True)
        ax.text(cx[j], 0.75, nm, ha="center", va="center", fontsize=8.8, fontweight="bold", color=tc, zorder=5)
        ax.add_patch(plt.Circle((cx[j], 1.34), 0.085, fc=g[1], ec="none", zorder=3))
        ax.text(cx[j], 1.34, str(j + 1), ha="center", va="center", fontsize=6.2, fontweight="bold", color=tc, zorder=4)

    def ycen(level, idx):
        span = 2 ** (L - level - 1)
        return -(idx * span + span / 2) * rh

    # connectors first (under chips), tinted by the child's state
    for lev in range(L - 1):
        for idx in range(2 ** (lev + 1)):
            y = ycen(lev, idx)
            for c in (0, 1):
                y2 = ycen(lev + 1, idx * 2 + c)
                x0, x1 = cx[lev] + cwid / 2, cx[lev + 1] - cwid / 2
                xm = (x0 + x1) / 2
                col = HIGH_G[1] if c else LOW_G[1]
                p = MPath([(x0, y), (xm, y), (xm, y2), (x1, y2)], [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4])
                ax.add_patch(PathPatch(p, fc="none", ec=col, lw=1.5, alpha=0.75, zorder=1, capstyle="round"))
    for lev in range(L):
        for idx in range(2 ** (lev + 1)):
            hi = idx % 2
            y = ycen(lev, idx)
            span = 2 ** (L - lev - 1)
            ch = min(rh * span * 0.82, 0.34) if lev < L - 1 else rh * 0.82
            gchip(cx[lev] - cwid / 2, y - ch / 2, cwid, ch, HIGH_G if hi else LOW_G)
            ax.text(cx[lev], y, "HIGH" if hi else "LOW", ha="center", va="center", fontsize=6.2 if lev == L - 1 else 7,
                    fontweight="bold", color=HIGH_T if hi else LOW_T, zorder=5)
    # leaf counters, soft gradient along the sequential ramp
    xl = cx[-1] + cwid / 2 + 0.1
    soft = [mix(c, "#ffffff", 0.18) for c in SEQ_RAMP]
    for leaf in range(nrow):
        nh = sum((leaf >> (L - 1 - b)) & 1 for b in range(L))
        y = -(leaf + 0.5) * rh
        g = (mix(SEQ_RAMP[nh], "#ffffff", 0.3), SEQ_RAMP[nh])
        gchip(xl, y - rh * 0.41, 0.55, rh * 0.82, g, shadow=False)
        ax.text(xl + 0.275, y, f"{nh}/5", ha="center", va="center", fontsize=5.8, fontweight="bold", zorder=5,
                color="white" if nh >= 4 else PAL["ink"])
    # risk triangle with soft gradient
    xt0, xt1 = xl + 0.95, xl + 1.75
    ytop, ybot = -0.2, -(nrow * rh) + 0.1
    tri = Polygon([((xt0 + xt1) / 2, ytop), (xt1, ybot), (xt0, ybot)], closed=True, fc="none", ec="none", zorder=2)
    ax.add_patch(tri)
    cm = LinearSegmentedColormap.from_list("risk", ["#fdeede", "#f6c5a1", "#ee9c63", "#d9732f", "#a8460f", "#6b2c0e"])
    im = ax.imshow(np.linspace(0, 1, 100).reshape(-1, 1), cmap=cm, extent=[xt0, xt1, ybot, ytop], aspect="auto", zorder=2)
    im.set_clip_path(tri)
    ax.text((xt0 + xt1) / 2, ytop + 0.4, "Hazard", ha="center", va="center", fontsize=10, fontweight="bold", color=PAL["ink"])
    ax.text(xl + 0.275, ytop + 0.27, "n. HIGH", ha="center", va="center", fontsize=7, color=PAL["ink2"])
    # formula strip
    yb = -(nrow * rh) - 0.45
    gchip(-1.0, yb - 0.3, W + 1.55, 0.6, ("#f7f3ea", "#ece6d8"), z=1, rad=0.15, shadow=False)
    ax.text((W - 0.45) / 2, yb, r"PORI = $\sqrt{ARG \times critical\ ARG}\;\times\;\sqrt{pathogens \times high\text{-}alert}\;\times\;(1 + mobility)$   rescaled 0–100",
            ha="center", va="center", fontsize=9, color=PAL["ink"], zorder=5)
    ax.set_xlim(*XL); ax.set_ylim(*YL)
    fig.text(0.03, 1 - 0.38 / fig_h, "The five PORI components and how they combine into a hazard score", fontsize=13, fontweight="bold", color=PAL["ink"], va="center")
    fig.text(0.03, 1 - 0.72 / fig_h, "Each component read as LOW or HIGH gives 32 profiles; hazard tends to grow with the number of HIGH components, and only profiles with HIGH resistome AND HIGH pathogens reach the top",
             fontsize=8.5, color=PAL["ink2"], va="center")
    fig.text(0.03, 0.28 / fig_h,
             "Conceptual scheme: the index itself uses continuous 0-1 values (log1p, min-max across stations) and a geometric structure, so one LOW pillar keeps the score low.",
             fontsize=7.2, color=PAL["ink2"], va="center")
    fig.savefig(path, dpi=200, facecolor="white"); plt.close(fig)


def fig_substrate(path, SUB):
    """Stations of each substrate as dots (coloured by tier, spread side by side), with the mean as a bar."""
    T, pts = SUB["table"], SUB["points"]
    order = list(T["Substrate"])
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for k, sub in enumerate(order):
        d = pts[pts.Substrate == sub].sort_values("PORI").reset_index(drop=True)
        xs = k + (np.arange(len(d)) - (len(d) - 1) / 2) * 0.075
        m = T.loc[T.Substrate == sub, "Mean_PORI"].values[0]
        ax.plot([k - 0.32, k + 0.32], [m, m], color=PAL["ink"], lw=2, zorder=2)
        ax.scatter(xs, d.PORI, s=70, c=[TIER_COL[t] for t in d.tier], edgecolor=PAL["ink"], linewidth=0.8, zorder=3)
        top = d.iloc[-1]
        ax.annotate(top.Location, (xs[-1], top.PORI), xytext=(6, 2), textcoords="offset points", fontsize=7, color=PAL["ink2"])
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{s}\n(n = {int(n)})" for s, n in zip(T.Substrate, T.n_stations)])
    ax.set_ylabel("PORI (0-100, relative to the 36 stations)"); ax.set_ylim(-3, 106)
    for v in TIER_EDGES:
        ax.axhline(v, color=PAL["grid"], lw=0.8, ls=":", zorder=1)
    ax.set_title("PORI by substrate: each dot is a station, the bar is the substrate mean")
    tier_legend(ax, "upper right")
    fig.tight_layout(); fig.savefig(path, dpi=200); plt.close(fig)


def fig_outlier_map(path, D, ZF):
    """
    Outlier map: stations at their coordinates (overlapping points nudged side by side),
    colour = number of variables flagged as unusually high, callout boxes list the flagged
    variables with robust Z and two-sided p. Boxes sit in two columns, leader lines point
    to the stations.
    """
    meta = D["meta"]
    lon = pd.to_numeric(D["phys"]["Long (E)"]).values.astype(float)
    lat = pd.to_numeric(D["phys"]["Lat (N)"]).values.astype(float)
    names = meta["Location"].values
    cnt = ZF["counts"].values
    Z = ZF["flags"]

    # nudge overlapping points apart (distance measured in a lon/lat space with lat stretched)
    xy = np.c_[lon, lat].copy()
    rng = np.random.default_rng(SEED)
    dmin = 0.13
    for _ in range(300):
        moved = False
        for i in range(len(xy)):
            for j in range(i + 1, len(xy)):
                dx, dy = (xy[j, 0] - xy[i, 0]), (xy[j, 1] - xy[i, 1]) * 1.0
                d = np.hypot(dx, dy)
                if d < dmin:
                    if d < 1e-9:
                        a = rng.uniform(0, 2 * np.pi); dx, dy, d = np.cos(a) * 1e-3, np.sin(a) * 1e-3, 1e-3
                    push = (dmin - d) / 2
                    xy[i] -= np.array([dx, dy]) / d * push
                    xy[j] += np.array([dx, dy]) / d * push
                    moved = True
        if not moved:
            break

    def colour(n):
        if n == 0:
            return "#ffffff"
        return next(c for lo, hi, c in OUT_BINS if lo <= n <= hi)

    fig = plt.figure(figsize=(17, 12.5))
    axm = fig.add_axes([0.31, 0.07, 0.38, 0.86])
    axm.set_xlim(lon.min() - 0.5, lon.max() + 0.5); axm.set_ylim(lat.min() - 0.4, lat.max() + 0.4)
    axm.set_xlabel("Longitude (E)"); axm.set_ylabel("Latitude (N)")
    axm.scatter(xy[:, 0], xy[:, 1], s=95, c=[colour(n) for n in cnt], edgecolor=PAL["ink"], linewidth=0.9, zorder=3)
    axm.set_title("Stations with unusually high values (robust z > 3.5): number of variables flagged", fontsize=10)
    axm.legend(handles=[Patch(facecolor="#ffffff", edgecolor=PAL["ink"], label="0")] +
               [Patch(facecolor=c, edgecolor=PAL["ink"], label=(str(lo) if lo == hi else f"{lo}-{hi if hi < 99 else 10}")) for lo, hi, c in OUT_BINS],
               title="Variables flagged", loc="lower left", fontsize=8, ncol=2)

    # callout boxes: two columns, ordered by latitude (north on top)
    flagged = [i for i in np.argsort(-lat) if cnt[i] > 0]
    trend = np.polyfit(lat, lon, 1)
    # geographic side (west of the coastal trend -> left), rebalanced so neither column overflows
    left, right, rows_l, rows_r = [], [], 0, 0
    tot_rows = sum(min(int(cnt[i]), 9) + 1 for i in flagged)
    for i in flagged:
        r = min(int(cnt[i]), 9) + 1
        want_left = lon[i] <= np.polyval(trend, lat[i]) + 0.2
        if want_left and rows_l + r > 0.56 * tot_rows:
            want_left = False
        elif (not want_left) and rows_r + r > 0.56 * tot_rows:
            want_left = True
        (left if want_left else right).append(i)
        if want_left: rows_l += r
        else: rows_r += r
    ov = fig.add_axes([0, 0, 1, 1]); ov.set_xlim(0, 1); ov.set_ylim(0, 1); ov.axis("off")
    MAXROWS = 9
    for side, idxs in (("L", left), ("R", right)):
        nrows = [min(int(cnt[i]), MAXROWS) + 1 + (1 if cnt[i] > MAXROWS else 0) for i in idxs]
        rh = min(0.0132, 0.90 / max(sum(nrows) + 0.9 * len(idxs), 1))
        gap, y = 0.9 * rh, 0.95
        x0, w = (0.012, 0.27) if side == "L" else (0.718, 0.27)
        for i, nr in zip(idxs, nrows):
            h = nr * rh + 0.2 * rh
            ov.add_patch(plt.Rectangle((x0, y - h), w, h, facecolor="#dbe8e8", edgecolor="#b6cccc", lw=0.8))
            sub = Z[Z.Location == names[i]].head(MAXROWS)
            ov.text(x0 + 0.006, y - 0.7 * rh, names[i], fontsize=7.2, fontweight="bold", va="center")
            ov.text(x0 + w * 0.60, y - 0.7 * rh, "Z score", fontsize=7, fontweight="bold", va="center", ha="right")
            ov.text(x0 + w - 0.006, y - 0.7 * rh, "p value", fontsize=7, fontweight="bold", va="center", ha="right")
            for k, (_, r) in enumerate(sub.iterrows()):
                yy = y - (k + 1.7) * rh
                ov.text(x0 + 0.006, yy, str(r.Variable).replace(" (nr)", ""), fontsize=6.6, va="center")
                ov.text(x0 + w * 0.60, yy, f"{r.Z:.2f}", fontsize=6.6, va="center", ha="right")
                ov.text(x0 + w - 0.006, yy, f"{r.p:.4f}" if r.p >= 0.0001 else "<0.0001", fontsize=6.6, va="center", ha="right")
            if cnt[i] > MAXROWS:
                ov.text(x0 + 0.006, y - (MAXROWS + 1.7) * rh, f"+{int(cnt[i]) - MAXROWS} more", fontsize=6.4, va="center", ha="left", color=PAL["ink2"])
            px, py = axm.transData.transform(xy[i]); px, py = fig.transFigure.inverted().transform((px, py))
            ov.plot([x0 + (w if side == "L" else 0), px], [y - h / 2, py], color="#8fb0b0", lw=0.8, zorder=0)
            y -= h + gap * 0.35
    # station names for unflagged stations next to the dots
    for i in range(len(names)):
        if cnt[i] == 0:
            axm.annotate(names[i], xy[i], xytext=(5, -3), textcoords="offset points", fontsize=6, color=PAL["ink2"])
    fig.savefig(path, dpi=200); plt.close(fig)


# ======================================================================
# MAIN
# ======================================================================
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    np.random.seed(SEED)
    D = load_data()
    C = build_outcomes(D)
    AB = build_ab_structure(D["ab"])

    S1 = step1_correlations(D, C, AB)
    G1b = step1b_gene_level(D, C, AB)
    S1["gene_gate"] = G1b["gate_n"]
    S2 = step2_outliers(D, C, AB)
    ZF = step2_zflags(D, C, AB)
    P = step3_pori(D, C, S1)
    SP = step3b_spatial(D, C, P)
    SUB = step3c_substrate(D, P)
    S4 = step4_sentinels(D, C, P)
    S6 = step6_ami(D, C, AB, S1, P, S4)
    S6b = step6b_assignment_sensitivity(D, C, AB)
    S6c = step6c_multidrug_ami(D, C, AB, S6)
    S6s = step6s_ami_source(D, C, AB, P, S4)
    S6v = step6v_ami_validation(D, C, P, S4)

    # Sensitivity of PORI and AMI to the definition of "critical" ARGs (curated <-> keyword)
    alt_mode = "keyword" if C["mode"] == "curated" else "curated"
    QUIET[0] = True
    C_alt = build_outcomes(D, alt_mode)
    P_alt = step3_pori(D, C_alt, S1)
    S6_alt = step6_ami(D, C_alt, AB, S1, P_alt, S4)
    QUIET[0] = False
    m = P["table"].merge(P_alt["table"], on="Location", suffixes=("", "_alt"))
    a = S6["table"].merge(S6_alt["table"], on="Location", suffixes=("", "_alt"))
    SENS = a[["Location", "Target_group", "Target_class", "AMI (0-100)", "AMI tier",
              "Target_group_alt", "Target_class_alt", "AMI (0-100)_alt", "AMI tier_alt"]].rename(
        columns={c: c.replace("_alt", f" [{alt_mode}]") for c in a.columns if c.endswith("_alt")})
    same_t = ((a.Target_group == a.Target_group_alt) & (a.Target_class == a.Target_class_alt)).mean()
    log(f"[sensitivity] critical-ARG definition {C['mode']} vs {alt_mode}: PORI rho = "
        f"{spearmanr(m['PORI (0-100)'], m['PORI (0-100)_alt'])[0]:.2f}, same PORI tier = {(m['PORI tier'] == m['PORI tier_alt']).mean() * 100:.0f}%; "
        f"AMI rho = {spearmanr(a['AMI (0-100)'], a['AMI (0-100)_alt'])[0]:.2f}; same AMI tier = {(a['AMI tier'] == a['AMI tier_alt']).mean() * 100:.0f}%; "
        f"same target (group and class) = {same_t * 100:.0f}%; same group = {(a.Target_group == a.Target_group_alt).mean() * 100:.0f}%")

    readme = [
        "Adriatic One Health AMR workflow - results workbook",
        "PORI and AMI are descriptive, relative (0-100 over 36 stations) indices; tiers: "
        + ", ".join(f"{l} ({a}-{b})" for l, a, b in zip(TIER_LABELS, [0] + TIER_EDGES, TIER_EDGES + [100])),
        "S1_*: class-matched and group-level correlations (basin level n=17 primary; BH-FDR within family).",
        "S1b_*: gene-level consumption-resistome test (every gene vs every class, watershed level; matched vs cross pairs; permutation checks).",
        "S2_*: variable distribution summary and outlier flags per station.",
        "S3_*: PORI per station (components, rank, tier, sensitivity), validation against ARPA indicators, structure checks.",
        "S4_*: origin fractions (W/E/U), source signals, sentinel-taxa screen (environmental taxa; and E+U sensitivity).",
        "S6s_*: SOURCE-BASED AMI (main): PORI points avoided if the anthropogenic-taxa signal falls to the level of the cleaner stations (max cut 50%), with elasticities estimated from the data (S6s_elasticity); where/how to act per station.",
        "S6v_*: out-of-sample check of the source-based AMI: leave-one-watershed-out prediction of the components, within-watershed slopes and a pseudo-intervention test on station pairs of the same watershed (it does not test a real intervention).",
        "S6_*: previous AMI per station. AMI_total / AMI_critical: best (host group, class) target to lower the TOTAL / the CRITICAL resistome (absolute reduction rescaled 0-100; % of the station's own load; evidence labels). AMI_PORI (columns 'AMI (0-100)', dPORI): composite scenario on the whole PORI.",
        f"Critical ARG definition: {C['mode']} (alternative re-run reported in S6_sens_critical_def). Seed: {SEED}.",
    ]
    sheets = {
        "S1_class_matched": S1["A"], "S1_group_outcomes": S1["B"], "S1_heatmap_population_vs_use": S1["H"], "S1_consumption_vs_ARG_matrix": S1["M"],
        "S1b_gene_summary": G1b["SUM"], "S1b_gene_by_class": G1b["CL"], "S1b_gene_prevalence_sens": G1b["SENS"],
        "S1b_gene_pairs_fam1": G1b["T1"], "S1b_gene_pairs_critical": G1b["T2"], "S1b_gene_pairs_q_lt_0_10": G1b["sig"],
        "S2_variables": S2["summary"], "S2_stations": S2["stations"], "S2_zscore_flags": ZF["flags"],
        "S3_PORI": P["table"], "S3_spatial_structure": SP, "S3_PORI_by_substrate": SUB["table"], "S3_validation": P["validation"], "S3_sensitivity_structure": P["kv"],
        "S3_component_corr": P["corr"].reset_index().rename(columns={"index": "Component"}),
        "S4_origin_fractions": S4["origin"], "S4_origin_summary": S4["origin_summary"],
        "S4_sentinels_E": S4["sent"], "S4_sentinels_E+U": S4["sent_EU"], "S4_sentinels_W_nonpath": S4["sent_W"],
        "S6_AMI_stations": S6["table"], "S6_AMI_top3": S6["top3"],
        "S6_by_outcome_pct_reduction": S6["by_outcome"],
        "S6s_AMI_source": S6s["table"], "S6s_elasticity": S6s["elasticity"],
        "S6v_AMI_LOWO": S6v["lowo"], "S6v_AMI_within_slopes": S6v["within"], "S6v_AMI_pairs": S6v["pairs"], "S6v_AMI_pairs_summary": S6v["summary"],
        "S6c_multidrug_AMI": S6c["table"], "S6_sens_assign_summary": S6b["summary"], "S6_sens_assign_stations": S6b["stations"], "S6_by_group": S6["by_group"], "S6_by_class": S6["by_class"], "S6_sens_critical_def": SENS,
    }
    write_excel(os.path.join(OUT_DIR, "AMR_results.xlsx"), sheets, readme)
    make_figures(os.path.join(OUT_DIR, "figures"), D, C, S1, S2, ZF, P, S4, S6)
    fig_outlier_map(os.path.join(OUT_DIR, "figures", "Fig2_outlier_map.png"), D, ZF)
    fig_substrate(os.path.join(OUT_DIR, "figures", "Fig3b_PORI_by_substrate.png"), SUB)
    fig_gene_level(os.path.join(OUT_DIR, "figures", "Fig1b_gene_level_test.png"), G1b)
    fig_pori_table(os.path.join(OUT_DIR, "figures", "Fig3c_PORI_table.png"), P["table"],
                   "Five components rescaled 0-100 across the 36 stations, grouped as resistome, pathogens and mobility; chips use the PORI tier colours",
                   "PORI = geometric combination of the resistome pillars and the pathogen pillars, amplified by mobility, rescaled 0-100 (a relative hazard score, not a measure of exposure or of the probability of infection).\n"
                   "Tiers: Low <20, Moderate 20-40, High 40-60, Very high 60-80, Extreme >=80. Component values are relative to this set of stations.")
    fig_pori_scheme(os.path.join(OUT_DIR, "figures", "Fig0_PORI_scheme.png"))
    note = ("AMI = best mitigation target per station, rescaled 0-100 across the 36 stations (relative priority, not an effect size). "
            "Scenario: 50% cut in the use of a class by a host group, beta = 1 (upper bound), single-class genes only.\n"
            "Evidence for a consumption-resistome link (class and gene level) is not supported: targets are hypotheses to test. PORI chip colours = PORI tiers.")
    for tag, ttl in (("total", "Where to act to lower the total resistome"), ("critical", "Where to act to lower the critical resistome")):
        Tt = S6["table"].sort_values(f"AMI_{tag} (0-100)", ascending=False).reset_index(drop=True)
        fig_ami_table(os.path.join(OUT_DIR, "figures", f"Fig7_AMI_{tag}_table.png"), Tt, f"AMI \u2013 {ttl}",
                      "Mitigation group and antibiotic class whose reduced use is expected to lower the resistome most, station by station", "class", note,
                      g_col=f"Target_group_{tag}", a_col=f"AMI_{tag} (0-100)", t_col=f"AMI_{tag} tier", c_col=f"Target_class_{tag}")
    fig_ami_source_table(os.path.join(OUT_DIR, "figures", "Fig7_AMI_table.png"), S6s["table"],
                         "AMI \u2013 where and how to act on the anthropogenic source",
                         "PORI points avoided if the signal of anthropogenic taxa is brought down to the level of the cleaner stations (maximum cut 50%), using elasticities estimated from the data",
                         ("AMI = PORI points avoided, rescaled 0-100 across the 36 stations (relative priority). Elasticities: log-log OLS with salinity, basin-cluster bootstrap; only components whose 95% CI excludes zero respond (ARGs, pathogens; not mobility).\n"
                          "Interval = bootstrap 2.5-97.5%. WWTP = wastewater treatment plant. Action focus is indicative (composition of the signal: waste-associated vs gut-associated taxa). Cross-sectional relationship: the gain holds only if the signal causally drives the resistome (out-of-sample check in Table 1).\n"
                          "Antibiotic-use context = the two groups (humans included) with the largest share of total antibiotic use in the basin; it is not part of the AMI (use and livestock abundance show no association with the resistome)."))
    fig_ami_table(os.path.join(OUT_DIR, "figures", "FigS3_AMI_multidrug_table.png"), S6c["table"],
                  "Supplementary AMI \u2013 multidrug resistance genes",
                  "Multidrug genes respond to several classes, so the target is the host group whose antibiotic use to reduce, not a class", "multi",
                  ("AMI = best host group per station, rescaled 0-100 across the 36 stations (relative priority, not an effect size). Scenario: 50% cut in the use of the group's antibiotics, beta = 1 (upper bound).\n"
                   "Multi-class genes only; each gene is assumed to respond equally to every class that affects it (mean share of the group over its classes).\n"
                   "Evidence for a consumption-resistome link is not supported: the target is a hypothesis to test. PORI chip colours = PORI tiers."),
                  g_col="Target_group_multidrug", a_col="AMI_multidrug (0-100)", t_col="AMI_multidrug tier")
    with open(os.path.join(OUT_DIR, "run_log.txt"), "w") as fh:
        fh.write("\n".join(LOG))
    log(f"Done. Results in '{OUT_DIR}'.")


if __name__ == "__main__":
    main()
