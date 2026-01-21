#!/usr/bin/env python3
import argparse
import os
import sys
from io import StringIO
from typing import Dict, List, Tuple, Optional
import pandas as pd
import yaml

CORE_PREFIXES = ["stx", "wzx", "wzy", "wzt", "wzm", "flic", "fli", "fl", "eae", "ehxa"]

FINAL_DETAIL_SPECS: List[Tuple[str, str]] = [
    ("toxin_details", "toxin"),
    ("O_type_details", "Otype"),
    ("H_type_details", "Htype"),
    ("type_details", "adhesin_virulence"),
    ("other_details", "other"),
]

# ------------------------- Parsing & Threshold resolution ------------------------- #

def parse_gene_from_template(template: str) -> Tuple[str, str]:
    """
    Extract (gene, allele) from KMA '#Template'.
    Supports '__' as primary delimiter, '_' as fallback.
    """
    NA = "-"
    parts = template.split("__") if "__" in template else template.split("_")
    gene = parts[1] if len(parts) >= 2 else parts[0]
    allele = parts[2] if len(parts) >= 3 else NA
    return gene, allele

def resolve_threshold_for_gene(gene: str, thresholds: Dict[str, List[float]]) -> List[float]:
    """
    Priority:
      1) exact key match
      2) case-insensitive exact
      3) prefix match (e.g., 'stx1' uses 'stx' thresholds)
      4) 'other' (case-insensitive)
    Returns a list: [coverage_min, identity_min, depth_min]
    """
    if gene in thresholds:
        return thresholds[gene]
    for k in thresholds:
        if k.lower() == gene.lower():
            return thresholds[k]
    for k in thresholds:
        if k.lower() != "other" and gene.lower().startswith(k.lower()):
            return thresholds[k]
    for k in thresholds:
        if k.lower() == "other":
            return thresholds[k]
    # if we get here, metafile lacks both a matching entry and 'other'
    raise ValueError(f"No threshold for gene '{gene}' and no 'other' fallback in metafile.")

def resolve_threshold_key_for_gene(gene: str, thresholds: Dict[str, List[float]]) -> str:
    """
    Same matching logic as resolve_threshold_for_gene(), but returns the *matched key*.
    This lets us filter out hits that only match via the 'other' fallback.
    """
    if gene in thresholds:
        return gene
    for k in thresholds:
        if k.lower() == gene.lower():
            return k
    for k in thresholds:
        if k.lower() != "other" and gene.lower().startswith(k.lower()):
            return k
    for k in thresholds:
        if k.lower() == "other":
            return k
    raise ValueError(f"No threshold for gene '{gene}' and no 'other' fallback in config.")

# ------------------------- Step 1: Read YAML gene config ------------------------- #

def read_geneconfig(config_path: str, organism_key: str) -> Dict[str, List[float]]:
    """
    Reads YAML config like:
      Escherichia:
        stx1: [98,98,1]
        ...
      Shigella:
        ipaH: [98,98,1]
        ...

    Returns: dict { gene_name: [template_cov_min, query_ident_min, depth_min] }
    """
    if yaml is None:
        raise RuntimeError("Missing dependency 'pyyaml'. Install with: pip install pyyaml")

    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as fh:
        cfg = yaml.safe_load(fh)

    if not isinstance(cfg, dict):
        raise ValueError("Config YAML must be a mapping at the top level (dict).")

    if organism_key not in cfg:
        raise ValueError(
            f"Organism key '{organism_key}' not found in config. Available: {list(cfg.keys())}"
        )

    section = cfg[organism_key]
    if not isinstance(section, dict):
        raise ValueError(f"Config section for '{organism_key}' must be a mapping of gene -> [cov,id,depth].")

    thresholds: Dict[str, List[float]] = {}
    for gene, vals in section.items():
        if not isinstance(gene, str):
            raise ValueError(f"Gene name must be a string. Got: {gene!r}")
        if not (isinstance(vals, list) and len(vals) == 3):
            raise ValueError(f"Threshold for gene '{gene}' must be a 3-item list: [cov,id,depth]. Got: {vals!r}")
        thresholds[gene] = [float(vals[0]), float(vals[1]), float(vals[2])]

    return thresholds

# ------------------------- Step 2: Read + Extract KMA .res ------------------------- #

def process_kma_res(res_path: str, thresholds: Dict[str, List[float]]) -> pd.DataFrame:
    """
    Read KMA .res, parse gene/allele from #Template, and KEEP ONLY rows where gene matches
    something in the config (excluding pure 'other' fallback).

    Returns DF with columns:
      gene, allele, Template_Coverage, Query_Identity, Depth
    """
    if not os.path.exists(res_path):
        raise FileNotFoundError(f"KMA .res file not found: {res_path}")

    with open(res_path, "r") as fh:
        lines = fh.readlines()

    header_idx: Optional[int] = None
    for i, line in enumerate(lines):
        if line.strip().startswith("#Template"):
            header_idx = i
            break

    if header_idx is None:
        raise ValueError("Could not find header line starting with '#Template' in the .res file.")

    content = "".join(lines[header_idx:])
    df = pd.read_csv(StringIO(content), sep=r"\s+", engine="python")

    required_cols = {"#Template", "Template_Coverage", "Query_Identity", "Depth"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in .res: {sorted(missing)}. Found: {df.columns.tolist()}")

    # Parse gene + allele
    parsed = df["#Template"].apply(parse_gene_from_template)
    df["gene"] = parsed.apply(lambda x: x[0])
    df["allele"] = parsed.apply(lambda x: x[1])

    """
    # Determine which config key this gene would match (or 'other')
    df["__matched_key"] = df["gene"].apply(lambda g: resolve_threshold_key_for_gene(g, thresholds))

    # Keep only real matches, NOT entries that only map to 'other'
    df = df[df["__matched_key"].str.lower() != "other"].copy()
    """

    # Determine which config key this gene would match (or 'other')
    df["__matched_key"] = df["gene"].apply(lambda g: resolve_threshold_key_for_gene(g, thresholds))
    # Return exactly the requested columns/order

    out = df[["gene", "allele", "Template_Coverage", "Query_Identity", "Depth"]].copy()
    return out

# ------------------------- Placeholder functions (empty for now) ------------------------- #
def threshold_triplet_for_key(matched_key: str, thresholds: Dict[str, List[float]]) -> Tuple[float, float, float]:
    """
    Given a key that exists in thresholds (e.g. 'wzx', 'stx1', or 'other'),
    return (cov_min, id_min, depth_min) as floats.
    """
    t = thresholds[matched_key]
    return float(t[0]), float(t[1]), float(t[2])

def filter_kma_res(df: pd.DataFrame, thresholds: Dict[str, List[float]]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    For each row:
      - match gene to a config key (exact / case-insensitive / prefix)
      - if no match -> treat as 'other'
      - apply [Template_Coverage, Query_Identity, Depth] thresholds

    Returns:
      pass_df, fail_df
    """
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"filter_kma_res: input df missing columns: {sorted(missing)}")

    work = df.copy()

    # 1) resolve gene -> matched threshold key (or 'other')
    matched_keys: List[str] = [
        resolve_threshold_key_for_gene(str(g), thresholds)
        for g in work["gene"].tolist()
    ]
    work["matched_key"] = matched_keys

    # 2) attach thresholds per row (based on matched_key, so unmatched -> 'other')
    triplets: List[Tuple[float, float, float]] = [
        threshold_triplet_for_key(k, thresholds)
        for k in matched_keys
    ]
    work["th_cov"] = [t[0] for t in triplets]
    work["th_id"] = [t[1] for t in triplets]
    work["th_depth"] = [t[2] for t in triplets]

    # 3) apply thresholds
    cov = work["Template_Coverage"].astype(float)
    ide = work["Query_Identity"].astype(float)
    dep = work["Depth"].astype(float)

    pass_mask = (cov >= work["th_cov"]) & (ide >= work["th_id"]) & (dep >= work["th_depth"])
    work["status"] = ["PASS" if x else "FAIL" for x in pass_mask.tolist()]

    out_cols = [
        "gene", "allele", "matched_key",
        "Template_Coverage", "Query_Identity", "Depth",
        "th_cov", "th_id", "th_depth", "status"
    ]
    pass_df = work.loc[pass_mask, out_cols].copy()
    fail_df = work.loc[~pass_mask, out_cols].copy()

    return pass_df, fail_df


# ------------------------- Placeholder functions (empty for now) ------------------------- #

def collect_alleles(df: pd.DataFrame, prefix: str) -> List[str]:
    """
    Collect ALL alleles (including duplicates) for rows where gene startswith(prefix),
    excluding '-' and empty.
    """
    genes_lower = df["gene"].astype(str).str.lower()
    alleles = df["allele"].astype(str)

    mask = genes_lower.str.startswith(prefix.lower())
    out: List[str] = []
    for a in alleles[mask].tolist():
        if a and a != "-":
            out.append(a)
    return out

def build_details_from_df(df: pd.DataFrame, prefixes: List[str], include_gene: bool = False) -> str:
    """
    Build details string for rows whose gene startswith ANY prefix in prefixes.

    Format per row:
      - include_gene=False: allele_cov_id_depth
      - include_gene=True:  gene_allele_cov_id_depth

    Joined by ';' in df order. If none -> '-'.
    """
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"build_details_from_df: df missing columns: {sorted(missing)}")

    if df.empty:
        return "-"

    genes = df["gene"].astype(str)
    genes_lower = genes.str.lower()

    mask = pd.Series([False] * len(df), index=df.index)
    for p in prefixes:
        mask = mask | genes_lower.str.startswith(p.lower())

    if not mask.any():
        return "-"

    subset = df.loc[mask, ["gene", "allele", "Template_Coverage", "Query_Identity", "Depth"]]
    parts: List[str] = []
    for row in subset.itertuples(index=False):
        gene = str(row.gene)
        allele = str(row.allele)
        cov = float(row.Template_Coverage)
        qid = float(row.Query_Identity)
        dep = float(row.Depth)

        if include_gene:
            parts.append(f"{gene}_{allele}_{cov:.2f}_{qid:.2f}_{dep:.2f}")
        else:
            parts.append(f"{allele}_{cov:.2f}_{qid:.2f}_{dep:.2f}")

    return "-" if len(parts) == 0 else ";".join(parts)

#---------------------------

def determine_stx_subtype(pass_df: pd.DataFrame, fail_df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    """
    Output columns:
      sample_id, stx1, stx2, Toxin, toxin_pass_details, toxin_fail_details, toxin_details

    stx1/stx2: ALL alleles from PASS lines (duplicates kept).
    Toxin: concordant PASS subtype calls only (1 unique allele per subtype).
    toxin_pass_details: per-line allele+metrics for PASS stx* lines.
    toxin_fail_details: per-line allele+metrics for FAIL stx* lines.
    toxin_details: 'pass:{toxin_pass_details}_fail:{toxin_fail_details}'
    """
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing_pass = required - set(pass_df.columns)
    missing_fail = required - set(fail_df.columns)

    if missing_pass:
        raise ValueError(f"determine_stx_subtype: pass_df missing columns: {sorted(missing_pass)}")
    if missing_fail:
        raise ValueError(f"determine_stx_subtype: fail_df missing columns: {sorted(missing_fail)}")

    # PASS allele lists (keep duplicates)
    stx1_alleles = collect_alleles(pass_df, "stx1")
    stx2_alleles = collect_alleles(pass_df, "stx2")

    stx1_col = "-" if len(stx1_alleles) == 0 else ";".join(stx1_alleles)
    stx2_col = "-" if len(stx2_alleles) == 0 else ";".join(stx2_alleles)

    # Concordance check (unique among PASS)
    stx1_unique = sorted(set(stx1_alleles))
    stx2_unique = sorted(set(stx2_alleles))

    toxin_parts: List[str] = []
    if len(stx1_unique) == 1:
        toxin_parts.append(stx1_unique[0])
    if len(stx2_unique) == 1:
        toxin_parts.append(stx2_unique[0])

    toxin_call = "-" if len(toxin_parts) == 0 else ";".join(toxin_parts)

    toxin_pass_details = build_details_from_df(pass_df, ["stx"], include_gene=False)
    toxin_fail_details = build_details_from_df(fail_df, ["stx"], include_gene=False)

    toxin_details = f"pass:{toxin_pass_details}|fail:{toxin_fail_details}"

    stx_out = pd.DataFrame([{
        "sample_id": sample_id,
        "stx1": stx1_col,
        "stx2": stx2_col,
        "Toxin": toxin_call,
        "toxin_pass_details": toxin_pass_details,
        "toxin_fail_details": toxin_fail_details,
        "toxin_details": toxin_details,
    }])

    return stx_out

def determine_O_type(pass_df: pd.DataFrame, fail_df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    """
    Output columns:
      sample_id | wzx | wzy | wzt | wzm | O_type | O_type_pass_details | O_type_fail_details | O_type_details

    Allele columns (wzx/wzy/wzt/wzm): ALL PASS alleles for that gene (duplicates kept).
    O_type is called ONLY if:
      1) wzx and wzy are concordant (single unique allele each, and equal), OR
      2) wzm and wzt are concordant (single unique allele each, and equal), OR
      3) all four are concordant to the same allele (single unique allele for each and all equal)

    If none of these hold (or conflicting pair calls), O_type = '-'.
    """
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing_pass = required - set(pass_df.columns)
    missing_fail = required - set(fail_df.columns)
    if missing_pass:
        raise ValueError(f"determine_O_type: pass_df missing columns: {sorted(missing_pass)}")
    if missing_fail:
        raise ValueError(f"determine_O_type: fail_df missing columns: {sorted(missing_fail)}")

    # Collect ALL PASS alleles per gene (keep duplicates)
    wzx_alleles = collect_alleles(pass_df, "wzx")
    wzy_alleles = collect_alleles(pass_df, "wzy")
    wzt_alleles = collect_alleles(pass_df, "wzt")
    wzm_alleles = collect_alleles(pass_df, "wzm")

    wzx_col = "-" if len(wzx_alleles) == 0 else ";".join(wzx_alleles)
    wzy_col = "-" if len(wzy_alleles) == 0 else ";".join(wzy_alleles)
    wzt_col = "-" if len(wzt_alleles) == 0 else ";".join(wzt_alleles)
    wzm_col = "-" if len(wzm_alleles) == 0 else ";".join(wzm_alleles)

    # Unique alleles (for concordance logic)
    wzx_unique = sorted(set(wzx_alleles))
    wzy_unique = sorted(set(wzy_alleles))
    wzt_unique = sorted(set(wzt_alleles))
    wzm_unique = sorted(set(wzm_alleles))

    candidate1 = None  # wzx/wzy
    candidate2 = None  # wzt/wzm
    candidate3 = None  # all four

    # (3) all four concordant (strongest) - 1 unique in each
    if len(wzx_unique) == 1 and len(wzy_unique) == 1 and len(wzt_unique) == 1 and len(wzm_unique) == 1:
        union = set([wzx_unique[0], wzy_unique[0], wzt_unique[0], wzm_unique[0]])
        if len(union) == 1: #simply pick one as they are all identical
            candidate3 = wzx_unique[0]

    # (1) wzx and wzy concordant -> if the wzx for all four hasn't been filled out yet
    if candidate3 is None:
        if len(wzx_unique) == 1 and len(wzy_unique) == 1 and wzx_unique[0] == wzy_unique[0]: 
            candidate1 = wzx_unique[0]

        # (2) wzt and wzm concordant
        if len(wzt_unique) == 1 and len(wzm_unique) == 1 and wzt_unique[0] == wzm_unique[0]:
            candidate2 = wzt_unique[0]

    # Decide O_type
    if candidate3 is not None:
        O_type = candidate3
    else:
        # If both candidates exist but differ -> no call
        if candidate1 is not None and candidate2 is not None and candidate1 != candidate2:
            O_type = "-"
        elif candidate1 is not None: # if only one candidate differ 
            O_type = candidate1
        elif candidate2 is not None:
            O_type = candidate2
        else:
            O_type = "-"

    prefixes = ["wzx", "wzy", "wzt", "wzm"]
    O_type_pass_details = build_details_from_df(pass_df, prefixes, include_gene=True)
    O_type_fail_details = build_details_from_df(fail_df, prefixes, include_gene=True)
    O_type_details = f"pass:{O_type_pass_details}|fail:{O_type_fail_details}"

    out = pd.DataFrame([{
        "sample_id": sample_id,
        "wzx": wzx_col,
        "wzy": wzy_col,
        "wzt": wzt_col,
        "wzm": wzm_col,
        "O_type": O_type,
        "O_type_pass_details": O_type_pass_details,
        "O_type_fail_details": O_type_fail_details,
        "O_type_details": O_type_details,
    }])

    return out

def determine_H_type(pass_df: pd.DataFrame, fail_df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    """
    Determine H type using:
      - Prefer fliC if present (PASS) and concordant
      - Else use other 'fl*' genes with allele like H\\d+ (PASS) and concordant

    Output columns:
      sample_id | flagellin | H_type | H_type_pass_details | H_type_fail_details | H_type_details
    """
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing_pass = required - set(pass_df.columns)
    missing_fail = required - set(fail_df.columns)
    if missing_pass:
        raise ValueError(f"determine_H_type: pass_df missing columns: {sorted(missing_pass)}")
    if missing_fail:
        raise ValueError(f"determine_H_type: fail_df missing columns: {sorted(missing_fail)}")

    # 1) Prefer fliC
    fliC_alleles = collect_alleles(pass_df, "fliC")

    use_prefixes: List[str]
    fl_alleles: List[str]

    if len(fliC_alleles) > 0:
        # Use fliC only (even if discordant, we do NOT fall back)
        use_prefixes = ["fliC"]
        fl_alleles = fliC_alleles
    else:
        # 2) Fall back to other fl* genes (includes fli*, flnA, flkA, etc.)
        use_prefixes = ["fl"]
        fl_alleles = collect_alleles(pass_df, "fl")

    # flagellin column: list ALL alleles seen in the chosen source (duplicates kept)
    flagellin = "-" if len(fl_alleles) == 0 else ";".join(fl_alleles)

    # H_type: only if concordant (exactly 1 unique allele)
    unique = sorted(set(fl_alleles))
    if len(unique) == 1:
        H_type = unique[0]          # e.g. "H51"
    else:
        H_type = "-"                # none or discordant

    # Details (include gene name)
    H_type_pass_details = build_details_from_df(pass_df, use_prefixes, include_gene=True)
    H_type_fail_details = build_details_from_df(fail_df, use_prefixes, include_gene=True)
    H_type_details = f"pass:{H_type_pass_details}|fail:{H_type_fail_details}"

    out = pd.DataFrame([{
        "sample_id": sample_id,
        "flagellin": flagellin,
        "H_type": H_type,
        "H_type_pass_details": H_type_pass_details,
        "H_type_fail_details": H_type_fail_details,
        "H_type_details": H_type_details,
    }])

    return out

def determine_adhesin_virulence(pass_df: pd.DataFrame, fail_df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    """
    Output columns:
      sample_id | Adhesin | Virulence | type_pass_details | type_fail_details | type_details

    Adhesin:  'positive' if eae present in PASS else '-'
    Virulence:'positive' if ehxA present in PASS else '-'

    Details include gene name + allele + cov + id + depth for eae/ehxA rows.
    """
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing_pass = required - set(pass_df.columns)
    missing_fail = required - set(fail_df.columns)
    if missing_pass:
        raise ValueError(f"determine_adhesin_virulence: pass_df missing columns: {sorted(missing_pass)}")
    if missing_fail:
        raise ValueError(f"determine_adhesin_virulence: fail_df missing columns: {sorted(missing_fail)}")

    eae_positive = collect_alleles(pass_df, "eae")
    ehxA_positive = collect_alleles(pass_df, "ehxA")

    Adhesin = "positive" if eae_positive else "-"
    Virulence = "positive" if ehxA_positive else "-"

    prefixes = ["eae", "ehxA"]
    type_pass_details = build_details_from_df(pass_df, prefixes, include_gene=True)
    type_fail_details = build_details_from_df(fail_df, prefixes, include_gene=True)
    type_details = f"pass:{type_pass_details}|fail:{type_fail_details}"

    out = pd.DataFrame([{
        "sample_id": sample_id,
        "Adhesin": Adhesin,
        "Virulence": Virulence,
        "type_pass_details": type_pass_details,
        "type_fail_details": type_fail_details,
        "type_details": type_details,
    }])

    return out

def is_missing_detail_value(v: object) -> bool:
    if v is None:
        return True
    if isinstance(v, float) and pd.isna(v):
        return True
    s = str(v).strip()
    return s == "" or s == "-"

def build_verbose_from_detail_columns(row: pd.Series) -> str:
    parts: List[str] = []
    for col, label in FINAL_DETAIL_SPECS:
        if col in row and not is_missing_detail_value(row[col]):
            parts.append(f"{label}:{row[col]}")
    return "-" if len(parts) == 0 else "||".join(parts)

def merge_final_outputs(
    stx_out: pd.DataFrame,
    o_out: pd.DataFrame,
    h_out: pd.DataFrame,
    av_out: pd.DataFrame,
    other_out: pd.DataFrame,
) -> pd.DataFrame:
    """
    Merge per-module outputs on sample_id and collapse the *_details (NOT pass/fail)
    into a single 'verbose' column.

    Final output columns ONLY:
      sample_id | Toxin | O_type | H_type | Adhesin | Virulence | other | verbose
    """
    final_df = stx_out.merge(o_out, on="sample_id", how="outer")
    final_df = final_df.merge(h_out, on="sample_id", how="outer")
    final_df = final_df.merge(av_out, on="sample_id", how="outer")
    final_df = final_df.merge(other_out, on="sample_id", how="outer") 

    # build combined verbose column from the final *_details fields
    final_df["verbose"] = final_df.apply(build_verbose_from_detail_columns, axis=1)

    # drop pass/fail detail columns
    drop_cols: List[str] = []
    for c in final_df.columns:
        if c.endswith("_pass_details") or c.endswith("_fail_details"):
            drop_cols.append(c)

    # drop the individual final *_details columns (they're now in 'verbose')
    for col, _label in FINAL_DETAIL_SPECS:
        if col in final_df.columns:
            drop_cols.append(col)

    drop_cols = sorted(set(drop_cols))
    if drop_cols:
        final_df = final_df.drop(columns=drop_cols)

    # Keep ONLY these columns in final output
    keep = ["sample_id", "Toxin", "O_type", "H_type", "Adhesin", "Virulence", "other", "verbose"]
    for c in keep:
        if c not in final_df.columns:
            final_df[c] = "-"

    final_df = final_df[keep].copy()

    # replace NaN/None with '-'
    final_df = final_df.fillna("-")

    return final_df

def build_details_excluding_prefixes(df: pd.DataFrame, excluded_prefixes: List[str]) -> str:
    """
    Details for rows whose gene does NOT start with any excluded prefix.
    Format per row: gene_allele_cov_id_depth
    Joined by ';' in df order. If none -> '-'
    """
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"build_details_excluding_prefixes: df missing columns: {sorted(missing)}")

    if df.empty:
        return "-"

    genes_lower = df["gene"].astype(str).str.lower()

    mask_excl = pd.Series([False] * len(df), index=df.index)
    for p in excluded_prefixes:
        mask_excl = mask_excl | genes_lower.str.startswith(p.lower())

    subset = df.loc[~mask_excl, ["gene", "allele", "Template_Coverage", "Query_Identity", "Depth"]]
    if subset.empty:
        return "-"

    parts: List[str] = []
    for row in subset.itertuples(index=False):
        gene = str(row.gene)
        allele = str(row.allele)
        cov = float(row.Template_Coverage)
        qid = float(row.Query_Identity)
        dep = float(row.Depth)
        parts.append(f"{gene}_{allele}_{cov:.2f}_{qid:.2f}_{dep:.2f}")

    return "-" if len(parts) == 0 else ";".join(parts)

def determine_other(pass_df: pd.DataFrame, fail_df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    """
    Output:
      sample_id | Other | other_pass_details | other_fail_details | other_details

    Other: unique gene names among PASS rows that are NOT part of the core typing prefixes.
    Details: gene+allele+metrics for PASS/FAIL rows in this 'other' group.
    """
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing_pass = required - set(pass_df.columns)
    missing_fail = required - set(fail_df.columns)
    if missing_pass:
        raise ValueError(f"determine_other: pass_df missing columns: {sorted(missing_pass)}")
    if missing_fail:
        raise ValueError(f"determine_other: fail_df missing columns: {sorted(missing_fail)}")

    # Which PASS genes belong to "other" group?
    genes_lower = pass_df["gene"].astype(str).str.lower()
    mask_excl = pd.Series([False] * len(pass_df), index=pass_df.index)
    for p in CORE_PREFIXES:
        mask_excl = mask_excl | genes_lower.str.startswith(p)

    other_pass_genes = sorted(set(pass_df.loc[~mask_excl, "gene"].astype(str).tolist()))
    Other = "-" if len(other_pass_genes) == 0 else ";".join(other_pass_genes)

    other_pass_details = build_details_excluding_prefixes(pass_df, CORE_PREFIXES)
    other_fail_details = build_details_excluding_prefixes(fail_df, CORE_PREFIXES)
    other_details = f"pass:{other_pass_details}|fail:{other_fail_details}"

    return pd.DataFrame([{
        "sample_id": sample_id,
        "Other": Other,
        "other_pass_details": other_pass_details,
        "other_fail_details": other_fail_details,
        "other_details": other_details,
    }])

# ------------------------- Main ------------------------- #

def normalize_organism_key(user_organism: str) -> str:
    org = user_organism.strip()

    ecoli_aliases = {"Escherichia", "Escherichia coli", "E. coli", "E.coli"}
    shigella_aliases = {
        "Shigella",
        "Shigella sonnei", "S. sonnei", "S.sonnei",
        "Shigella flexneri", "S. flexneri", "S.flexneri",
        "Shigella boydii", "S. boydii", "S.boydii",
        "Shigella dysenteriae", "S. dysenteriae", "S.dysenteriae",
    }

    if org in ecoli_aliases:
        return "Escherichia"
    if org in shigella_aliases:
        return "Shigella"

    raise ValueError(f"Unrecognized organism '{user_organism}'.")


def main(args: argparse.Namespace) -> None:
    try:
        organism_key = normalize_organism_key(args.organism)
        thresholds = read_geneconfig(args.config, organism_key)

        extracted_df = process_kma_res(args.KMA_res,thresholds)
        pass_df, fail_df = filter_kma_res(extracted_df, thresholds)

        # --output is a PREFIX
        prefix = args.output
        out_dir = os.path.dirname(prefix)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

        pass_path = f"{prefix}_pass.tsv"
        fail_path = f"{prefix}_fail.tsv"

        pass_df.to_csv(pass_path, sep="\t", index=False)
        fail_df.to_csv(fail_path, sep="\t", index=False)

        print(f"Organism input: {args.organism!r} -> config key: {organism_key!r}", file=sys.stderr)
        print(f"Loaded {len(thresholds)} thresholds from {args.config}", file=sys.stderr)
        print(f"Read {len(extracted_df)} hits from {args.KMA_res}", file=sys.stderr)
        print(f"PASS: {len(pass_df)} -> {pass_path}", file=sys.stderr)
        print(f"FAIL: {len(fail_df)} -> {fail_path}", file=sys.stderr)

        # ----- store allele information in gene specific files -------
        stx_out = determine_stx_subtype(pass_df,fail_df,args.sample_id)
        o_out = determine_O_type(pass_df, fail_df, args.sample_id)
        h_out = determine_H_type(pass_df, fail_df, args.sample_id)
        av_out = determine_adhesin_virulence(pass_df, fail_df, args.sample_id)
        other_out = determine_other(pass_df, fail_df, args.sample_id)

        if args.store_tmp:
            # Store d
            stx_path = f"{prefix}_stx.tsv"
            stx_out.to_csv(stx_path, sep="\t", index=False)
            print(f"STX: wrote {stx_path}", file=sys.stderr)

            o_path = f"{prefix}_O.tsv"
            o_out.to_csv(o_path, sep="\t", index=False)
            print(f"O: wrote {o_path}", file=sys.stderr)

            h_path = f"{prefix}_H.tsv"
            h_out.to_csv(h_path, sep="\t", index=False)
            print(f"H: wrote {h_path}", file=sys.stderr)
        
            av_path = f"{prefix}_adhesin_virulence.tsv"
            av_out.to_csv(av_path, sep="\t", index=False)
            print(f"Adhesin/Virulence: wrote {av_path}", file=sys.stderr)

            other_path = f"{prefix}_other.tsv"
            other_out.to_csv(other_path, sep="\t", index=False)
            print(f"Other: wrote {other_path}", file=sys.stderr)
        else:
            print("store_tmp not set: skipping *_stx.tsv, *_O.tsv, *_H.tsv, *_adhesin_virulence.tsv", file=sys.stderr)

        if args.verbose:
            final_df = merge_final_outputs(stx_out, o_out, h_out, av_out, other_out)
            final_path = f"{prefix}_final.tsv"
            final_df.to_csv(final_path, sep="\t", index=False)
            print(f"FINAL: wrote {final_path}", file=sys.stderr)
        else:
            print("verbose not set: skipping *_final.tsv", file=sys.stderr)

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Read YAML gene thresholds and filter a KMA .res into PASS/FAIL files. '--output' is a prefix."
    )
    parser.add_argument("--KMA_res", required=True, help="Input KMA .res file")
    parser.add_argument("--config", required=True, help="YAML config with per-organism gene thresholds")
    parser.add_argument("--organism", required=True, help="Organism string (e.g., 'Escherichia coli', 'Shigella sonnei')")
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--store_tmp",action="store_true",help="If set, write extra summary TSVs: _stx.tsv, _O.tsv, _H.tsv, _adhesin_virulence.tsv")
    parser.add_argument("--verbose",action="store_true",help="If set, write {prefix}_final.tsv by merging summary outputs and collapsing *_details into one 'verbose' column.")
    parser.add_argument("--output", required=True, help="Output prefix (writes {prefix}_pass.tsv and {prefix}_fail.tsv)")
    args = parser.parse_args()
    main(args)

#python ecolityping_filter.py --KMA_res Result/SRR5023683_test.res --config GeneFilter.yaml --organism Escherichia --sample_id SRR5023683 --output FilterTest