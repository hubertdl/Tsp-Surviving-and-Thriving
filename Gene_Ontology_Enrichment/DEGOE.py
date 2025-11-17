#!/usr/bin/env python3
"""
DEGOE: Differential Expression -> Gene Ontology Enrichment
v1.0 (2025)
Written by E. J. Bentz

DEGOE takes differential-expression results tables from DESeq or generic scores 
from either single files or whole folders with multiple files to process in batch.
It performs Gene Ontology enrichment using flaxible and tunable methods, and creates 
summaries and usable output in one step, without the need to manually create score files.
It is designed for quick, repeatable runs. Most jobs finish in under a minute.

Key capabilities:
 - Load DESeq2 results and extract scores automatically
 - Batch-process every score file in a given directory
 - Expand GO annotations with ancestor terms from the go.obo definitions file
 - Apply optional GO-level and category-size filters
 - Optionally collapse redundant GO terms by shared genes and keep only the
     most informative representative per cluster.
 - Run either Mann–Whitney U (rank-based) or Fisher’s exact tests to find GO
     categories enriched and apply Benjamini–Hochberg FDR 
 - Log parameters used for every run

Tests available:
  - Mann–Whitney U (MWU) enrichment (Non-thresholded test)
  - Fisher's exact test enrichment (Thresholded input)

(Future versions are planned to include methods such as GSEA and ORA)

"""

import os
import sys
import argparse
import logging
from pathlib import Path
from datetime import datetime
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster

# Limit threads for BLAS/OpenMP libraries before importing anything heavy
def _get_cpu_arg_for_env():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--cpus', type=int, default=1)
    args, _ = parser.parse_known_args()
    return args.cpus

_n_cpus = _get_cpu_arg_for_env()
os.environ['OPENBLAS_NUM_THREADS'] = str(_n_cpus)
os.environ['MKL_NUM_THREADS'] = str(_n_cpus)
os.environ['OMP_NUM_THREADS'] = str(_n_cpus)
os.environ['BLIS_NUM_THREADS'] = str(_n_cpus)
os.environ['NUMEXPR_NUM_THREADS'] = str(_n_cpus)
os.environ['VECLIB_MAXIMUM_THREADS'] = str(_n_cpus)

# Logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


############################
# GO DATABASE AND ANNOTATIONS
############################

def ensure_go_obo(go_file: str = "go.obo") -> str:
    """
    Ensure that a GO OBO file exists and is reasonably recent.
    Downloads the current version if missing or older than 30 days.
    """
    go_path = Path(go_file)
    if go_path.exists():
        age_days = (datetime.now().timestamp() - go_path.stat().st_mtime) / (24 * 3600)
        if age_days < 30:
            logger.info(f"Using existing GO OBO file '{go_file}' (age ~{age_days:.1f} days)")
            return str(go_path)
        else:
            logger.info(f"Existing GO OBO file '{go_file}' is old (~{age_days:.1f} days), refreshing...")
    else:
        logger.info(f"GO OBO file '{go_file}' not found, downloading...")

    try:
        import urllib.request
        url = "http://current.geneontology.org/ontology/go.obo"
        urllib.request.urlretrieve(url, go_file)
        logger.info("Successfully downloaded go.obo")
    except Exception as e:
        logger.error(f"Failed to download go.obo: {e}")
        if not go_path.exists():
            logger.error("No local go.obo file available. Exiting.")
            sys.exit(1)
        logger.warning("Using existing go.obo despite download failure")

    return str(go_path)


class GOHierarchy:
    """
    Parse GO OBO file and provide:
      - names: GO ID -> term name
      - parents: GO ID -> set of parent GO IDs
      - levels: GO ID -> hierarchy level (root = 0)
    Only one namespace (BP/MF/CC) is kept.
    """

    def __init__(self, obo_file: str, namespace_code: str = "BP"):
        code_to_name = {
            "BP": "biological_process",
            "MF": "molecular_function",
            "CC": "cellular_component"
        }
        if namespace_code not in code_to_name:
            raise ValueError(f"Unknown namespace code '{namespace_code}'. Use BP, MF, or CC.")
        self.namespace_code = namespace_code
        self.namespace = code_to_name[namespace_code]

        self.names = {}           # GO ID -> name
        self.parents = defaultdict(set)  # GO ID -> parents
        self.levels = {}          # GO ID -> level
        self._ancestor_cache = {} # GO ID -> set of ancestors

        self._parse_obo(obo_file)
        self._compute_levels()

    def _parse_obo(self, obo_file: str):
        logger.info(f"Parsing GO OBO file '{obo_file}' for namespace '{self.namespace}'")
        current_id = None
        current_name = None
        current_ns = None

        with open(obo_file, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue

                if line == "[Term]":
                    current_id = None
                    current_name = None
                    current_ns = None
                    continue

                if line.startswith("id: GO:"):
                    current_id = line.split()[1]
                    continue

                if line.startswith("name:"):
                    current_name = line[6:]
                    continue

                if line.startswith("namespace:"):
                    current_ns = line.split()[1]
                    continue

                if line.startswith("is_a:") and current_id is not None:
                    parent_id = line.split()[1]
                    # Only store parent relationships within the chosen namespace
                    if current_ns == self.namespace:
                        self.parents[current_id].add(parent_id)

                # At the end of a term block, if namespace matches, record the term
                # This is done by checking when we see the next "[Term]" or EOF.
                # The simple implementation above depends on the last values seen.

        # Second pass to collect names and filter to namespace
        # Some OBOs do not repeat namespace information near parents, so we re-parse for names.
        current_id = None
        current_name = None
        current_ns = None

        with open(obo_file, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue

                if line == "[Term]":
                    if current_id and current_ns == self.namespace and current_name:
                        self.names[current_id] = current_name
                    current_id = None
                    current_name = None
                    current_ns = None
                    continue

                if line.startswith("id: GO:"):
                    current_id = line.split()[1]
                    continue

                if line.startswith("name:"):
                    current_name = line[6:]
                    continue

                if line.startswith("namespace:"):
                    current_ns = line.split()[1]
                    continue

            # Last term
            if current_id and current_ns == self.namespace and current_name:
                self.names[current_id] = current_name

        logger.info(f"Loaded {len(self.names)} GO terms in namespace '{self.namespace}'")

        # Filter parents to keep only terms in this namespace
        filtered_parents = defaultdict(set)
        for term, parents in self.parents.items():
            if term in self.names:
                filtered_parents[term] = {p for p in parents if p in self.names}
        self.parents = filtered_parents

    def _compute_levels(self):
        """
        Compute approximate hierarchy levels using BFS from roots.
        A root is any term with no parents in this namespace.
        """
        roots = [t for t in self.names if not self.parents.get(t)]
        level = {}
        from collections import deque
        q = deque()

        for r in roots:
            q.append((r, 0))

        while q:
            term, lvl = q.popleft()
            if term in level:
                continue
            level[term] = lvl
            # Children are any terms that list this term as a parent
            for child, parents in self.parents.items():
                if term in parents and child not in level:
                    q.append((child, lvl + 1))

        self.levels = level
        logger.info(f"Computed levels for {len(self.levels)} GO terms")

    def ancestors(self, term: str):
        """
        Return all ancestor terms (parents, grandparents, etc.) for a GO term.
        Cached for speed.
        """
        if term in self._ancestor_cache:
            return self._ancestor_cache[term]

        seen = set()
        stack = [term]
        while stack:
            t = stack.pop()
            for p in self.parents.get(t, []):
                if p not in seen:
                    seen.add(p)
                    stack.append(p)
        self._ancestor_cache[term] = seen
        return seen


def load_annotations(annotation_file: str, go_hierarchy: GOHierarchy) -> pd.DataFrame:
    """
    Load gene-to-GO annotations from a tab-delimited file:
      gene_id  GO:0000001;GO:0000002;...
    Expand each term to include its ancestors in the chosen namespace.
    Returns a DataFrame with columns:
      gene, go_term, go_name, level
    """
    logger.info(f"Loading annotations from '{annotation_file}'")
    raw = pd.read_csv(
        annotation_file,
        sep="\t",
        header=None,
        names=["gene", "go_terms"],
        dtype=str
    )

    records = []
    gene_to_expanded_terms = {}

    for _, row in raw.iterrows():
        gene = row["gene"]
        terms_str = row["go_terms"]
        if pd.isna(terms_str) or not terms_str or terms_str == "unknown":
            continue

        base_terms = {t.strip() for t in terms_str.split(";") if t.strip()}
        expanded_terms = set()

        for t in base_terms:
            if t in go_hierarchy.names:
                expanded_terms.add(t)
                expanded_terms.update(go_hierarchy.ancestors(t))

        expanded_terms = {t for t in expanded_terms if t in go_hierarchy.names}
        if not expanded_terms:
            continue

        gene_to_expanded_terms[gene] = sorted(expanded_terms)

        for t in expanded_terms:
            records.append({
                "gene": gene,
                "go_term": t,
                "go_name": go_hierarchy.names.get(t, t),
                "level": go_hierarchy.levels.get(t, -1)
            })

    if not records:
        logger.error("No valid gene-GO annotations found after expansion")
        sys.exit(1)

    df = pd.DataFrame.from_records(records)
    logger.info(
        f"Expanded annotations: {df['gene'].nunique()} genes, "
        f"{df['go_term'].nunique()} GO terms, {len(df)} gene-term pairs"
    )
    return df, gene_to_expanded_terms


def save_expanded_annotations(
    gene_to_terms: dict,
    output_dir: Path,
    original_annotation_file: str
) -> None:
    """
    Save expanded annotations in a simple format:
      gene_id  GO:...;GO:...
    The filename includes the original annotation filename stem.
    """
    input_stem = Path(original_annotation_file).stem
    out_file = output_dir / f"GO_Ancestors_{input_stem}.tsv"
    logger.info(f"Saving expanded annotations to '{out_file}'")

    with open(out_file, "w") as fh:
        for gene, terms in gene_to_terms.items():
            fh.write(f"{gene}\t{';'.join(terms)}\n")


############################
# SCORES LOADING
############################

class ScoreFileProcessor:
    """
    Utilities for loading gene-level scores from DESeq2 results or generic tables.
    """

    @staticmethod
    def _detect_format(filepath: str):
        """
        Heuristics:
          - Separator: first of [tab, comma, space] that splits the first line.
          - has_header: assume True.
          - is_deseq: header has N columns, second line has N+1 (rownames).
        """
        with open(filepath, "r") as fh:
            first = fh.readline().strip()
            second = fh.readline().strip()

        if not first:
            raise ValueError(f"Score file '{filepath}' is empty")

        sep = "\t"
        for candidate in ["\t", ",", " "]:
            if candidate in first:
                sep = candidate
                break

        first_cols = first.split(sep)
        second_cols = second.split(sep) if second else []

        has_header = True
        is_deseq = False
        if second_cols and len(second_cols) == len(first_cols) + 1:
            is_deseq = True

        return {
            "sep": sep,
            "has_header": has_header,
            "is_deseq": is_deseq,
            "n_cols_header": len(first_cols),
            "n_cols_second": len(second_cols)
        }

    @staticmethod
    def load_scores(
        filepath: str,
        score_type: str,
        score_column: str = None,
        pvalue_column: str = "pvalue",
        fold_change_column: str = "log2FoldChange",
        respect_independent_filtering: bool = False
    ) -> pd.DataFrame:
        """
        Load gene scores and return a DataFrame with columns:
          gene, value

        score_type:
          - "raw": take raw values from score_column (for Fisher, raw p)
          - "pvalue": value = -log10(pvalue)
          - "wald": value = Wald statistic
          - "log2foldchange": value = log2FoldChange
          - "signed_pvalue": value = -log10(p) * sign(log2FoldChange)
          - "pval_signedlog2fc": value = -log10(p) * log2FoldChange
        """
        info = ScoreFileProcessor._detect_format(filepath)
        sep = info["sep"]

        if info["is_deseq"]:
            df = pd.read_csv(filepath, sep=sep, header=0, index_col=0)
            df = df.reset_index()
            df.rename(columns={df.columns[0]: "gene"}, inplace=True)
            gene_col = "gene"
        else:
            df = pd.read_csv(filepath, sep=sep, header=0)
            gene_col = df.columns[0]
            df.rename(columns={gene_col: "gene"}, inplace=True)
            gene_col = "gene"

        if respect_independent_filtering:
            initial = len(df)
            cols_to_check = []
            for col in df.columns:
                cl = col.lower()
                if cl == "pvalue" or "pvalue" in cl:
                    cols_to_check.append(col)
                if cl == "padj" or "padj" in cl:
                    cols_to_check.append(col)
            if cols_to_check:
                mask = np.ones(len(df), dtype=bool)
                for col in cols_to_check:
                    mask &= df[col].notna().values
                df = df.loc[mask].copy()
                removed = initial - len(df)
                logger.info(
                    f"Independent filtering: removed {removed}/{initial} genes "
                    f"with NA in {cols_to_check}"
                )
            else:
                logger.warning(
                    "respect_independent_filtering requested, "
                    "but pvalue/padj columns were not found"
                )

        # Now derive a single numeric "value" column
        df_values = None

        if score_type == "raw":
            col = score_column or pvalue_column
            if col not in df.columns:
                raise ValueError(f"Score column '{col}' not found in score file")
            df_values = df[["gene", col]].rename(columns={col: "value"})

        elif score_type == "pvalue":
            if pvalue_column not in df.columns:
                raise ValueError(f"P-value column '{pvalue_column}' not found")
            pvals = df[pvalue_column].astype(float).replace(0, 1e-300)
            values = -np.log10(pvals)
            df_values = pd.DataFrame({"gene": df["gene"], "value": values})
            logger.info(
                f"score_type='pvalue': "
                f"-log10(p) range [{values.min():.2f}, {values.max():.2f}]"
            )

        elif score_type == "wald":
            if fold_change_column not in df.columns:
                raise ValueError(f"Wald statistic column '{fold_change_column}' not found")
            vals = df[fold_change_column].astype(float)
            df_values = pd.DataFrame({"gene": df["gene"], "value": vals})
            logger.info(
                f"score_type='wald': range [{vals.min():.2f}, {vals.max():.2f}]"
            )

        elif score_type == "log2foldchange":
            if fold_change_column not in df.columns:
                raise ValueError(f"Fold change column '{fold_change_column}' not found")
            vals = df[fold_change_column].astype(float)
            df_values = pd.DataFrame({"gene": df["gene"], "value": vals})
            logger.info(
                f"score_type='log2foldchange': range [{vals.min():.2f}, {vals.max():.2f}]"
            )

        elif score_type == "signed_pvalue":
            if pvalue_column not in df.columns:
                raise ValueError(f"P-value column '{pvalue_column}' not found")
            if fold_change_column not in df.columns:
                raise ValueError(f"Fold change column '{fold_change_column}' not found")
            pvals = df[pvalue_column].astype(float).replace(0, 1e-300)
            l10 = -np.log10(pvals)
            fc = df[fold_change_column].astype(float)
            vals = l10 * np.sign(fc)
            df_values = pd.DataFrame({"gene": df["gene"], "value": vals})
            logger.info(
                f"score_type='signed_pvalue': range [{vals.min():.2f}, {vals.max():.2f}]"
            )

        elif score_type == "pval_signedlog2fc":
            if pvalue_column not in df.columns:
                raise ValueError(f"P-value column '{pvalue_column}' not found")
            if fold_change_column not in df.columns:
                raise ValueError(f"Fold change column '{fold_change_column}' not found")
            pvals = df[pvalue_column].astype(float).replace(0, 1e-300)
            l10 = -np.log10(pvals)
            fc = df[fold_change_column].astype(float)
            vals = l10 * fc
            df_values = pd.DataFrame({"gene": df["gene"], "value": vals})
            logger.info(
                f"score_type='pval_signedlog2fc': range [{vals.min():.2f}, {vals.max():.2f}]"
            )

        else:
            raise ValueError(f"Unsupported score_type '{score_type}'")

        df_values = df_values.dropna()
        logger.info(f"Loaded {len(df_values)} gene scores from '{filepath}'")
        return df_values


############################
# COMMON UTILITIES
############################

def filter_by_go_level(annotations: pd.DataFrame,
                       min_level: int | None,
                       max_level: int | None) -> pd.DataFrame:
    if min_level is None and max_level is None:
        return annotations

    before_terms = annotations["go_term"].nunique()
    before_rows = len(annotations)

    mask = np.ones(len(annotations), dtype=bool)
    if min_level is not None:
        mask &= annotations["level"].values >= min_level
    if max_level is not None:
        mask &= annotations["level"].values <= max_level

    filtered = annotations.loc[mask].copy()
    after_terms = filtered["go_term"].nunique()
    after_rows = len(filtered)

    logger.info(
        f"GO level filter: kept {after_terms}/{before_terms} terms, "
        f"{after_rows}/{before_rows} rows (min={min_level}, max={max_level})"
    )
    if after_terms == 0:
        logger.warning("All GO terms removed by level filter")
    return filtered


def filter_categories_by_size(merged: pd.DataFrame,
                              min_genes: int,
                              max_genes: int | None,
                              max_genes_fraction: float | None) -> pd.DataFrame:
    """
    Filter GO categories by number of unique genes per term.
    If max_genes_fraction is given, max_genes = fraction * total_genes.
    If both max_genes and max_genes_fraction are given, use the smaller.
    If neither is given, default max_genes = 10% of total genes.
    """
    total_genes = merged["gene"].nunique()
    if max_genes_fraction is not None:
        if not (0 < max_genes_fraction <= 1.0):
            raise ValueError("max_genes_fraction must be in (0, 1]")
        frac_limit = int(total_genes * max_genes_fraction)
        if max_genes is None:
            max_genes = frac_limit
        else:
            max_genes = min(max_genes, frac_limit)
    elif max_genes is None:
        max_genes = int(0.1 * total_genes)

    if max_genes <= 0:
        raise ValueError("max_genes must be positive")

    sizes = merged.groupby("go_term")["gene"].nunique()
    keep_terms = sizes[(sizes >= min_genes) & (sizes <= max_genes)].index

    filtered = merged[merged["go_term"].isin(keep_terms)].copy()
    logger.info(
        f"Category size filter: kept {len(keep_terms)} GO terms "
        f"(min={min_genes}, max={max_genes}, total genes={total_genes})"
    )
    return filtered


def collapse_redundant_terms(
    annotations: pd.DataFrame,
    cut_height: float
) -> pd.DataFrame:
    """Cluster GO categories based on shared genes and merge redundant ones."""

    if cut_height <= 0:
        logger.warning("cluster cut height must be positive; skipping redundancy reduction")
        return annotations

    terms = annotations["go_term"].unique()
    if len(terms) < 2:
        logger.info("Redundancy reduction skipped (fewer than 2 GO terms)")
        return annotations

    gene_sets = (
        annotations.groupby("go_term")["gene"].apply(lambda genes: set(genes))
    )
    gene_sets = gene_sets.to_dict()

    metadata = (
        annotations.drop_duplicates("go_term")
        .set_index("go_term")[["go_name", "level"]]
    )

    terms = list(gene_sets.keys())
    condensed = []
    for i in range(len(terms) - 1):
        genes_i = gene_sets[terms[i]]
        len_i = len(genes_i)
        for j in range(i + 1, len(terms)):
            genes_j = gene_sets[terms[j]]
            len_j = len(genes_j)
            smaller = min(len_i, len_j)
            if smaller == 0:
                dist = 1.0
            else:
                shared = len(genes_i & genes_j)
                dist = 1.0 - (shared / smaller)
            condensed.append(dist)

    if not condensed:
        logger.info("Redundancy reduction skipped (insufficient term pairs)")
        return annotations

    condensed = np.clip(condensed, 0.0, 1.0)

    try:
        link = linkage(condensed, method="complete")
        cluster_ids = fcluster(link, t=cut_height, criterion="distance")
    except Exception as exc:
        logger.warning(f"GO term clustering failed ({exc}); skipping redundancy reduction")
        return annotations

    clusters = defaultdict(list)
    for term, cluster_id in zip(terms, cluster_ids):
        clusters[cluster_id].append(term)

    def level_value(term: str) -> float:
        if term not in metadata.index:
            return -1
        lvl = metadata.at[term, "level"]
        return float(lvl) if pd.notna(lvl) else -1.0

    replacements = {}
    merge_details = []
    for cluster_terms in clusters.values():
        if len(cluster_terms) < 2:
            continue
        cluster_terms.sort(
            key=lambda term: (len(gene_sets[term]), level_value(term)),
            reverse=True
        )
        representative = cluster_terms[0]
        rep_name = metadata.at[representative, "go_name"] if representative in metadata.index else representative
        rep_size = len(gene_sets[representative])
        for term in cluster_terms[1:]:
            replacements[term] = representative
            term_name = metadata.at[term, "go_name"] if term in metadata.index else term
            merge_details.append(
                f"  - '{term_name}' ({term}) -> '{rep_name}' ({representative}) "
                f"[{len(gene_sets[term])}→{rep_size} genes]"
            )

    if not replacements:
        logger.info("Redundancy reduction found no clusters below cut height %.2f", cut_height)
        return annotations

    if len(merge_details) <= 10:
        for detail in merge_details:
            logger.info(detail)
    else:
        for detail in merge_details[:10]:
            logger.info(detail)
        logger.info("  ... %d additional merges", len(merge_details) - 10)

    collapsed = annotations.copy()
    collapsed["go_term"] = collapsed["go_term"].map(lambda term: replacements.get(term, term))

    name_lookup = metadata["go_name"].to_dict()
    level_lookup = metadata["level"].to_dict()
    collapsed["go_name"] = collapsed["go_term"].map(lambda term: name_lookup.get(term, term))
    collapsed["level"] = collapsed["go_term"].map(lambda term: level_lookup.get(term, -1))

    before_terms = annotations["go_term"].nunique()
    collapsed = collapsed.drop_duplicates(subset=["gene", "go_term"])
    after_terms = collapsed["go_term"].nunique()
    logger.info(
        f"Redundancy reduction collapsed {len(replacements)} GO terms: "
        f"{before_terms} -> {after_terms} unique categories at cut height {cut_height}"
    )

    return collapsed


def add_bh_fdr(results: pd.DataFrame) -> pd.DataFrame:
    """
    Add Benjamini–Hochberg FDR column 'p.adj' to a results DataFrame
    that already has a 'pval' column.
    """
    if results.empty:
        results["p.adj"] = []
        return results

    try:
        from statsmodels.stats.multitest import multipletests
    except ImportError:
        logger.warning("statsmodels not available; using raw p-values as 'p.adj'")
        results["p.adj"] = results["pval"].values
        return results

    pvals = results["pval"].values
    _, p_adj, _, _ = multipletests(pvals, method="fdr_bh")
    results["p.adj"] = p_adj
    return results


############################
# MANN–WHITNEY U ENRICHMENT
############################

def run_mwu(
    scores: pd.DataFrame,
    annotations: pd.DataFrame,
    score_type: str,
    test_direction: str,
    min_genes: int,
    max_genes: int | None,
    max_genes_fraction: float | None,
    source_label: str | None = None
) -> pd.DataFrame:
    """
    MWU enrichment on merged gene scores and GO annotations.

    scores: DataFrame with columns gene, value
    annotations: DataFrame with gene, go_term, go_name, level
    test_direction options:
      - 'enriched' (non-directional, for pvalue scores)
      - 'up', 'down', 'both' (for directional scores)
    """

    # Merge score and annotation tables
    merged = annotations.merge(scores, on="gene", how="inner")
    if merged.empty:
        raise ValueError("No overlapping genes between scores and annotations")

    merged = filter_categories_by_size(merged, min_genes, max_genes, max_genes_fraction)

    # Aggregate per gene to a single score (if duplicated)
    gene_values = merged.groupby("gene")["value"].mean()
    all_genes = gene_values.index.values
    all_values = gene_values.values

    # Precompute ranks for delta.rank
    all_ranks = stats.rankdata(all_values)
    gene_ranks = dict(zip(all_genes, all_ranks))

    categories = merged["go_term"].unique()
    label = source_label or getattr(scores, "source_file", "input")
    logger.info(f"MWU: testing {len(categories)} GO categories from '{label}'")

    def test_category(term: str,
                      alternative: str) -> dict | None:
        genes_in_term = set(merged.loc[merged["go_term"] == term, "gene"])
        genes_outside = set(all_genes) - genes_in_term
        if len(genes_in_term) < 2 or len(genes_outside) < 2:
            return None

        vals_in = [gene_values[g] for g in genes_in_term]
        vals_out = [gene_values[g] for g in genes_outside]

        stat, pval = stats.mannwhitneyu(vals_in, vals_out, alternative=alternative)

        ranks_in = [gene_ranks[g] for g in genes_in_term]
        ranks_out = [gene_ranks[g] for g in genes_outside]
        delta_rank = np.mean(ranks_in) - np.mean(ranks_out)

        row = merged.loc[merged["go_term"] == term].iloc[0]
        return {
            "term": term,
            "name": row["go_name"],
            "level": row["level"],
            "nseqs": len(genes_in_term),
            "delta.rank": delta_rank,
            "pval": pval
        }

    results_rows = []

    if score_type == "pvalue":
        # Non-directional: larger -log10(p) means more significant.
        # We only support 'enriched' here.
        if test_direction != "enriched":
            raise ValueError(
                "score_type='pvalue' must be used with test_direction='enriched'"
            )
        for term in categories:
            res = test_category(term, alternative="greater")
            if res is not None:
                results_rows.append(res)

    else:
        # Directional score types
        if test_direction not in {"up", "down", "both"}:
            raise ValueError(
                f"score_type='{score_type}' must be used with test_direction "
                "'up', 'down', or 'both'"
            )

        if test_direction in {"up", "both"}:
            for term in categories:
                res = test_category(term, alternative="greater")
                if res is not None:
                    res["direction"] = "up"
                    results_rows.append(res)

        if test_direction in {"down", "both"}:
            for term in categories:
                res = test_category(term, alternative="less")
                if res is not None:
                    res["direction"] = "down"
                    results_rows.append(res)

    if not results_rows:
        logger.warning("MWU: no valid categories tested")
        return pd.DataFrame(
            columns=["term", "name", "level", "nseqs", "delta.rank", "pval", "p.adj"]
        )

    results = pd.DataFrame.from_records(results_rows)
    results = add_bh_fdr(results)

    # Column order
    cols = ["term", "name", "level", "nseqs"]
    if "direction" in results.columns:
        cols.append("direction")
    cols.extend(["delta.rank", "pval", "p.adj"])
    results = results[cols].sort_values("p.adj").reset_index(drop=True)
    return results


############################
# FISHER'S EXACT TEST ENRICHMENT
############################

def run_fisher(
    scores: pd.DataFrame,
    annotations: pd.DataFrame,
    p_threshold: float,
    min_genes: int,
    max_genes: int | None,
    max_genes_fraction: float | None,
    source_label: str | None = None
) -> pd.DataFrame:
    """
    Fisher's exact test enrichment using a p-value threshold to classify
    genes as significant vs non-significant.

    scores: DataFrame with columns gene, value (value = raw p, 0..1)
    annotations: DataFrame with gene, go_term, go_name, level
    """

    # Check that values look like p-values
    vals = scores["value"].astype(float)
    if (vals < 0).any() or (vals > 1).any():
        raise ValueError(
            "Fisher's exact test requires raw p-values in [0, 1]. "
            "Use --score-type raw with a p-value column."
        )

    merged = annotations.merge(scores, on="gene", how="inner")
    if merged.empty:
        raise ValueError("No overlapping genes between scores and annotations")

    merged = filter_categories_by_size(merged, min_genes, max_genes, max_genes_fraction)

    gene_pvals = merged.groupby("gene")["value"].mean()
    all_genes = set(gene_pvals.index)
    sig_genes = set(gene_pvals[gene_pvals <= p_threshold].index)
    nonsig_genes = all_genes - sig_genes

    logger.info(
        f"Fisher: p-threshold={p_threshold}, "
        f"{len(sig_genes)} significant genes, {len(nonsig_genes)} non-significant"
    )
    if len(sig_genes) == 0 or len(nonsig_genes) == 0:
        logger.warning("Fisher: all genes fall into one category; no test possible")
        return pd.DataFrame(
            columns=[
                "term", "name", "level", "nseqs",
                "sig_genes", "total_genes",
                "enrichment_ratio", "odds_ratio", "delta.rank",
                "pval", "p.adj"
            ]
        )

    categories = merged["go_term"].unique()
    label = source_label or getattr(scores, "source_file", "input")
    logger.info(f"Fisher: testing {len(categories)} GO categories from '{label}'")

    from scipy.stats import fisher_exact

    rows = []
    for term in categories:
        genes_in_term = set(merged.loc[merged["go_term"] == term, "gene"])
        if not genes_in_term:
            continue

        cat_sig = len(genes_in_term & sig_genes)
        cat_nonsig = len(genes_in_term & nonsig_genes)
        noncat_sig = len(sig_genes - genes_in_term)
        noncat_nonsig = len(nonsig_genes - genes_in_term)

        # 2x2 table:
        #          sig   non-sig
        # in term   a      b
        # not term  c      d
        table = [[cat_sig, cat_nonsig],
                 [noncat_sig, noncat_nonsig]]

        # fisher_exact can handle zeros; odds_ratio may be inf.
        odds_ratio, pval = fisher_exact(table, alternative="greater")

        total_in_term = cat_sig + cat_nonsig
        frac_sig_in_term = cat_sig / total_in_term if total_in_term > 0 else 0.0
        frac_sig_overall = len(sig_genes) / len(all_genes)
        enrichment_ratio = (
            frac_sig_in_term / frac_sig_overall if frac_sig_overall > 0 else 0.0
        )

        row = merged.loc[merged["go_term"] == term].iloc[0]
        rows.append({
            "term": term,
            "name": row["go_name"],
            "level": row["level"],
            "nseqs": len(genes_in_term),
            "sig_genes": cat_sig,
            "total_genes": len(genes_in_term),
            "enrichment_ratio": enrichment_ratio,
            "odds_ratio": odds_ratio,
            # delta.rank here is a placeholder to keep a similar interface to MWU
            "delta.rank": enrichment_ratio - 1.0,
            "pval": pval
        })

    if not rows:
        logger.warning("Fisher: no valid categories tested")
        return pd.DataFrame(
            columns=[
                "term", "name", "level", "nseqs",
                "sig_genes", "total_genes",
                "enrichment_ratio", "odds_ratio", "delta.rank",
                "pval", "p.adj"
            ]
        )

    results = pd.DataFrame.from_records(rows)
    results = add_bh_fdr(results)

    cols = [
        "term", "name", "level", "nseqs",
        "sig_genes", "total_genes",
        "enrichment_ratio", "odds_ratio",
        "delta.rank", "pval", "p.adj"
    ]
    results = results[cols].sort_values("p.adj").reset_index(drop=True)
    return results


############################
# RESULTS SAVING
############################

def save_results(
    results: pd.DataFrame,
    output_dir: Path,
    base_name: str,
    *,
    up_dir: Path | None = None,
    down_dir: Path | None = None,
    write_combined: bool = True
) -> None:
    """
    Save main result table as TSV.
    If a 'direction' column exists (bidirectional MWU), also save _UP and _DOWN subsets.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    df = results.copy()
    for col in df.select_dtypes(include=["object"]).columns:
        df[col] = df[col].astype(str).str.replace("\t", " ", regex=False)

    if write_combined:
        main_file = output_dir / f"{base_name}_results.tsv"
        logger.info(f"Saving results to '{main_file}'")
        df.to_csv(main_file, sep="\t", index=False)

    if "direction" in df.columns:
        for direction_value, subset in df.groupby("direction"):
            suffix = direction_value.upper()
            if suffix == "UP" and up_dir is not None:
                target_dir = up_dir
            elif suffix == "DOWN" and down_dir is not None:
                target_dir = down_dir
            else:
                target_dir = output_dir
            out_file = target_dir / f"{base_name}_{suffix}_results.tsv"
            subset.to_csv(out_file, sep="\t", index=False)
            logger.info(f"Saved {direction_value} results to '{out_file}'")

    # Quick summary
    if not df.empty:
        n_10 = (df["p.adj"] < 0.1).sum()
        n_05 = (df["p.adj"] < 0.05).sum()
        n_01 = (df["p.adj"] < 0.01).sum()
        logger.info(f"Terms with FDR < 0.10: {n_10}")
        logger.info(f"Terms with FDR < 0.05: {n_05}")
        logger.info(f"Terms with FDR < 0.01: {n_01}")


############################
# CLI
############################

class VerticalHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Formatter that prints each option on its own line in the usage block."""

    def _format_usage(self, usage, actions, groups, prefix):
        if prefix is None:
            prefix = "usage: "
        if usage:
            return super()._format_usage(usage, actions, groups, prefix)

        actions = [a for a in (actions or []) if a.help != argparse.SUPPRESS]
        if not actions:
            return f"{prefix}{self._prog}\n"

        required = []
        optional = []
        for action in actions:
            invocation = self._format_action_invocation(action)
            if action.required or not action.option_strings:
                required.append(invocation)
            else:
                optional.append(invocation)

        lines = [f"{prefix}{self._prog}"]
        for inv in required:
            lines.append(f"\n  {inv}")
        for inv in optional:
            lines.append(f"\n  [{inv}]")
        lines.append("\n")
        return "".join(lines)


def build_parser():
    parser = argparse.ArgumentParser(
        description="DEGOE: GO enrichment using MWU or Fisher's test",
        formatter_class=VerticalHelpFormatter
    )

    # Required
    parser.add_argument(
        "--annotations", required=True, metavar="FILE",
        help="Gene-to-GO annotation file (gene<TAB>GO:...;GO:...)"
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--scores", metavar="FILE",
        help="Single score file (e.g., DESeq2 results)"
    )
    input_group.add_argument(
        "--scores-dir", metavar="DIR",
        help="Directory containing multiple score files to analyze (non-recursive)"
    )

    # Method
    parser.add_argument(
        "--method", choices=["mwu", "fisher"], default="mwu",
        help="Enrichment method (default: mwu)"
    )
    parser.add_argument(
        "--namespace", choices=["BP", "MF", "CC"], default="BP",
        help="GO namespace (default: BP)"
    )

    # Score loading
    parser.add_argument(
        "--score-type",
        choices=[
            "raw", "pvalue", "wald",
            "log2foldchange", "signed_pvalue", "pval_signedlog2fc"
        ],
        default="pvalue",
        help="Score transformation (default: pvalue for MWU, raw for Fisher)"
    )
    parser.add_argument(
        "--score-column", metavar="COL",
        help="Specific score column (used for score_type='raw')"
    )
    parser.add_argument(
        "--pvalue-column", default="pvalue", metavar="COL",
        help="P-value column (default: pvalue)"
    )
    parser.add_argument(
        "--fold-change-column", default="log2FoldChange", metavar="COL",
        help="Fold-change column (default: log2FoldChange)"
    )
    parser.add_argument(
        "--respect-independent-filtering",
        action="store_true",
        help="Drop rows with NA in pvalue/padj columns (DESeq2 filtering)"
    )

    # MWU options
    parser.add_argument(
        "--test-direction",
        choices=["enriched", "up", "down", "both"],
        default="enriched",
        help="Direction for MWU (default: enriched)"
    )

    # Fisher options
    parser.add_argument(
        "--threshold", type=float, default=0.05,
        help="P-value threshold for Fisher's test (default: 0.05)"
    )

    # GO database and levels
    parser.add_argument(
        "--go-file", default="go.obo", metavar="FILE",
        help="GO OBO file (downloaded if missing, default: go.obo)"
    )
    parser.add_argument(
        "--min-level", type=int, default=None,
        help="Minimum GO hierarchy level to include"
    )
    parser.add_argument(
        "--max-level", type=int, default=None,
        help="Maximum GO hierarchy level to include"
    )

    # Category size filters
    parser.add_argument(
        "--min-genes", type=int, default=5,
        help="Minimum number of genes per GO term (default: 5)"
    )
    parser.add_argument(
        "--max-genes", type=int, default=None,
        help="Maximum number of genes per GO term (absolute)"
    )
    parser.add_argument(
        "--max-genes-fraction", type=float, default=None,
        help="Maximum fraction of total genes per GO term (e.g. 0.1 = 10%%)"
    )
    parser.add_argument(
        "--reduce-redundancy",
        action="store_true",
        help="Cluster GO categories based on shared genes and collapse redundant terms"
    )
    parser.add_argument(
        "--cluster-cut-height", type=float, default=0.25,
        help="Distance threshold for GO-term clustering when reducing redundancy (default: 0.25)"
    )

    # Output
    parser.add_argument(
        "--output-dir", default="./GO_Enrichment_results", metavar="DIR",
        help="Output directory (default: ./GO_Enrichment_results)"
    )
    parser.add_argument(
        "--output-prefix", default=None, metavar="PREFIX",
        help="Output prefix (default: derived from score filename)"
    )

    # Processing
    parser.add_argument(
        "--cpus", type=int, default=_n_cpus,
        help="Number of CPU cores for internal use (currently only affects BLAS env)"
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    params_file = Path(args.output_dir) / "parameters.txt"

    # Some basic consistency checks
    if args.method == "fisher":
        if args.score_type != "raw":
            logger.info("Overriding score_type to 'raw' for Fisher's exact test")
            args.score_type = "raw"
        if args.test_direction != "enriched":
            logger.info("Fisher's test ignores test_direction and always tests enrichment")
            args.test_direction = "enriched"

    if args.min_level is not None and args.max_level is not None:
        if args.min_level > args.max_level:
            logger.error("min-level cannot be greater than max-level")
            sys.exit(1)

    # Ensure GO file
    go_file = ensure_go_obo(args.go_file)
    go_hierarchy = GOHierarchy(go_file, args.namespace)

    # Output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.method == "mwu":
        if args.test_direction == "both":
            combined_dir = output_dir / "BOTH"
            up_dir = output_dir / "UP"
            down_dir = output_dir / "DOWN"
            for d in (combined_dir, up_dir, down_dir):
                d.mkdir(parents=True, exist_ok=True)
        elif args.test_direction == "up":
            up_dir = output_dir / "UP"
            up_dir.mkdir(parents=True, exist_ok=True)
            combined_dir = up_dir
            down_dir = None
        elif args.test_direction == "down":
            down_dir = output_dir / "DOWN"
            down_dir.mkdir(parents=True, exist_ok=True)
            combined_dir = down_dir
            up_dir = None
        else:
            combined_dir = output_dir / args.test_direction.upper()
            combined_dir.mkdir(parents=True, exist_ok=True)
            up_dir = down_dir = None
    else:
        combined_dir = output_dir / "FISHER"
        combined_dir.mkdir(parents=True, exist_ok=True)
        up_dir = down_dir = None

    # Load annotations and expand
    annotations, expanded_map = load_annotations(args.annotations, go_hierarchy)
    save_expanded_annotations(expanded_map, output_dir, args.annotations)

    # Level filtering
    annotations = filter_by_go_level(annotations, args.min_level, args.max_level)
    if args.reduce_redundancy:
        annotations = collapse_redundant_terms(annotations, args.cluster_cut_height)

    if args.scores:
        score_files = [Path(args.scores)]
    else:
        scores_dir = Path(args.scores_dir)
        if not scores_dir.is_dir():
            logger.error(f"Scores directory '{scores_dir}' does not exist or is not a directory")
            sys.exit(1)
        score_files = sorted(p for p in scores_dir.iterdir() if p.is_file())
        if not score_files:
            logger.error(f"Scores directory '{scores_dir}' contains no files")
            sys.exit(1)

    multiple_runs = len(score_files) > 1

    with open(params_file, "w") as pf:
        pf.write("DEGOE parameter log\n")
        pf.write(f"Timestamp: {datetime.now().isoformat()}\n")
        for key, value in vars(args).items():
            pf.write(f"{key}: {value}\n")
        pf.write("Score files:\n")
        for path in score_files:
            pf.write(f"  - {path}\n")

    for idx, score_path in enumerate(score_files, start=1):
        logger.info(
            f"Processing score file {idx}/{len(score_files)}: '{score_path}'"
        )
        scores = ScoreFileProcessor.load_scores(
            filepath=str(score_path),
            score_type=args.score_type,
            score_column=args.score_column,
            pvalue_column=args.pvalue_column,
            fold_change_column=args.fold_change_column,
            respect_independent_filtering=args.respect_independent_filtering
        )

        if args.method == "mwu":
            results = run_mwu(
                scores=scores,
                annotations=annotations,
                score_type=args.score_type,
                test_direction=args.test_direction,
                min_genes=args.min_genes,
                max_genes=args.max_genes,
                max_genes_fraction=args.max_genes_fraction,
                source_label=score_path.name
            )
        else:
            results = run_fisher(
                scores=scores,
                annotations=annotations,
                p_threshold=args.threshold,
                min_genes=args.min_genes,
                max_genes=args.max_genes,
                max_genes_fraction=args.max_genes_fraction,
                source_label=score_path.name
            )

        if args.output_prefix:
            if multiple_runs:
                base_name = f"{args.output_prefix}_{score_path.stem}"
            else:
                base_name = args.output_prefix
        else:
            base_name = f"{score_path.stem}_{args.method}_{args.namespace}"

        if args.method == "mwu" and args.test_direction == "both":
            save_results(
                results,
                combined_dir,
                base_name,
                up_dir=up_dir,
                down_dir=down_dir,
                write_combined=True
            )
        elif args.method == "mwu" and args.test_direction == "up" and up_dir is not None:
            save_results(
                results,
                up_dir,
                base_name,
                up_dir=up_dir,
                write_combined=False
            )
        elif args.method == "mwu" and args.test_direction == "down" and down_dir is not None:
            save_results(
                results,
                down_dir,
                base_name,
                down_dir=down_dir,
                write_combined=False
            )
        else:
            save_results(results, combined_dir, base_name)

        label = f" for {score_path.name}" if multiple_runs else ""
        if results.empty:
            print(f"\nNo enriched GO terms found{label}.")
        else:
            print(f"\nTop GO terms{label}:")
            display_cols = ["term", "name", "level", "nseqs", "pval", "p.adj"]
            if "direction" in results.columns:
                display_cols.insert(4, "direction")
            print(results[display_cols].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
