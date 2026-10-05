#!/usr/bin/env python3
"""Declarative tables describing the study's clinical attributes and assay profiles.

Kept beside `cbioportal_study` so the assembly module stays readable; these are
data, not logic, and changing a description here changes only the emitted files.
"""

MAX_CLONE_ENTITY = 63
N_TREES = 5

CLINICAL_ATTRIBUTES = [
    ("PATIENT_ID", "Patient Identifier", "Patient identifier", "STRING", "1"),
    ("SAMPLE_ID", "Sample Identifier", "Sample identifier", "STRING", "1"),
    (
        "CCF_WEIGHTED_NEOANTIGEN_LOAD",
        "CCF-weighted neoantigen load",
        "Expected neoantigen load of a tumor cell, top-scoring tree",
        "NUMBER",
        "1",
    ),
    (
        "CCF_WEIGHTED_TMB",
        "CCF-weighted missense count",
        "Expected missense count of a tumor cell, top-scoring tree",
        "NUMBER",
        "1",
    ),
    (
        "CCF_WEIGHTED_FITNESS",
        "CCF-weighted clone fitness",
        "Prevalence-weighted clone fitness, top-scoring tree",
        "NUMBER",
        "1",
    ),
    (
        "DOMINANT_CLONE_FITNESS",
        "Dominant clone fitness",
        "Fitness of the most prevalent clone, top-scoring tree",
        "NUMBER",
        "1",
    ),
    (
        "DOMINANT_CLONE_NEOANTIGEN_LOAD",
        "Dominant clone neoantigen load",
        "Neoantigen load of the most prevalent clone, top-scoring tree",
        "NUMBER",
        "1",
    ),
    ("N_CLONES", "Number of clones", "Tumor clones in the top-scoring tree", "NUMBER", "1"),
    (
        "HAS_TRUNCAL_CLONE",
        "Has truncal clone",
        "Whether the top-scoring tree has a shared trunk",
        "STRING",
        "1",
    ),
    (
        "TRUNCAL_NEOANTIGEN_LOAD",
        "Truncal neoantigen load",
        "Neoantigen load of the MRCA, top-scoring tree",
        "NUMBER",
        "1",
    ),
    ("TOTAL_NEOANTIGEN_LOAD", "Total neoantigen load", "Neoantigens called in this sample", "NUMBER", "1"),
    (
        "UNASSIGNED_NEOANTIGEN_LOAD",
        "Unassigned neoantigen load",
        "Neoantigens whose mutation is in no clone, so clonality is unknown",
        "NUMBER",
        "1",
    ),
    (
        "SUBCLONAL_NEOANTIGEN_LOAD",
        "Subclonal neoantigen load",
        "Total minus truncal minus unassigned",
        "NUMBER",
        "1",
    ),
    ("MAX_CLONE_FITNESS", "Max clone fitness", "Highest clone fitness, top-scoring tree", "NUMBER", "1"),
    ("MAX_CLONE_F_P", "Max clone F_P", "Highest F_P value, top-scoring tree", "NUMBER", "1"),
    # Priority 0 hides the Study View chart: log-likelihood scales with mutation
    # count (-2420 at 255 mutations vs -30433 at 336), so a cross-sample
    # histogram of it is meaningless. Kept as an inspectable per-sample value.
    (
        "TREE_LOGLIK",
        "Tree log-likelihood (top-scoring tree)",
        "PhyloWGS log-likelihood; not comparable across samples",
        "NUMBER",
        "0",
    ),
    ("EFFECTIVE_N", "Effective N", "Effective population size from the pipeline", "NUMBER", "1"),
]

# Entity meta properties: one row per neoantigen, ordered for cBioPortal portal file.
NEOANTIGEN_META_PROPERTIES = (
    "NAME", "DESCRIPTION", "GENE", "MUTATION_ID", "SEQUENCE", "WT_SEQUENCE", "HLA", "MUTATED_POSITION",
)

NEOANTIGEN_CLONE_PROFILE = (
    "neoantigen_clone",
    "Neoantigen clone",
    "Clone carrying the neoantigen's mutation in this candidate tree; NA if unassigned.",
)

# Trailing field is value_sort_order. Kd and KdWT are affinities in nM where
# LOWER is the stronger binder, so they sort ASC; DESC there would rank the
# weakest binders as most important in the waterfall plot and OncoPrint tooltip.
NEOANTIGEN_PROFILES = [
    ("neoantigen_quality", "quality", "Neoantigen quality", "Neoantigen fitness from NeoantigenEditing.", "DESC"),
    (
        "neoantigen_kd",
        "Kd",
        "Neoantigen Kd",
        "netMHCpan binding affinity of the mutant peptide (nM).",
        "ASC",
    ),
    (
        "neoantigen_kdwt",
        "KdWT",
        "Neoantigen Kd (WT)",
        "netMHCpan binding affinity of the wild-type peptide (nM).",
        "ASC",
    ),
    ("neoantigen_r", "R", "Neoantigen R", "Similarity of the mutant peptide to IEDB peptides.", "DESC"),
    ("neoantigen_logc", "logC", "Neoantigen logC", "Log cross-reactivity.", "DESC"),
    ("neoantigen_loga", "logA", "Neoantigen logA", "Log amplitude.", "DESC"),
]

# The tidytree TSVs are the durable artifact: treeio/ggtree read them directly
# and they land in Databricks as flat tables, with every portal file a projection
# of them. Each entry is (filename, sample key, columns). Columns absent from a
# row render empty: `tilde_x` exists only on recurrent samples.
TIDY_TABLES = [
    (
        "tree_nodes.tsv",
        "nodes",
        [
            "sample_id",
            "tree_idx",
            "clone_id",
            "parent",
            "X",
            "x",
            "TMB",
            "neoantigen_load",
            "NA_Mut",
            "F_I",
            "F_P",
            "new_x",
            "tilde_x",
        ],
    ),
    ("tree_scores.tsv", "scores", ["sample_id", "tree_idx", "loglik"]),
    (
        "neoantigens.tsv",
        "neoantigens",
        [
            "sample_id",
            "neoantigen_id",
            "mutation_id",
            "gene",
            "HLA_gene_id",
            "sequence",
            "WT_sequence",
            "mutated_position",
            "Kd",
            "KdWT",
            "R",
            "logC",
            "logA",
            "quality",
        ],
    ),
    ("mutation_clones.tsv", "mutation_clones", ["sample_id", "tree_idx", "mutation_id", "clone_id"]),
    ("mutations.tsv", "mutations", ["sample_id", "mutation_id", "gene", "missense"]),
]

# Hidden like the other clone profiles, and sorted DESC: these are log-likelihoods,
# where the best-scoring tree is the largest value.
TREE_SCORE_PROFILE = (
    "clone_tree_score",
    "Candidate tree log-likelihood",
    "PhyloWGS log-likelihood per candidate tree, best first.",
)

CLONE_PROFILES = [
    ("clone_parent", "parent", "Clone parent", "Parent clone id; -1 marks the root."),
    ("clone_prevalence", "X", "Clone prevalence", "Inclusive clonal prevalence."),
    ("clone_prevalence_excl", "x", "Clone prevalence (exclusive)", "Exclusive clonal prevalence."),
    ("clone_tmb", "TMB", "Clone missense count", "Cumulative missense count along the lineage."),
    (
        "clone_neoantigen_load",
        "neoantigen_load",
        "Clone neoantigen load",
        "Cumulative neoantigen count along the lineage.",
    ),
    ("clone_na_mut", "NA_Mut", "Clone neoantigenic mutations", "Cumulative neoantigen-producing mutations."),
    ("clone_fitness", "F_I", "Clone fitness", "Clone fitness from NeoantigenEditing."),
    ("clone_f_p", "F_P", "Clone F_P", "F_P from NeoantigenEditing; integer, not a flag."),
]
