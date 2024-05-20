"""``snakemake`` file that runs pipeline."""

# some configuration
reference = "MN908947.3"  # Wuhan-Hu-1 reference
outgroup = "EPI_ISL_406030"  # Lv et al roots at this A0 lineage in Fig 3A
lineages_to_drop = ["CC", "T-"]  # Lv et al drop sequences of these lineages (2.6 of methods)
lineages_to_keep_cc_or_tt = ["B0"]  # Lv et al drop CC / TT genotypes not in B0 (2.6 of methods)
mask_sites = list(range(1, 223)) + list(range(29700, 29904))  # mask sites 1-222 and 29,700-29,903 as done by Lv et al
max_unmasked_n = 1000 - len(mask_sites)  # max number of N nucleotides outside masked sites
drop_accessions = []  # problematic accessions not caught by other filters
color_lineages = {
    # same colors used in Fig 3A of Lv et al 
    "A0": "#800080",
    "A": "#0000FF",
    "B0": "#028A0F",
    "B": "#90EE90",
    "B-B.1": "#FFFF00",
    "B.1": "#FFA500",
    "B.4": "#FFC0CB",
}


rule all:
    """Target rule w final outputs."""
    input:
        "auspice/rebuild_Lv_SARS2_tree.json",


rule prep_sequences_and_metadata:
    """Prepare sequences and metadata for ``augur``."""
    input:
        seqs=["data/genbank_sequences.fa", "data/gisaid_sequences.fa"],
        accession_info="data/all_accessions.csv",
    output:
        sequences="results/unaligned_sequences.fa",
        metadata="results/metadata_w_fixed_accessions.tsv",
    log:
        notebook="results/prep_sequences_and_metadata.ipynb",
    notebook:
        "notebooks/prep_sequences_and_metadata.py.ipynb"


rule align:
    """Make an alignment including the outgroup."""
    input:
        unaligned_seqs=rules.prep_sequences_and_metadata.output.sequences,
    output:
        alignment="results/alignment_unfiltered.fa",
    params:
        reference=reference,
    threads: 4
    shell:
        """
        augur align \
            -s {input.unaligned_seqs} \
            -o {output.alignment} \
            --reference-name {params.reference} \
            --fill-gaps \
            --nthreads {threads}
        """


rule filter_alignment:
    """Filter the alignment as done by Lv et al."""
    input:
        alignment=rules.align.output.alignment,
        metadata=rules.prep_sequences_and_metadata.output.metadata,
    output:
        alignment="results/alignment_filtered.fa",
    params:
        lineages_to_drop=lineages_to_drop,
        lineages_to_keep_cc_or_tt=lineages_to_keep_cc_or_tt,
        mask_sites=mask_sites,
        max_unmasked_n=max_unmasked_n,
        drop_accessions=drop_accessions,
    script:
        "scripts/filter_alignment.py"


rule build_tree:
    """Build the phylogenetic tree."""
    input:
        alignment=rules.filter_alignment.output.alignment,
    output:
        tree="results/initial_tree.nwk",
        exclude_sites="results/exclude_sites.txt",
    params:
        outgroup=outgroup,
        exclude_sites="\n".join(map(str, mask_sites)),
        blmin=0.0000001,  # minimum branch length for zero-length branches
    shell:
        """
        printf "{params.exclude_sites}" > {output.exclude_sites}
        augur tree \
            -a {input.alignment} \
            --method iqtree \
            --tree-builder-args "\-o {params.outgroup} \-ninit 100 \-n 5 \-me 0.01 \-seed 1 \-blmin {params.blmin}" \
            --output {output.tree}
        """


rule collapse_zero_length_branches:
    """Collapse zero-length branches."""
    input:
        tree=rules.build_tree.output.tree,
    output:
        tree="results/collapsed_zero_branches_tree.nwk",
    params:
        blmin=rules.build_tree.params.blmin,
    script:
        "scripts/collapse_zero_length_branches.py"


rule refine_tree:
    """Run the ``augur`` refine command."""
    input:
        alignment=rules.filter_alignment.output.alignment,
        tree=rules.collapse_zero_length_branches.output.tree,
        metadata=rules.prep_sequences_and_metadata.output.metadata,
    output:
        tree="results/refined_tree.nwk",
        node_data="results/branch_lengths.json",
    shell:
        """
        augur refine \
            -a {input.alignment} \
            -t {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns accession \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --keep-root \
            --verbosity 2
        """

rule ancestral:
    """Mutations and ancestral sequences."""
    input:
        tree=rules.refine_tree.output.tree,
        alignment=rules.filter_alignment.output.alignment,
    output:
        node_data="results/nt-muts.json",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --keep-ambiguous
        """


rule translate:
    """Call protein mutations."""
    input:
        tree=rules.refine_tree.output.tree,
        ancestral_sequences=rules.ancestral.output.node_data,
        reference_sequence="data/reference_seq.gb",
    output:
        node_data="results/aa-muts.json",
    params:
        genes=["ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b"],
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.ancestral_sequences} \
            --reference-sequence {input.reference_sequence} \
            --output-node-data {output.node_data} \
            --genes {params.genes}
        """


rule add_jointWHO_metadata:
    """Add metadata about sequences from 2019 onset per joint WHO-China report."""
    input:
        metadata=rules.prep_sequences_and_metadata.output.metadata,
        jointWHO_metadata="data/jointWHO_sequences.csv",
    output:
        metadata="results/metadata_w_jointWHO.tsv",
    log:
        notebook="results/add_jointWHO_metadata.ipynb",
    notebook:
        "notebooks/add_jointWHO_metadata.py.ipynb"
        

rule export_tree:
    """Export the tree to JSON for Nextstrain auspice."""
    input:
        tree=rules.refine_tree.output.tree,
        node_data=rules.refine_tree.output.node_data,
        nt_muts=rules.ancestral.output.node_data,
        aa_muts=rules.translate.output.node_data,
        metadata=rules.add_jointWHO_metadata.output.metadata,
        description="data/description.md"
    output:
        colors="results/export_tree_colors.tsv",
        tree="auspice/rebuild_Lv_SARS2_tree.json",
    params:
        colors="\n".join(f"lineage\t{lineage}\t{color}" for (lineage, color) in color_lineages.items()),
    shell:
        """
        printf "{params.colors}" > {output.colors}
        augur export v2 \
            --tree {input.tree} \
            --output {output.tree} \
            --node-data {input.node_data} {input.nt_muts} {input.aa_muts} \
            --metadata {input.metadata} \
            --color-by-metadata lineage joint_WHO_report_2019_case \
            --metadata-columns strain date \
            --metadata-id-columns accession \
            --include-root-sequence-inline \
            --colors {output.colors} \
            --panels tree \
            --title 'Phylogeny from "Evolutionary trajectory of diverse SARS-CoV-2 variants at the beginning of COVID-19 outbreak" study by group of Yong-Zhen Zhang' \
            --description {input.description}
        """
