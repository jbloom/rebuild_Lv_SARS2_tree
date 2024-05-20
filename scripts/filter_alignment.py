"""Filter alignment to match Lv et al filtering."""


import Bio.SeqIO
import Bio.Seq

import pandas as pd

lineages_to_drop = snakemake.params.lineages_to_drop
lineages_to_keep_cc_or_tt = snakemake.params.lineages_to_keep_cc_or_tt
mask_sites = snakemake.params.mask_sites
max_unmasked_n = snakemake.params.max_unmasked_n
drop_accessions = snakemake.params.drop_accessions

lineages = pd.read_csv(snakemake.input.metadata, sep="\t").set_index("accession")["lineage"].to_dict()

retained_seqs = []
for seq in Bio.SeqIO.parse(snakemake.input.alignment, "fasta"):
    if seq.id in drop_accessions:
        print(f"Dropping {seq.id} as a problematic accession")
        continue
    s = list(str(seq.seq))
    for i in mask_sites:
        s[i - 1] = "N"
    seq.seq = Bio.Seq.Seq("".join(s))
    lineage = lineages[seq.id]
    if lineage in lineages_to_drop:
        print(f"Dropping {seq.id} because it is {lineage=}")
        continue
    nt8782 = seq[8781]
    nt28144 = seq[28143]
    nts = f"{nt8782}{nt28144}"
    if nts in {"CC", "TT"}:
        if lineage not in lineages_to_keep_cc_or_tt:
            print(f"Dropping {seq.id} as {nts} of {lineage=}")
            continue
    elif nts not in {"CT", "TC"}:
        print(f"Dropping {seq.id} because it has {nts} and {lineage=}")
        continue
    unmasked_n = sum(nt in {"N", "n"} for (i, nt) in enumerate(seq.seq) if i + 1 not in mask_sites)
    if unmasked_n > max_unmasked_n:
        print(f"Dropping {seq.id} ({lineage=}) because it has {unmasked_n} ambiguous nucleotides at unmasked sites")
        continue
    retained_seqs.append(seq)

print(f"Retaining {len(retained_seqs)} sequences.")
Bio.SeqIO.write(retained_seqs, snakemake.output.alignment, "fasta")
