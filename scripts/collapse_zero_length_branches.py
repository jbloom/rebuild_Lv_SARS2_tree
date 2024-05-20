"""Collapse zero length branches."""


import Bio.Phylo


tree = Bio.Phylo.read(snakemake.input.tree, "newick")

to_prune = []
for node in tree.find_clades(order='postorder'):
    if node.branch_length <= snakemake.params.blmin and node != tree.root and (node.is_terminal() == False):
        to_prune.append(node)
print(f"Collapsing {len(to_prune)} nodes")
for node in to_prune:
    tree.collapse(node)

Bio.Phylo.write(tree, snakemake.output.tree, "newick")
