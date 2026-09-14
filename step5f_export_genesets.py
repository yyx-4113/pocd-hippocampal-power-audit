"""Export the four a priori mechanism gene sets (A-D) used in the bulk
meta-analysis to a citable supplementary CSV, so readers can verify the
"four lines down-regulated" conclusion. The gene members are defined once in
step2_bulk_integrate_pathway.py (SET_A..SET_D / PATHWAY_DEF); we exec only the
definition block (before the data-loading section) to avoid re-running the
pipeline.
"""
import os
import csv

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "step2_bulk_integrate_pathway.py")

src = open(SRC, encoding="utf-8").read()
# keep only the SET_A..SET_D / PATHWAY_DEF definition block.
# The definitions end right before "NDUF_PREFIX" (Complex I prefix note), which is
# the last statement before the data-loading section.
start = src.index("SET_A = [")
end = src.index("# Complex I 用前缀扩展")
block = src[start:end]
# PATHWAY_DEF is declared just after the Complex-I note; append it explicitly.
pdef_start = src.index("PATHWAY_DEF = {")
pdef_end = src.index("]", pdef_start) + 1
block = block + "\n" + src[pdef_start:pdef_end] + "\n"
ns = {}
exec(block, ns)

SETS = [
    ("A_complement_synapse_pruning", "complement / synaptic pruning"),
    ("B_DAM_microglia", "damage-associated microglia (DAM) / microglial activation"),
    ("C_mitochondrial_OXPHOS", "mitochondrial OXPHOS"),
    ("D_myelin_oligodendrocyte", "myelin / oligodendrocyte"),
]

OUT = os.path.join(HERE, "sc_out", "mechanism_genesets.csv")
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["set_id", "set_label", "gene"])
    for sid, label in SETS:
        for g in ns["PATHWAY_DEF"][sid]:
            w.writerow([sid, label, g])

print("wrote", OUT)
total = 0
for sid, label in SETS:
    n = len(ns["PATHWAY_DEF"][sid])
    total += n
    print(f"  {sid}: {n} genes")
print("  total:", total)
