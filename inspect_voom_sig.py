import pandas as pd, glob, json, os
sigs = []
for f in sorted(glob.glob('sc_out/pbcounts/voom_*.csv')):
    ct = os.path.basename(f)[5:-4]
    d = pd.read_csv(f, index_col=0)
    s = d[d['adj.P.Val'] < 0.05]
    for g, r in s.iterrows():
        sigs.append((ct, g, round(r['logFC'], 2), r['P.Value'], r['adj.P.Val']))
print('=== significant genes (limma-voom, adj.P.Val<0.05) : %d ===' % len(sigs))
for ct, g, fc, p, a in sigs:
    print('  %-16s %-10s logFC=%7.2f  P=%.2e  adjP=%.2e' % (ct, g, fc, p, a))

pw = json.load(open('step2_pathway_sets.json', encoding='utf-8'))
genes = [g for _, g, _, _, _ in sigs]
print()
print('=== overlap with 4 mechanism lines ===')
hit = False
for k, v in pw.items():
    ov = [g for g in genes if g in set(v)]
    if ov:
        print('  ', k, ov); hit = True
print('  none' if not hit else '')

print()
print('=== Ttr / high-abundance choroid-plexus markers across cell types ===')
for f in sorted(glob.glob('sc_out/pbcounts/voom_*.csv')):
    ct = os.path.basename(f)[5:-4]
    d = pd.read_csv(f, index_col=0)
    for g in ['Ttr', 'Nov', 'Wfdc17']:
        if g in d.index:
            r = d.loc[g]
            print('  %-16s %-7s logFC=%7.2f  P=%.2e  adjP=%.3g' % (ct, g, r['logFC'], r['P.Value'], r['adj.P.Val']))

print()
print('=== summary: how many cell types have >=1 signal, and are they concordant? ===')
import collections
c = collections.Counter(ct for ct, _, _, _, _ in sigs)
print('  per cell type:', dict(c))
