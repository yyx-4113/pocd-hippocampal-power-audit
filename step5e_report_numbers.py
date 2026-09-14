import pandas as pd, glob, os, numpy as np
cts = sorted(os.path.basename(f)[5:-4] for f in glob.glob('sc_out/pbcounts/voom_*.csv'))
# 各类型 Ttr 表达水平 (mean log2 CPM) + logFC + adjP
print('=== Ttr panel (TMM voom) ===')
print('%-16s %10s %10s %10s %10s' % ('celltype', 'mean_log2CPM', 'logFC', 'P', 'adjP'))
for ct in cts:
    d = pd.read_csv('sc_out/pbcounts/voom_%s.csv' % ct, index_col=0)
    if 'Ttr' in d.index:
        r = d.loc['Ttr']
        # 平均表达 = 平均 logCPM (limma 输出 AveExpr)
        av = r['AveExpr'] if 'AveExpr' in d.columns else np.nan
        print('%-16s %10.2f %10.2f %10.2e %10.3g' % (ct, av, r['logFC'], r['P.Value'], r['adj.P.Val']))

print()
print('=== 11 significant genes summary (TMM voom) ===')
sigs = []
for ct in cts:
    d = pd.read_csv('sc_out/pbcounts/voom_%s.csv' % ct, index_col=0)
    s = d[d['adj.P.Val'] < 0.05]
    for g, r in s.iterrows():
        sigs.append(dict(celltype=ct, gene=g, logFC=round(float(r['logFC']), 2),
                         AveExpr=round(float(r['AveExpr']), 1),
                         P=float(r['P.Value']), adjP=float(r['adj.P.Val'])))
sg = pd.DataFrame(sigs).sort_values(['celltype', 'adjP'])
sg.to_csv('sc_out/step5e_voom_significant_genes.csv', index=False)
print(sg.to_string(index=False))
print()
print('all logFC negative?', bool((sg.logFC < 0).all()))
print('Ttr count:', int((sg.gene == 'Ttr').sum()), '/', len(sg))
print('per celltype:', sg.celltype.value_counts().to_dict())

# 汇总
st = pd.read_csv('sc_out/step5e_limma_voom_stats.csv')
cal = pd.read_csv('sc_out/step5c_genomewide_calibration.csv').set_index('celltype')
mg = st.merge(cal[['p_lt_0.05']].rename(columns={'p_lt_0.05': 'welch_p_lt_0.05'}),
              left_on='celltype', right_index=True, how='left')
mg.to_csv('sc_out/step5e_method_comparison.csv', index=False)
print()
print('=== method comparison (limma-voom vs Welch/lognorm) ===')
print(mg[['celltype', 'n_genes_tested', 'n_p_lt_0.05', 'exp_p_lt_0.05',
          'n_adjp_lt_0.05', 'min_adjp', 'welch_p_lt_0.05']].to_string(index=False))
print()
print('TOTAL voom: tested=%d adjp<0.05=%d ; Welch: tested=173421 adjp<0.05=0'
      % (st['n_genes_tested'].sum(), int(st['n_adjp_lt_0.05'].sum())))
