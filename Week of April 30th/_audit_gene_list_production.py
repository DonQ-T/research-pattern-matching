"""End-to-end audit of how the shipped gene-list CSVs are produced.

Walks through the data structures and the production logic step-by-step,
printing concrete shapes and values at each stage. Intended as a one-pass
read-along to verify nothing surprising is happening.
"""
import os, pickle, numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.abspath(__file__))
os.chdir(ROOT)

# ===================================================================
# STAGE 1 — INPUTS:  what we load from disk
# ===================================================================
print('=' * 70)
print('STAGE 1 — INPUTS')
print('=' * 70)

# 1A. The pickle from Week of April 9 holds the 6-method scores per cluster
results = pickle.load(open('../Week of April 9th/metrics_results.pkl', 'rb'))
print(f'\n1A. metrics_results.pkl top-level keys: {list(results.keys())}')

constraint_patternsLT = results['constraint_patternsLT']
all_methodsLT = results['all_methodsLT']

# constraint_patternsLT: dict of cluster -> {"constraint": [...], possibly others}
print(f'\n    constraint_patternsLT keys: {list(constraint_patternsLT.keys())}')
print(f'    Example: constraint_patternsLT["clusterOneLT"] = {dict(constraint_patternsLT["clusterOneLT"])}')

# all_methodsLT: dict of method -> (scores_dict, ascending_bool)
print(f'\n    all_methodsLT keys: {list(all_methodsLT.keys())}')
pearson_scores_obj, pearson_asc = all_methodsLT['Pearson']
print(f'    Example: all_methodsLT["Pearson"] = (scores_dict, ascending={pearson_asc})')
print(f'    pearson_scores_obj is a dict with keys: {list(pearson_scores_obj.keys())}')
print(f'    pearson_scores_obj["clusterOneLT"] keys: {list(pearson_scores_obj["clusterOneLT"].keys())}')

pearson_clusterOne = pearson_scores_obj['clusterOneLT']['constraint']
print(f'\n    pearson_scores_obj["clusterOneLT"]["constraint"] is a {type(pearson_clusterOne).__name__}')
print(f'    shape: {pearson_clusterOne.shape}')
print(f'    head:')
print(pearson_clusterOne.head().to_string())
print(f'    => one Pearson r value per gene; index is gene name. ascending=False means higher r = better.')

# 1B. The gene expression matrices (log-transformed counts)
x = [0, 3, 6, 9]
gene_dataLT = {}
for c in constraint_patternsLT:
    df_ = pd.read_csv(f'../analysis data/gene_countsLT/{c}_annotated.csv', index_col=0)
    df_.columns = x
    gene_dataLT[c] = df_

clusterOne_genes = gene_dataLT['clusterOneLT']
print(f'\n1B. gene_dataLT["clusterOneLT"]:  shape = {clusterOne_genes.shape}')
print(f'    columns = {list(clusterOne_genes.columns)}  (timepoints)')
print(f'    head:')
print(clusterOne_genes.head().to_string())
print(f'    => rows are genes, columns are timepoints, values are LT counts.')

# ===================================================================
# STAGE 2 — DERIVED:  rawMAE is computed in the notebook
# ===================================================================
print('\n' + '=' * 70)
print('STAGE 2 — DERIVED:  rawMAE per gene per cluster')
print('=' * 70)

raw_mae = {}
for c in gene_dataLT:
    pat = np.array(constraint_patternsLT[c]['constraint'], dtype=float)
    mae_per_gene = np.mean(np.abs(gene_dataLT[c].values - pat), axis=1)
    raw_mae[c] = pd.Series(mae_per_gene, index=gene_dataLT[c].index)

print(f'\nrawMAE math (per gene g, pattern p):')
print(f'    rawMAE(g) = mean( |g[t] - p[t]| )  averaged over t in [0, 3, 6, 9]')
print(f'\nraw_mae["clusterTwoLT"]: shape = {raw_mae["clusterTwoLT"].shape}, head:')
print(raw_mae['clusterTwoLT'].head().to_string())
print(f'    => one MAE value per gene; lower = closer to pattern in raw LT space.')

# Hand-verify rawMAE for one gene
g0 = raw_mae['clusterTwoLT'].index[0]
gene_row = gene_dataLT['clusterTwoLT'].loc[g0].values
pat = np.array(constraint_patternsLT['clusterTwoLT']['constraint'])
manual = float(np.mean(np.abs(gene_row - pat)))
stored = float(raw_mae['clusterTwoLT'].loc[g0])
print(f'\n    Hand-verify gene {g0!r}:')
print(f'        gene values:    {gene_row}')
print(f'        pattern values: {pat}')
print(f'        |g - p|:        {np.abs(gene_row - pat)}')
print(f'        mean(|g - p|) = {manual:.6f}')
print(f'        stored value  = {stored:.6f}')
print(f'        match: {np.isclose(manual, stored)}')

# ===================================================================
# STAGE 3 — RANKING:  how the top-N is produced from the scores
# ===================================================================
print('\n' + '=' * 70)
print('STAGE 3 — RANKING')
print('=' * 70)

print('\nThe notebook helpers:')
print('''
    def get_scores(method, cluster):
        if method == "rawMAE":
            return raw_mae[cluster], True            # ascending=True (lower better)
        s, asc = all_methodsLT[method]
        return s[cluster]["constraint"], asc          # asc as set in the pkl

    def get_top_genes_plus(method, cluster, n=20):
        s, asc = get_scores(method, cluster)
        return s.nsmallest(n).index.tolist() if asc else s.nlargest(n).index.tolist()
''')

# Demo: top-5 Pearson on clusterOne
s, asc = all_methodsLT['Pearson']
sc = s['clusterOneLT']['constraint']
print(f'Demo: top-5 Pearson on clusterOne (ascending={asc}, so we take nlargest):')
top5 = sc.nlargest(5)
print(top5.to_string())

# Demo: top-5 rawMAE on clusterTwo
print(f'\nDemo: top-5 rawMAE on clusterTwo (ascending=True, so we take nsmallest):')
print(raw_mae['clusterTwoLT'].nsmallest(5).to_string())

# ===================================================================
# STAGE 4 — CSV PRODUCTION:  winner_full_ranking() builds the DataFrame
# ===================================================================
print('\n' + '=' * 70)
print('STAGE 4 — CSV PRODUCTION')
print('=' * 70)

print('''
For each cluster, winner_full_ranking(cluster, winner) does this:
    1. Look up the winner's score Series for that cluster (via get_scores).
    2. Sort by score (ascending if "lower better", descending if "higher better").
    3. Build a DataFrame:
         cluster   = short cluster name (e.g. "clusterOne")
         gene      = gene name from the sorted index
         rank      = 1..N (1 = best)
         score     = the metric value
         method    = winner method name
         direction = "min" or "max"
         avg_rank  = NaN for base methods (used for ensembles)
    4. Write to deliverables/<cluster>_ranking.csv.
''')

WINNERS = {
    'clusterOneLT':   'Pearson',
    'clusterTwoLT':   'rawMAE',
    'clusterThreeLT': 'rawMAE',
    'clusterFourLT':  'rawMAE',
}

# Read back what's actually in the CSVs and show
for c, w in WINNERS.items():
    short = c.replace('LT', '')
    csv = pd.read_csv(f'deliverables/{short}_ranking.csv')
    print(f'\n--- {short}_ranking.csv  ({w}) ---')
    print(f'  total rows: {len(csv):,}  (expected 17,856 = all genes)')
    print(f'  columns:    {list(csv.columns)}')
    print(f'  top 5 rows:')
    print(csv.head().to_string(index=False))
    print(f'  bottom 3 rows:')
    print(csv.tail(3).to_string(index=False))

# ===================================================================
# STAGE 5 — INVARIANTS:  things that must be true if the export is correct
# ===================================================================
print('\n' + '=' * 70)
print('STAGE 5 — INVARIANTS (sanity)')
print('=' * 70)

passed = []

for c, w in WINNERS.items():
    short = c.replace('LT', '')
    csv = pd.read_csv(f'deliverables/{short}_ranking.csv')

    # (1) row count
    n_genes = len(gene_dataLT[c])
    inv1 = len(csv) == n_genes
    passed.append((f'{short}: row count == 17,856', inv1))

    # (2) ranks are 1..N with no gaps
    inv2 = csv['rank'].tolist() == list(range(1, n_genes + 1))
    passed.append((f'{short}: ranks are 1..N exactly', inv2))

    # (3) gene set matches the input gene list
    inv3 = set(csv['gene']) == set(gene_dataLT[c].index)
    passed.append((f'{short}: gene set matches input cluster', inv3))

    # (4) scores monotone in the right direction
    if w == 'rawMAE':
        inv4 = (csv['score'].diff().dropna() >= -1e-12).all()  # non-decreasing (ascending)
        passed.append((f'{short}: rawMAE scores non-decreasing', inv4))
    else:  # Pearson, higher better -> non-increasing
        inv4 = (csv['score'].diff().dropna() <= 1e-12).all()
        passed.append((f'{short}: Pearson scores non-increasing', inv4))

    # (5) rank-1 gene matches the actual best by the metric
    if w == 'rawMAE':
        expected_best = raw_mae[c].idxmin()
    else:
        s, asc = all_methodsLT[w]
        sc = s[c]['constraint']
        expected_best = sc.idxmin() if asc else sc.idxmax()
    inv5 = csv.iloc[0]['gene'] == expected_best
    passed.append((f'{short}: rank-1 gene == idxmax/min of metric ({expected_best})', inv5))

print()
all_pass = True
for name, ok in passed:
    flag = 'PASS' if ok else 'FAIL'
    print(f'  [{flag}]  {name}')
    if not ok:
        all_pass = False

print(f'\nOVERALL: {"all invariants hold" if all_pass else "FAILURES detected"}')
