# Which protein drives NCAN incorporation? Exploratory analyses on existing images.
# Called automatically by run_all.py
# (1) Is the BCAN effect on NCAN explained by bigger ACAN condensates? -> compare at matched condensate size
# (2) Within each image, does NCAN content per condensate track ACAN or BCAN?  -> per-FOV Spearman / partial correlation
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import recruitment_quant as rq
from scipy import ndimage as ndi
from scipy.stats import spearmanr, rankdata
from skimage import measure
from skimage.registration import phase_cross_correlation as pcc
p = dict(rq.DEFAULTS)
sheet_path = sys.argv[1] if len(sys.argv) > 1 else 'triple_sheet.csv'
out_dir = sys.argv[2] if len(sys.argv) > 2 else 'results_triple'
sheet = pd.read_csv(sheet_path)
def sh(ref, x): return ndi.shift(x, pcc(ref, x, upsample_factor=10, normalization=None)[0], order=1, mode='nearest')
rows = []
for _, r in sheet.iterrows():
    A, _ = rq.load_channel(r.scaffold_path); C, _ = rq.load_channel(r.client_path); C = sh(A, C)
    other = None
    if r.condition == 'Triple|NCAN':
        bpath = r.other_path if 'other_path' in r and isinstance(r.other_path, str) else r.client_path.replace('ch02', 'ch00')
        other, _ = rq.load_channel(bpath); other = sh(A, other)
    lab, _, _ = rq.segment(A, p, 'otsu'); _, fullC, _ = rq.segment(C, p, 'otsu')
    for o in measure.regionprops(lab):
        y, x = o.coords[:, 0], o.coords[:, 1]
        rows.append(dict(condition=r.condition, replicate=r.replicate, fov=r.image_id, obj=o.label, area=o.area,
                         acan=A[y, x].mean(), acan_sat=(A[y, x] >= 255).mean(),
                         client=C[y, x].mean(), client_cover=fullC[y, x].mean(),
                         bcan=(other[y, x].mean() if other is not None else np.nan)))
d = pd.DataFrame(rows); d.to_csv(os.path.join(out_dir, 'driver_analysis_objects.csv'), index=False)

print('=== (1) NCAN in ACAN condensates of matched size: ACAN+NCAN vs three-protein mix ===')
n = d[d.condition.isin(['ACAN+NCAN', 'Triple|NCAN'])].copy()
n['size_bin'] = pd.cut(n.area, [9, 50, 150, 400, 10000], labels=['9-50', '50-150', '150-400', '>400'])
tab = n.groupby(['size_bin', 'condition'], observed=True).agg(n=('obj', 'size'), ncan_area_cover=('client_cover', 'mean')).unstack()
print(tab.round(2))
fovbin = n.groupby(['condition', 'fov', 'size_bin'], observed=True).client_cover.mean().reset_index()
for b in ['50-150', '150-400']:
    a = fovbin[(fovbin.size_bin == b) & (fovbin.condition == 'ACAN+NCAN')].client_cover
    t = fovbin[(fovbin.size_bin == b) & (fovbin.condition == 'Triple|NCAN')].client_cover
    from scipy.stats import mannwhitneyu
    print(f'  size {b} px: FOV means double {a.mean():.2f} (n={len(a)}) vs triple {t.mean():.2f} (n={len(t)}), MWU P = {mannwhitneyu(a, t).pvalue:.3g}')

print('\n=== (2) Within-image correlation across condensates (per-FOV Spearman, then median over FOVs) ===')
def fov_rho(df, a, b):
    out = []
    for f, g in df.groupby('fov'):
        if len(g) >= 8: out.append(spearmanr(g[a], g[b])[0])
    return np.array(out)
def partial_rho(g, x, y, z):          # Spearman partial correlation of x,y given z
    rx, ry, rz = rankdata(g[x]), rankdata(g[y]), rankdata(g[z])
    res = lambda u, v: u - np.polyval(np.polyfit(v, u, 1), v)
    return spearmanr(res(rx, rz), res(ry, rz))[0]
big = lambda df: df[df.area >= 50]
for cond, lab in [('ACAN+BCAN', 'BCAN vs ACAN (ACAN+BCAN)'), ('ACAN+NCAN', 'NCAN vs ACAN (ACAN+NCAN)')]:
    r = fov_rho(big(d[d.condition == cond]), 'client', 'acan'); print(f'  {lab}: median rho {np.median(r):.2f}  (FOV range {r.min():.2f} to {r.max():.2f})')
tri = big(d[d.condition == 'Triple|NCAN'])
uns = tri[tri.acan_sat < 0.05]
print(f'  Three-protein mix, condensates >= 50 px: {len(tri)}; with <5% saturated ACAN pixels: {len(uns)}')
for label, df in [('all', tri), ('unsaturated only', uns)]:
    res = []
    for f, g in df.groupby('fov'):
        if len(g) < 8: continue
        res.append(dict(fov=f, n=len(g), rho_N_A=spearmanr(g.client, g.acan)[0], rho_N_B=spearmanr(g.client, g.bcan)[0],
                        rho_B_A=spearmanr(g.bcan, g.acan)[0],
                        partial_N_A_given_B=partial_rho(g, 'client', 'acan', 'bcan'),
                        partial_N_B_given_A=partial_rho(g, 'client', 'bcan', 'acan')))
    res = pd.DataFrame(res)
    if len(res): print(f'  [{label}]\n', res.round(2).to_string(index=False), '\n  medians:', res.drop(columns=['fov', 'n']).median().round(2).to_dict())
