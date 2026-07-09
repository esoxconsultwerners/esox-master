import json, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'analysis' / 'c1_degeneration'))
sys.path.insert(0, str(ROOT / 'analysis' / 'c2_composition'))
import c1_common as cc
import c2_composition as c2

OUT = Path(__file__).resolve().parent
BAND_COLS = [f'refl_{b}' for b in cc.BANDS]
ri = cc.BANDS.index(cc.REF_NM)
RNG = np.random.default_rng(0)

c2._ensure_matchers()
clf, classes = c2._FULL
classes = list(classes)
_, meta, _ = cc.load_relab_resampled()
Xt = c2.relab_features(meta, cc.BANDS, cc.REF_NM)
ok = np.isfinite(Xt).all(1)
Xt, meta = Xt[ok], meta[ok].reset_index(drop=True)
yt = meta['relab_group'].map(c2.merge_label).to_numpy()
mnames = meta['meteorite_name'].to_numpy()
sep, deg, sepm = c2.separable_sets()
deg_merged = {c2.merge_label(x) for x in deg}

core = pd.read_parquet(ROOT / 'data/final/esox_master_core.parquet')
g = pd.read_parquet(cc.GASP_CAT)
X = g.set_index('number_mp').loc[core['number_mp'], BAND_COLS].to_numpy()
status = core['analog_status'].to_numpy()
top = core['analog_group_top'].astype(object).to_numpy()
cx = core['taxon_esox'].astype(object).to_numpy()
conf = core['analog_group_conf'].to_numpy()
okm = status == 'ok'

raw = clf.predict_proba(X)
raw_top = np.array([classes[i] for i in raw.argmax(1)])

res = {}

print('==== CHECK 1: analog_top (confident) vs taxon_esox complex ====')
ct = pd.crosstab(pd.Series(top[okm], name='analog_top'),
                 pd.Series(cx[okm], name='taxon_esox'))
print(ct.to_string())
ct.to_csv(OUT / 'check1_contingency.csv')
res['n_confident'] = int(okm.sum())
res['confident_CM'] = int(((top == 'CM') & okm).sum())
res['confident_CM_Scomplex'] = int(((top == 'CM') & okm & (cx == 'S')).sum())
res['confident_CM_Ccomplex'] = int(((top == 'CM') & okm & (cx == 'C')).sum())

print()
print('==== raw (flat-prior) argmax across all', len(core), 'objects ====')
rawvc = pd.Series(raw_top).value_counts()
print(rawvc.to_string())
res['raw_argmax_CM'] = int((raw_top == 'CM').sum())
res['raw_argmax_CM_frac'] = round(float((raw_top == 'CM').mean()), 4)

print()
print('==== CHECK 2: template count vs wins ====')
tcount = pd.Series(yt).value_counts()
nmet = pd.Series(mnames).groupby(yt).nunique() if False else        meta.assign(m=yt).groupby('m')['meteorite_name'].nunique()
rows = []
for c in sorted(set(yt)):
    rows.append({'group': c,
                 'n_templates': int(tcount.get(c, 0)),
                 'n_meteorites': int(nmet.get(c, 0)),
                 'times_top': int((top == c).sum()),
                 'times_confident': int(((top == c) & okm).sum()),
                 'in_separable': c in sepm})
t2 = pd.DataFrame(rows).sort_values('n_templates', ascending=False)
print(t2.to_string(index=False))
t2.to_csv(OUT / 'check2_template_vs_wins.csv', index=False)
sub = t2[t2['in_separable']]
r_all = float(np.corrcoef(t2['n_templates'], t2['times_confident'])[0, 1])
r_sep = float(np.corrcoef(sub['n_templates'], sub['times_confident'])[0, 1])
res['corr_templates_vs_confident_all'] = round(r_all, 3)
res['corr_templates_vs_confident_sep'] = round(r_sep, 3)
print('Pearson r (n_templates vs times_confident): all %.3f  separable %.3f' % (r_all, r_sep))
fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(t2['n_templates'], t2['times_confident'], color='#4059ad')
for _, rr in t2.iterrows():
    ax.annotate(rr['group'], (rr['n_templates'], rr['times_confident']), fontsize=7)
ax.set_xlabel('n_templates (RELAB spectra)')
ax.set_ylabel('times chosen confident')
ax.set_title('C2 template-density bias check (r_all=%.2f)' % r_all)
fig.tight_layout()
fig.savefig(OUT / 'check2_templates_vs_wins.png', dpi=140, metadata={'Date': None})
plt.close(fig)

print()
print('==== CHECK 3: CM winning on centroid or on nearest-template? ====')
cents = {c: Xt[yt == c].mean(0) for c in set(yt)}
smask = okm & (cx == 'S') & (top == 'CM')
idx = np.where(smask)[0][:50]
res['n_Scomplex_confCM_sample'] = int(len(idx))
cm_centroid_closest = 0
cm_nn_only = 0
nc_pred_cm = 0
rows3 = []
for j in idx:
    x = X[j]
    cd = {c: float(np.linalg.norm(x - cents[c])) for c in cents}
    nnd = {c: float(np.min(np.linalg.norm(Xt[yt == c] - x, axis=1))) for c in cents}
    cent_closest = min(cd, key=cd.get)
    nn_closest = min(nnd, key=nnd.get)
    if cent_closest == 'CM':
        cm_centroid_closest += 1
    elif nn_closest == 'CM':
        cm_nn_only += 1
    if cent_closest == 'CM':
        nc_pred_cm += 1
    rows3.append({'number_mp': int(core['number_mp'].iloc[j]),
                  'raw_p_CM': round(float(raw[j][classes.index('CM')]), 3),
                  'cent_closest': cent_closest, 'cent_d_CM': round(cd['CM'], 3),
                  'cent_d_min': round(cd[cent_closest], 3),
                  'nn_closest': nn_closest, 'nn_d_CM': round(nnd['CM'], 3),
                  'nn_d_min': round(nnd[nn_closest], 3)})
d3 = pd.DataFrame(rows3)
print(d3.head(20).to_string(index=False))
d3.to_csv(OUT / 'check3_mechanism.csv', index=False)
res['sample_CM_centroid_closest'] = cm_centroid_closest
res['sample_CM_nearest_template_only'] = cm_nn_only
res['sample_nearestcentroid_predicts_CM'] = nc_pred_cm
print('of %d S-complex confident-CM: CM-centroid-closest %d ; CM-only-nearest-single-template %d'
      % (len(idx), cm_centroid_closest, cm_nn_only))

print()
print('==== CHECK 4: feature contrast / normalization ====')
contrast = np.nanstd(X, axis=1)
res['contrast_confCM_median'] = round(float(np.nanmedian(contrast[(top == 'CM') & okm])), 4)
res['contrast_other_median'] = round(float(np.nanmedian(contrast[~((top == 'CM') & okm)])), 4)
res['contrast_all_median'] = round(float(np.nanmedian(contrast)), 4)
print('median band-std (contrast): confident-CM %.4f  rest %.4f  all %.4f'
      % (res['contrast_confCM_median'], res['contrast_other_median'], res['contrast_all_median']))
tmean = Xt.mean(0)
gmean = np.nanmean(X, axis=0)
print('per-band RELAB-template mean vs Gaia-core mean:')
for b, tm, gm in zip(cc.BANDS, tmean, gmean):
    print('  %4d  relab %.3f  gaia %.3f  diff %+.3f' % (b, tm, gm, gm - tm))
res['relab_1034_mean'] = round(float(tmean[-1]), 3)
res['gaia_1034_mean'] = round(float(gmean[-1]), 3)
cmtemp = Xt[yt == 'CM']
res['CM_template_1034_mean'] = round(float(cmtemp.mean(0)[-1]), 3)
res['CM_template_contrast_median'] = round(float(np.median(np.std(cmtemp, axis=1))), 4)

print()
print('==== CHECK 5: prior realized effect on CM ====')
cmS = okm & (cx == 'S') & (top == 'CM')
ci = classes.index('CM')
raw_cm = raw[cmS][:, ci]
post_cm = []
for j in np.where(cmS)[0]:
    post, used = c2.apply_prior(raw[j], np.array(classes), cx[j])
    post_cm.append(post[ci])
post_cm = np.array(post_cm)
res['Scomplex_confCM_raw_pCM_mean'] = round(float(raw_cm.mean()), 3)
res['Scomplex_confCM_post_pCM_mean'] = round(float(post_cm.mean()), 3)
res['Scomplex_confCM_raw_pCM_median'] = round(float(np.median(raw_cm)), 3)
print('S-complex confident-CM: raw p(CM) mean %.3f  post-prior p(CM) mean %.3f (prior downweights CM x%.2f but still wins)'
      % (raw_cm.mean(), post_cm.mean(), c2.PRIOR_LOW))
allS = okm & (cx == 'S')
res['frac_confident_Scomplex_that_are_CM'] = round(float(((top == 'CM') & allS).sum() / max(allS.sum(), 1)), 3)

(OUT / 'diag_numbers.json').write_text(json.dumps(res, indent=2))
print()
print('SAVED diag_numbers.json:')
print(json.dumps(res, indent=2))
