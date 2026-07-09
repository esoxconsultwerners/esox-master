"""C2.5 per-group aggregate-likelihood matcher (CM-bias fix).

Replaces the RandomForest density catch-all with a shared-covariance Gaussian
classifier (per-group centroid + one pooled Ledoit-Wolf within-group covariance
-> Mahalanobis score), i.e. a regularized nearest-shrunk-centroid. Uniform class
priors mean template count cannot win; a single pooled covariance is stable with
the small per-group samples (11-383 templates) where per-group QDA was not.
Scores are temperature-scaled (T fit on the hold-out) so the C4 prior remains
effective and confidence calibrates. Edge bands 374 and 1034 nm (diagnosed as a
20x-concentrated systematic Gaia-SSO-vs-RELAB offset, see
diag/PASSBAND_FINDING.md) are dropped from the feature basis.
"""

import numpy as np
from numpy.linalg import inv
from sklearn.covariance import LedoitWolf
from sklearn.isotonic import IsotonicRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import balanced_accuracy_score

import c1_common as cc

DROP_BANDS = {374, 1034}
MATCH_BANDS = [b for b in cc.BANDS if b not in DROP_BANDS]
REF_NM = cc.REF_NM
FEAT_BANDS = [b for b in MATCH_BANDS if b != REF_NM]
_ALL = cc.BANDS
_FEAT_IDX = [_ALL.index(b) for b in FEAT_BANDS]
_REF_IDX = _ALL.index(REF_NM)
EPS = 1e-4


def features_from_16(V16):
    V = np.asarray(V16, float)
    if V.ndim == 1:
        V = V[None, :]
    ref = V[:, _REF_IDX]
    Vn = V / ref[:, None]
    return np.log(np.clip(Vn[:, _FEAT_IDX], EPS, None))


def _softmax(S, T):
    Z = S / T
    Z -= Z.max(1, keepdims=True)
    P = np.exp(Z)
    return P / P.sum(1, keepdims=True)


class GroupGaussianMatcher:
    def __init__(self):
        self.classes_ = None
        self.mu_ = None
        self.prec_ = None
        self.T_ = 1.0
        self.iso_ = None

    def fit(self, F, y):
        y = np.asarray(y)
        self.classes_ = np.array(sorted(set(y)))
        resid = np.empty_like(F)
        mus = []
        for c in self.classes_:
            m = F[y == c].mean(0)
            mus.append(m)
            resid[y == c] = F[y == c] - m
        self.mu_ = np.vstack(mus)
        self.prec_ = inv(LedoitWolf().fit(resid).covariance_)
        return self

    def _scores(self, F):
        F = np.atleast_2d(F)
        S = np.empty((F.shape[0], len(self.classes_)))
        for j in range(len(self.classes_)):
            d = F - self.mu_[j]
            S[:, j] = -0.5 * np.einsum('ij,jk,ik->i', d, self.prec_, d)
        return S

    def predict_proba(self, F):
        return _softmax(self._scores(F), self.T_)

    def fit_calibration(self, F, y, groups, n_splits=5):
        y = np.asarray(y)
        oofS = np.full((len(y), len(self.classes_)), -1e9)
        cidx = {c: i for i, c in enumerate(self.classes_)}
        for tr, te in GroupKFold(n_splits).split(F, y, groups):
            sub = GroupGaussianMatcher().fit(F[tr], y[tr])
            s = sub._scores(F[te])
            for jj, c in enumerate(sub.classes_):
                oofS[te, cidx[c]] = s[:, jj]
        yi = np.array([cidx[v] for v in y])
        grid = np.linspace(0.5, 30, 60)
        best, bT = 1e18, 1.0
        for T in grid:
            P = _softmax(oofS, T)
            nll = -np.log(np.clip(P[np.arange(len(yi)), yi], 1e-12, None)).mean()
            if nll < best:
                best, bT = nll, T
        self.T_ = float(bT)
        P = _softmax(oofS, self.T_)
        pred = self.classes_[P.argmax(1)]
        conf = P.max(1)
        corr = (pred == y).astype(float)
        self.iso_ = IsotonicRegression(out_of_bounds='clip').fit(conf, corr)
        bins = np.linspace(0, 1, 9)
        idx = np.clip(np.digitize(conf, bins) - 1, 0, len(bins) - 2)
        ece = 0.0
        for b in range(len(bins) - 1):
            mm = idx == b
            if mm.sum() >= 5:
                ece += mm.mean() * abs(conf[mm].mean() - corr[mm].mean())
        return balanced_accuracy_score(y, pred), float(ece), self.T_, conf, corr

    def calibrate(self, conf):
        return self.iso_.predict(np.asarray(conf)) if self.iso_ is not None else conf
