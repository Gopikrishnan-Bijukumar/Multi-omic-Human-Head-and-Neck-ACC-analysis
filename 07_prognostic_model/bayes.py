"""Evidence-informed Bayesian multi-level shrinkage model (torch / GPU).

A single vector of gene effects beta is shared by two cohorts:

    JSE  (n=20, binary Good/Poor)   ->  logistic likelihood,  logit = a_j + k_j * (Z_jse @ beta)
    DK   (n<=54, right-censored OS) ->  Cox partial likelihood, eta   =        Z_dk  @ beta

Multi-level shrinkage on beta (regularized horseshoe, Piironen & Vehtari 2017):

    beta_g ~ Normal( d_g * mu_dir * s_g ,  (tau * s_g * lam_tilde_g)^2 )
    lam_g  ~ HalfCauchy(1)                     gene level
    s_g    = evidence multiplier (FIXED, data-independent of DK)   evidence level
    tau    ~ HalfCauchy(tau0)                  global level

`s_g` is where the external evidence enters: genes carried by more modalities and
stronger discovery statistics get a wider prior, i.e. less shrinkage. `d_g` is the
evidence direction (+1 = up in Poor); `mu_dir >= 0` is a learned global pull toward
that direction, so the prior is informative in direction but the data may overrule it.

Inference: mean-field ADVI (reparameterised, non-centred). Point predictions use the
posterior mean of beta.

Nothing here reads or writes any project input file.
"""
import numpy as np
import torch
import torch.nn.functional as F

DEV = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def _sp(x):
    return F.softplus(x) + 1e-6


class JointShrinkageModel:
    def __init__(self, P, s_evidence, d_dir, tau0=0.05, slab_scale=1.0, slab_df=4.0,
                 use_jse=True, use_dk=True, use_dir_prior=True, seed=0):
        self.P = P
        self.s = torch.as_tensor(s_evidence, dtype=torch.float32, device=DEV)
        self.d = torch.as_tensor(d_dir, dtype=torch.float32, device=DEV)
        self.tau0 = tau0
        self.slab_scale = slab_scale
        self.slab_df = slab_df
        self.use_jse, self.use_dk, self.use_dir_prior = use_jse, use_dk, use_dir_prior
        g = torch.Generator(device="cpu").manual_seed(seed)

        def par(shape, val=0.0, scale=0.01):
            t = torch.full(shape, val) + scale * torch.randn(shape, generator=g)
            return t.to(DEV).requires_grad_(True)

        # variational params: q(z) = N(mu, softplus(rho)^2) for the non-centred params
        self.q_bz_m, self.q_bz_r = par((P,)), par((P,), -3.0)          # beta raw ~ N(0,1)
        self.q_ll_m, self.q_ll_r = par((P,), -1.0), par((P,), -3.0)     # log lambda
        self.q_lt_m, self.q_lt_r = par((1,), -2.0), par((1,), -3.0)     # log tau
        self.q_lc_m, self.q_lc_r = par((1,), 0.0), par((1,), -3.0)      # log c (slab)
        self.q_md_m, self.q_md_r = par((1,), -2.0), par((1,), -3.0)     # log mu_dir
        self.q_aj_m, self.q_aj_r = par((1,)), par((1,), -3.0)           # JSE intercept
        self.q_kj_m, self.q_kj_r = par((1,), 0.0), par((1,), -3.0)      # log JSE scale
        self.params = [getattr(self, n) for n in dir(self) if n.startswith("q_")]

    def _sample(self, S):
        out = {}
        for name in ["bz", "ll", "lt", "lc", "md", "aj", "kj"]:
            m = getattr(self, f"q_{name}_m")
            sd = _sp(getattr(self, f"q_{name}_r"))
            eps = torch.randn((S,) + m.shape, device=DEV)
            out[name] = m + sd * eps
            out[f"lq_{name}"] = (-0.5 * eps ** 2 - torch.log(sd) - 0.9189385).sum(-1)
        return out

    def _beta(self, z):
        lam = torch.exp(z["ll"])                      # (S,P)
        tau = torch.exp(z["lt"])                      # (S,1)
        c = self.slab_scale * torch.exp(z["lc"])      # (S,1)
        lam_t = c * lam / torch.sqrt(c ** 2 + (tau * lam) ** 2 + 1e-12)
        sd = tau * lam_t * self.s
        mu = self.d * torch.exp(z["md"]) * self.s if self.use_dir_prior else 0.0
        return mu + sd * z["bz"], lam, tau, c

    def _log_prior(self, z, beta, lam, tau, c):
        # non-centred: beta_raw ~ N(0,1)
        lp = (-0.5 * z["bz"] ** 2 - 0.9189385).sum(-1)
        # half-Cauchy on lam and tau, in log space (+ log|jacobian| = log x)
        lp = lp + (torch.log(2.0 / (np.pi * (1 + lam ** 2))) + z["ll"]).sum(-1)
        lp = lp + (torch.log(2.0 / (np.pi * self.tau0 * (1 + (tau / self.tau0) ** 2))) + z["lt"]).sum(-1)
        # inverse-gamma(df/2, df/2) slab on c^2 -> weakly informative; use log-normal proxy
        lp = lp + (-0.5 * (z["lc"] / 0.7) ** 2).sum(-1)
        lp = lp + (-0.5 * ((z["md"] + 2.0) / 1.5) ** 2).sum(-1)   # log mu_dir ~ N(-2, 1.5)
        lp = lp + (-0.5 * (z["aj"] / 2.0) ** 2).sum(-1)
        lp = lp + (-0.5 * (z["kj"] / 1.0) ** 2).sum(-1)
        return lp

    def _ll_jse(self, beta, z, Xj, yj):
        eta = z["aj"] + torch.exp(z["kj"]) * (beta @ Xj.T)         # (S,n)
        return -F.binary_cross_entropy_with_logits(eta, yj.expand_as(eta), reduction="none").sum(-1)

    def _ll_cox(self, beta, Xd, order, event):
        """Breslow partial likelihood. Rows of Xd already sorted by DESCENDING time."""
        eta = beta @ Xd.T                                          # (S,n)
        lse = torch.logcumsumexp(eta, dim=-1)                      # risk set = prefix (desc time)
        return ((eta - lse) * event).sum(-1)

    def fit(self, Xj=None, yj=None, Xd=None, event=None, steps=1500, S=8, lr=0.05, verbose=False):
        opt = torch.optim.Adam(self.params, lr=lr)
        sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, steps)
        nll_scale = 1.0
        for it in range(steps):
            opt.zero_grad()
            z = self._sample(S)
            beta, lam, tau, c = self._beta(z)
            lp = self._log_prior(z, beta, lam, tau, c)
            if self.use_jse and Xj is not None:
                lp = lp + self._ll_jse(beta, z, Xj, yj)
            if self.use_dk and Xd is not None:
                lp = lp + self._ll_cox(beta, Xd, None, event)
            lq = sum(z[f"lq_{n}"] for n in ["bz", "ll", "lt", "lc", "md", "aj", "kj"])
            loss = -(lp - lq).mean() * nll_scale
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.params, 10.0)
            opt.step(); sched.step()
            if verbose and it % 300 == 0:
                print(f"  step {it:5d} elbo {-loss.item():.2f}")
        return self

    @torch.no_grad()
    def beta_mean(self, S=400):
        z = self._sample(S)
        beta, *_ = self._beta(z)
        return beta.mean(0).cpu().numpy(), beta.std(0).cpu().numpy(), beta.cpu().numpy()


def prep_cox(X, time, event):
    """Sort by descending time for the Breslow prefix-risk-set trick."""
    o = np.argsort(-np.asarray(time, float), kind="stable")
    return (torch.as_tensor(np.asarray(X)[o], dtype=torch.float32, device=DEV),
            torch.as_tensor(np.asarray(event, float)[o], dtype=torch.float32, device=DEV), o)


def to_t(a):
    return torch.as_tensor(np.asarray(a), dtype=torch.float32, device=DEV)
