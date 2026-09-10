"""Evidence-informed Bayesian MULTI-LEVEL shrinkage model (torch / GPU, ADVI).

Three levels of shrinkage on the gene effect vector beta (log-hazard per rank-normal SD):

  level 3 (source) : each module SOURCE (curated / wgcna / single_cell / evidence)
                     has its own scale  omega_s ~ HalfCauchy(0.1)
  level 2 (module) : theta_k ~ Normal(0, (omega_{s(k)} * lam_k)^2),  lam_k ~ HalfCauchy(1)
  level 1 (gene)   : delta_g ~ Normal(d_g * mu_dir * s_g, (tau * s_g * lamt_g)^2)   [reg. horseshoe]

      beta = M @ theta + delta

M is the (P x K) module design matrix; s_g the evidence multiplier and d_g the evidence
direction, both from the discovery cohort only. Genes therefore borrow strength through
the biological programme they belong to, which is what makes 38 events informative.

Likelihoods (shared beta):
    DK   : Cox partial likelihood (Breslow)
    JSE  : logistic,  logit = a + exp(k) * (Z beta)

Inference: reparameterised mean-field ADVI. Predictions use the posterior mean of beta.
"""
import numpy as np
import torch
import torch.nn.functional as F

DEV = torch.device("cuda" if torch.cuda.is_available() else "cpu")
LOG2PI2 = 0.9189385332046727


def _sp(x):
    return F.softplus(x) + 1e-6


class MultiLevelModel:
    def __init__(self, M, src_ids, n_src, s_evidence, d_dir, tau0=0.02, omega0=0.1,
                 slab=1.0, use_jse=True, use_dk=True, use_gene_level=True,
                 use_module_level=True, use_dir_prior=True, seed=0):
        self.M = torch.as_tensor(np.ascontiguousarray(M), dtype=torch.float32, device=DEV)
        self.P, self.K = self.M.shape
        self.src = torch.as_tensor(np.asarray(src_ids), dtype=torch.long, device=DEV)
        self.n_src = n_src
        self.s = torch.as_tensor(np.ascontiguousarray(s_evidence), dtype=torch.float32, device=DEV)
        self.d = torch.as_tensor(np.ascontiguousarray(d_dir), dtype=torch.float32, device=DEV)
        self.tau0, self.omega0, self.slab = tau0, omega0, slab
        self.use_jse, self.use_dk = use_jse, use_dk
        self.use_gene, self.use_mod, self.use_dir = use_gene_level, use_module_level, use_dir_prior
        g = torch.Generator().manual_seed(seed)

        def par(shape, val=0.0, sc=0.01):
            return (torch.full(shape, val) + sc * torch.randn(shape, generator=g)).to(DEV).requires_grad_(True)

        self.names = []
        def add(n, shape, val=0.0, rval=-3.0):
            setattr(self, f"q_{n}_m", par(shape, val))
            setattr(self, f"q_{n}_r", par(shape, rval))
            self.names.append(n)

        add("tz", (self.K,))                 # theta raw
        add("lk", (self.K,), -1.0)           # log lambda_k
        add("lo", (n_src,), np.log(omega0))  # log omega_source
        add("dz", (self.P,))                 # delta raw
        add("lg", (self.P,), -1.0)           # log lambda_g
        add("lt", (1,), np.log(tau0))        # log tau
        add("lc", (1,), 0.0)                 # log slab
        add("md", (1,), -2.0)                # log mu_dir
        add("aj", (1,))                      # JSE intercept
        add("kj", (1,), 0.0)                 # log JSE scale
        self.params = [getattr(self, f"q_{n}_{t}") for n in self.names for t in ("m", "r")]

    def _sample(self, S):
        z, lq = {}, 0.0
        for n in self.names:
            m = getattr(self, f"q_{n}_m"); sd = _sp(getattr(self, f"q_{n}_r"))
            eps = torch.randn((S,) + m.shape, device=DEV)
            z[n] = m + sd * eps
            lq = lq + (-0.5 * eps ** 2 - torch.log(sd) - LOG2PI2).sum(-1)
        return z, lq

    def _beta(self, z):
        lam_k = torch.exp(z["lk"])
        omega = torch.exp(z["lo"])[:, self.src] if False else torch.exp(z["lo"]).index_select(-1, self.src)
        theta = z["tz"] * lam_k * omega if self.use_mod else torch.zeros_like(z["tz"])
        lam_g = torch.exp(z["lg"]); tau = torch.exp(z["lt"]); c = self.slab * torch.exp(z["lc"])
        lamt = c * lam_g / torch.sqrt(c ** 2 + (tau * lam_g) ** 2 + 1e-12)
        mu = self.d * torch.exp(z["md"]) * self.s if self.use_dir else 0.0
        delta = mu + tau * lamt * self.s * z["dz"] if self.use_gene else torch.zeros_like(z["dz"])
        beta = delta + (theta @ self.M.T if self.use_mod else 0.0)
        return beta, theta, lam_k, lam_g, tau, omega

    def _log_prior(self, z, lam_k, lam_g, tau, omega):
        lp = (-0.5 * z["tz"] ** 2 - LOG2PI2).sum(-1) + (-0.5 * z["dz"] ** 2 - LOG2PI2).sum(-1)
        lp = lp + (torch.log(2.0 / (np.pi * (1 + lam_k ** 2))) + z["lk"]).sum(-1)
        lp = lp + (torch.log(2.0 / (np.pi * (1 + lam_g ** 2))) + z["lg"]).sum(-1)
        om = torch.exp(z["lo"])
        lp = lp + (torch.log(2.0 / (np.pi * self.omega0 * (1 + (om / self.omega0) ** 2))) + z["lo"]).sum(-1)
        lp = lp + (torch.log(2.0 / (np.pi * self.tau0 * (1 + (tau / self.tau0) ** 2))) + z["lt"]).sum(-1)
        lp = lp + (-0.5 * (z["lc"] / 0.7) ** 2).sum(-1)
        lp = lp + (-0.5 * ((z["md"] + 2.0) / 1.5) ** 2).sum(-1)
        lp = lp + (-0.5 * (z["aj"] / 2.0) ** 2).sum(-1) + (-0.5 * (z["kj"] / 1.0) ** 2).sum(-1)
        return lp

    def fit(self, Xj=None, yj=None, Xd=None, event=None, steps=1500, S=8, lr=0.05, verbose=False,
            trace=False):
        # `trace` is purely additive instrumentation (stage 33 inference diagnostics):
        # it records the per-step ELBO and consumes no random numbers, so fits with
        # trace=True are bit-identical to the runs that produced the locked results.
        self.elbo_trace = [] if trace else None
        opt = torch.optim.Adam(self.params, lr=lr)
        sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, steps)
        for it in range(steps):
            opt.zero_grad()
            z, lq = self._sample(S)
            beta, theta, lam_k, lam_g, tau, omega = self._beta(z)
            lp = self._log_prior(z, lam_k, lam_g, tau, omega)
            if self.use_dk and Xd is not None:
                eta = beta @ Xd.T
                lp = lp + ((eta - torch.logcumsumexp(eta, -1)) * event).sum(-1)
            if self.use_jse and Xj is not None:
                ej = z["aj"] + torch.exp(z["kj"]) * (beta @ Xj.T)
                lp = lp - F.binary_cross_entropy_with_logits(ej, yj.expand_as(ej), reduction="none").sum(-1)
            loss = -(lp - lq).mean()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.params, 10.0)
            opt.step(); sch.step()
            if trace:
                self.elbo_trace.append(-loss.item())
            if verbose and it % 400 == 0:
                print(f"   {it:5d} elbo {-loss.item():.1f}")
        return self

    @torch.no_grad()
    def posterior(self, S=500):
        z, _ = self._sample(S)
        beta, theta, *_ = self._beta(z)
        return dict(beta=beta.mean(0).cpu().numpy(), beta_sd=beta.std(0).cpu().numpy(),
                    theta=theta.mean(0).cpu().numpy(), theta_sd=theta.std(0).cpu().numpy(),
                    theta_pdir=(theta > 0).float().mean(0).cpu().numpy())


def prep_cox(X, time, event):
    o = np.argsort(-np.asarray(time, float), kind="stable")
    return (torch.as_tensor(np.ascontiguousarray(np.asarray(X)[o]), dtype=torch.float32, device=DEV),
            torch.as_tensor(np.ascontiguousarray(np.asarray(event, float)[o]), dtype=torch.float32, device=DEV))


def to_t(a):
    return torch.as_tensor(np.ascontiguousarray(np.asarray(a)), dtype=torch.float32, device=DEV)
