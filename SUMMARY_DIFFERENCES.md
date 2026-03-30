# Summary differences: ESRI.jl vs Economic-Systemic-Risk

Notation: firms \(i,j \in \{1,\ldots,n\}\); industries \(k\); directed input–output weights \(W_{ij}\) (supplier \(i\), customer \(j\)); row sums \(s^{\mathrm{out}}_i = \sum_j W_{ij}\), column sums \(s^{\mathrm{in}}_j = \sum_i W_{ij}\). For a shock scenario, let \(h^{u}_i,h^{d}_i \in [0,1]\) be upstream/downstream “health” and \(\psi_i\) the exogenous cap from the shock.

---

## 1) ESRI scalar (both implementations)

After convergence, define production factor \(f_i = \min(h^{u}_i, h^{d}_i)\) and disruption \(1 - f_i\).

**Both** aggregate with supplier weights \(s^{\mathrm{out}}_i\) and normalize by total output \(Y = \sum_i s^{\mathrm{out}}_i\):

\[
\mathrm{ESRI} = \frac{1}{Y} \sum_{i=1}^{n} s^{\mathrm{out}}_i \, \bigl(1 - f_i\bigr).
\]

(Small differences can still arise from different \(h^{u}, h^{d}\) dynamics, not from this aggregation formula.)

---

## 2) Economic-Systemic-Risk — downstream block

### 2.1 Market shares (sector relative to \(h^d\))

Let \(V_k(h^d) = \sum_{\ell : n_\ell = k} s^{\mathrm{out}}_\ell \, h^{d}_\ell\). For firm \(i\) with sector \(n_i\) and \(s^{\mathrm{out}}_i > 0\):

\[
m_i = \min\!\left(1,\; \frac{s^{\mathrm{out}}_i \, h^{d}_i}{V_{n_i}(h^d)}\right)
\quad (\text{and } m_i = 0 \text{ if } s^{\mathrm{out}}_i = 0).
\]

### 2.2 Downstream supply indices

For customer \(j\), build essential contributions by **supplier industry** \(k\) (NACE), and one **non-essential aggregate** \(D^{\mathrm{ne}}\):

\[
D^{\mathrm{ess}}_{jk} = \sum_{\substack{i \in \mathcal{N}^+(j)\\ i \text{ type essential},\, n_i=k}}
 m_i \, [\Lambda_{d1}]_{ij} \, (1 - h^{d}_i),
\]

\[
D^{\mathrm{ne}}_j = \sum_{\substack{i \in \mathcal{N}^+(j)\\ i \text{ type non-essential}}}
 m_i \, [\Lambda_{d2}]_{ij} \, (1 - h^{d}_i).
\]

Then

\[
h^{\mathrm{ess}}_j = 1 - \max_k D^{\mathrm{ess}}_{jk},
\qquad
h^{\mathrm{ne}}_j = 1 - D^{\mathrm{ne}}_j,
\]

\[
h^{d,\mathrm{new}}_j = \min\!\bigl(h^{\mathrm{ess}}_j,\, h^{\mathrm{ne}}_j,\, \psi_j\bigr).
\]

**Economic meaning:** essential inputs tie across industries like a **multi-leontief bottleneck** (\(\max\) over shortages); non-essential inputs enter as a **single pooled** constraint; the firm’s downstream health is the **minimum** of those channels and the shock cap.

### 2.3 Upstream (reference)

\[
h^{u,\mathrm{new}}_j = \min\!\left(
\sum_{i \in \mathcal{N}^-(j)} [\Lambda_u]_{ij}\, h^{u}_i,\; \psi_j
\right)
\]
(with a boundary case: no customers \(\Rightarrow h^u_j = 1\)).

Matrices \(\Lambda_u,\Lambda_{d1},\Lambda_{d2}\) come from `buildArrays` (edge lists and types).

---

## 3) ESRI.jl (current `src/`) — downstream block

### 3.1 Downstream impact matrix \(D\) (`impact.jl`)

For column \(j\) (customer \(j\)), let \(S_j = s^{\mathrm{in}}_j = \sum_i W_{ij}\).  
Let \(\mathcal{E}\) be firms in industries marked essential; \(g(i)\) industry of firm \(i\).

Define partial sums over **essential suppliers** only (and only when **customer \(j\)** is essential):

\[
P_{jg} =
\begin{cases}
\displaystyle\sum_{\substack{i:\, W_{ij}>0,\, i\in\mathcal{E}}} W_{ij}\,\mathbf{1}\{g(i)=g\} & j \in \mathcal{E} \\[0.5em]
0 & j \notin \mathcal{E}
\end{cases}
\]

Then entries \(D_{ij}\) of the downstream impact matrix satisfy:

\[
D_{ij} =
\begin{cases}
W_{ij} / P_{j\,g(i)} & j \in \mathcal{E},\; i \in \mathcal{E},\; P_{j\,g(i)}>0 \\[0.5em]
W_{ij} / S_j & \text{otherwise.}
\end{cases}
\]

**Difference from reference:** normalization uses “customer is essential” as a gate for sector-wise essential sums; reference uses **edge type** (essential vs non-essential links) and builds separate \(\Lambda_{d1},\Lambda_{d2}\).

### 3.2 Iteration: \(\sigma\) and product matrix (`downstream.jl`)

Let \(h^d\) be current downstream health. Sector masses:

\[
T_k(h^d) = \sum_{\ell:\, g(\ell)=k} s^{\mathrm{out}}_\ell \, h^{d}_\ell.
\]

Then

\[
\sigma_i = \min\!\left(1,\; \frac{s^{\mathrm{out}}_i}{T_{g(i)}(h^d)}\right)
\quad (\text{or } 1 \text{ if } T_{g(i)}=0).
\]

Initialize \(P \leftarrow \mathbf{1}\,\mathbf{1}^\top\) (matrix of ones). For each firm \(c\), with \(a_c = \sigma_c (1-h^d_c)\) and industry \(g(c)\):

\[
P_{r\,g(c)} \leftarrow P_{r\,g(c)} - a_c \, D_{cr}
\quad \forall r.
\]

So updates are **grouped by supplier firm’s own industry column** \(g(c)\), not by “essential sector columns + one non-essential bucket” like the reference.

### 3.3 Downstream step: essential vs non-essential **firm**

Let \(P\) have \(K\) columns (one per industry). Let \(\eta_i = \mathbf{1}\{i \in \mathcal{E}\}\) (firm essentiality from its industry).

\[
h^{d,\mathrm{new}}_i =
\begin{cases}
\displaystyle\min_{k=1,\ldots,K} P_{ik} & \eta_i = 1 \\[0.75em]
P_{iK} & \eta_i = 0
\end{cases}
\]

Then the shocked firm index is set back to \(0\) inside `downstream_shock!`.

**Difference from reference:**

| Feature | Economic-Systemic-Risk | ESRI.jl (`downstream_step!`) |
|--------|-------------------------|------------------------------|
| Who gets \(\min(\cdot,\cdot,\psi)\)? | **Every** firm: \(\min(h^{\mathrm{ess}}, h^{\mathrm{ne}}, \psi)\). | **Essential firms:** min over **all** industry columns. **Non-essential firms:** only **last** column \(P_{iK}\). |
| Non-essential bottleneck | Explicit \(h^{\mathrm{ne}} = 1 - D^{\mathrm{ne}}\). | Last column of \(P\); not the same object as reference \(h^{\mathrm{ne}}\). |

So the production-feasibility / bottleneck logic is **not** the same map from \((W,\text{types},h)\) to \(h^{d,\mathrm{new}}\).

---

## 4) Shock and outer loop

**Reference** (newer branch): scenarios from `psi_mat`; column \(s\) sets \(\psi\) on a subset of firms; fixed-point until \(\| \Delta \|_\infty < 10^{-2}\) or \(t < t_{\max}\).

**ESRI.jl:** one scenario per shocked firm; \(\psi\) implied by upstream/downstream routines; same broad tolerance idea in `downstream_shock!` / `upstream_shock!` but **not** identical to the reference loop structure.

---

## 5) Can reference CLI flags make old `src/` match?

On the newer branch, **`--psi_mat`**, **`--tmax`**, **`--timeseries`** change *which scenarios* are run and *how long* the fixed point iterates — not the **functional form** of \(h^{\mathrm{ess}}, h^{\mathrm{ne}}, h^{d,\mathrm{new}}\).

So **no**: there is no flag that turns reference downstream into “essential firm = min over industries, non-essential firm = last column only” without changing `functions.jl` / `buildArrays`.

---

## 6) Bottom line

Disagreement between **reverted** `src/downstream.jl` + `src/impact.jl` and **Economic-Systemic-Risk** is expected because:

1. **Downstream impact** \(D\) vs \(\Lambda_{d1},\Lambda_{d2}\) use different conditioning (customer essential vs link type; single vs split channels).
2. **Downstream update** applies a **firm-type split** (min all columns vs last column only), whereas reference always uses \(\min\) over **essential bottleneck, non-essential aggregate, and \(\psi\)**.

Aligning outputs requires aligning these maps (either by changing ESRI.jl or by changing the reference implementation), not by tweaking `tmax` / `psi` alone.
