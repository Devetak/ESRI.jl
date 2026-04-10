# Mathematical Model

Let `W ∈ ℝ^{N×N}` be the firm-to-firm weight matrix where `W_{ij}` is the supply from firm `i` to firm `j`.

## Symbols

- `N`: number of firms.
- `W_{ij}`: supplier `i` to customer `j` weight.
- `r_i`: total output of firm `i` (`sum_j W_{ij}`).
- `u_i`: upstream health of firm `i`.
- `d_i`: downstream health of firm `i`.
- `\psi_i`: exogenous cap in `[0,1]` (shock scenario).
- `w_i`: final ESRI aggregation weights (default `w_i = r_i`).

## Upstream Impact Operator

Define firm output totals

```math
r_i = \sum_{j=1}^N W_{ij}.
```

The upstream impact matrix `U` is

```math
U_{ji} =
\begin{cases}
\dfrac{W_{ij}}{r_i}, & r_i > 0,\\
0, & r_i = 0.
\end{cases}
```

The upstream health update for scenario cap `\psi ∈ [0,1]^N` is

```math
u_i^{(t+1)} = \min\!\left(\sum_{j=1}^N U_{ji}u_j^{(t)},\ \psi_i\right).
```

## Downstream Impact Operators

For customer `j`, split suppliers into essential and non-essential groups:

```math
E_j(k)=\sum_{i:\, \text{industry}(i)=k,\ \text{essential}(i)} W_{ij},\quad
A_j=\sum_{i=1}^N W_{ij}.
```

Then

```math
D^{(e)}_{ij}=
\begin{cases}
\dfrac{W_{ij}}{E_j(\text{industry}(i))}, & \text{if }i\text{ is essential and }E_j(\cdot)>0,\\
0, & \text{otherwise},
\end{cases}
```

```math
D^{(n)}_{ij}=
\begin{cases}
\dfrac{W_{ij}}{A_j}, & \text{if }i\text{ is non-essential and }A_j>0,\\
0, & \text{otherwise}.
\end{cases}
```

## Coupled Downstream Iteration

Given downstream health `d^{(t)}`, define industry-normalization

```math
\sigma_i^{(t)}=
\begin{cases}
0, & r_i=0,\\
1, & \sum_{m:\,\text{industry}(m)=\text{industry}(i)} r_m d_m^{(t)}=0,\\
\min\!\left(
\dfrac{r_i}{
\sum_{m:\,\text{industry}(m)=\text{industry}(i)} r_m d_m^{(t)}
},1\right), & \text{otherwise}.
\end{cases}
```

Industry shortages and non-essential shortage for customer `j`:

```math
S_{jk}^{(t)}=\sum_{i=1}^N \sigma_i^{(t)}(1-d_i^{(t)})D^{(e)}_{ij}\,\mathbf{1}_{\text{industry}(i)=k},
```

```math
N_j^{(t)}=\sum_{i=1}^N \sigma_i^{(t)}(1-d_i^{(t)})D^{(n)}_{ij}.
```

Update:

```math
d_j^{(t+1)}=\min\!\bigl(1-\max_k S_{jk}^{(t)},\ 1-N_j^{(t)},\ \psi_j\bigr).
```

The implementation performs a joint fixed-point iteration over `(u, d)` and stops when the infinity-norm change across consecutive iterates is below `tol` (or `maxiter` is reached).

## ESRI Reduction

After convergence, define final firm health `f_i`:

- `combine = :min`: `f_i = min(u_i, d_i)`
- `combine = :upstream`: `f_i = u_i`
- `combine = :downstream`: `f_i = d_i`

For final weights `w_i` (default `w_i = r_i`):

```math
\mathrm{ESRI} = \frac{\sum_{i=1}^N w_i (1-f_i)}{\sum_{i=1}^N r_i}.
```

The package computes this value for one firm shock, a subset of firms, or custom shock vectors.

## Equation-to-Code Mapping

- `U` is built by `create_upstream_impact_matrix`.
- `D^{(e)}, D^{(n)}` are built by `compute_downstream_impact_matrices`.
- `\sigma` and downstream shortages are updated by `compute_sigmas!`, `_accumulate_downstream_components!`, and `downstream_step!`.
- upstream updates apply `u^{(t+1)} = min(Uu^{(t)}, \psi)` via `upstream_step!`.
- final scalar ESRI reduction uses `_reduce_esri` with `combine` selection.
