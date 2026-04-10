# Mathematical Model

Let `W \in \mathbb{R}_{\ge 0}^{N \times N}` be the firm-to-firm weight matrix with entry `W_{ij}` equal to supply from firm `i` to firm `j`. Let `g(i)` be the industry of firm `i`. Let `e_k \in \{0,1\}` be the essentiality flag of industry `k`. In the package input, firm `i` is essential exactly when `e_{g(i)} = 1`.

Define firm output and input totals by

```math
r_i = \sum_{j=1}^N W_{ij},
\qquad
c_j = \sum_{i=1}^N W_{ij}.
```

Let the exogenous scenario be the vector `\psi \in [0,1]^N`. The package solves for upstream health `u \in [0,1]^N` and downstream health `d \in [0,1]^N`.

## Upstream operator

The upstream operator `U` is

```math
U_{ji} =
\begin{cases}
\dfrac{W_{ij}}{r_i}, & r_i > 0, \\
0, & r_i = 0.
\end{cases}
```

The upstream update is

```math
u_i^{(t+1)} =
\begin{cases}
\min(1,\psi_i), & r_i = 0, \\
\min\!\left(\sum_{j=1}^N U_{ji} u_j^{(t)}, \psi_i\right), & r_i > 0.
\end{cases}
```

## Downstream operators

For each customer `j`, define the total weight of essential suppliers from industry `k` by

```math
E_{jk} = \sum_{i=1}^N W_{ij} \, \mathbf{1}_{e_{g(i)} = 1} \, \mathbf{1}_{g(i)=k}.
```

Define the essential downstream operator `D^{(e)}` and the non-essential downstream operator `D^{(n)}` by

```math
D^{(e)}_{ij} =
\begin{cases}
\dfrac{W_{ij}}{E_{j,g(i)}}, & e_{g(i)} = 1 \text{ and } E_{j,g(i)} > 0, \\
0, & \text{otherwise},
\end{cases}
```

```math
D^{(n)}_{ij} =
\begin{cases}
\dfrac{W_{ij}}{c_j}, & e_{g(i)} = 0 \text{ and } c_j > 0, \\
0, & \text{otherwise}.
\end{cases}
```

## Industry scaling

Given `d^{(t)}`, define the current industry mass

```math
M_k^{(t)} = \sum_{m:g(m)=k} r_m d_m^{(t)}.
```

The package then defines

```math
\sigma_i^{(t)} =
\begin{cases}
0, & r_i = 0, \\
1, & r_i > 0 \text{ and } M_{g(i)}^{(t)} = 0, \\
\min\!\left(\dfrac{r_i}{M_{g(i)}^{(t)}}, 1\right), & r_i > 0 \text{ and } M_{g(i)}^{(t)} > 0.
\end{cases}
```

## Downstream update

For each customer `j` and industry `k`, define the essential shortage term

```math
S_{jk}^{(t)} =
\sum_{i=1}^N
\sigma_i^{(t)}
\left(1 - d_i^{(t)}\right)
D^{(e)}_{ij}
\mathbf{1}_{g(i)=k}.
```

Define the non-essential shortage term

```math
N_j^{(t)} =
\sum_{i=1}^N
\sigma_i^{(t)}
\left(1 - d_i^{(t)}\right)
D^{(n)}_{ij}.
```

Set

```math
h_{j,\mathrm{ess}}^{(t)} = 1 - \max_k S_{jk}^{(t)},
\qquad
h_{j,\mathrm{non}}^{(t)} = 1 - N_j^{(t)}.
```

The downstream update is

```math
d_j^{(t+1)} = \min\!\left(h_{j,\mathrm{ess}}^{(t)}, h_{j,\mathrm{non}}^{(t)}, \psi_j\right).
```

## Fixed point and reduction

The package starts from

```math
u^{(0)} = \mathbf{1},
\qquad
d^{(0)} = \mathbf{1}.
```

It iterates the coupled upstream and downstream map until

```math
\max\!\left(
\|u^{(t+1)} - u^{(t)}\|_{\infty},
\|d^{(t+1)} - d^{(t)}\|_{\infty}
\right)
< \mathrm{tol},
```

or until `maxiter` iterations have been reached.

Define the final firm health `f_i` by

```math
f_i =
\begin{cases}
\min(u_i, d_i), & \texttt{combine} = :\mathrm{min}, \\
u_i, & \texttt{combine} = :\mathrm{upstream}, \\
d_i, & \texttt{combine} = :\mathrm{downstream}.
\end{cases}
```

Let `w_i` be the final weights. The package default is `w_i = r_i`. The reported scalar is

```math
\mathrm{ESRI} =
\frac{\sum_{i=1}^N w_i (1 - f_i)}{\sum_{i=1}^N r_i}.
```

If `\sum_i r_i = 0`, the package returns the unnormalized numerator.

## Relation to the paper

The package follows the same ESRI setup used in Diem et al., Scientific Reports 12, 6214 (2022), and the accompanying supplementary information in `paper/`. The package input contract is narrower. It takes firm industry labels and one Boolean essentiality flag per industry. The original reference code in `Economic-Systemic-Risk/` classifies edge types from a sector-pair table and constructs scenarios through `psi_mat`. The formulas on this page state the current package implementation.

## References

Christian Diem, Andras Borsos, Tobias Reisch, Janos Kertesz, Stefan Thurner. Quantifying firm-level economic systemic risk from nation-wide supply networks. Scientific Reports 12, 6214, 2022. DOI `10.1038/s41598-022-11522-z`.

Supporting information is bundled in `paper/41598_2022_11522_MOESM1_ESM.pdf`.
