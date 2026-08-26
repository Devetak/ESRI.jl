# Mathematical Model

Let ``W \in \mathbb{R}_{\ge 0}^{N \times N}`` be the firm-to-firm weight
matrix. Entry ``W_{ij}`` is the supply from firm ``i`` to firm ``j``. Let
``g(i)`` be the industry of firm ``i``. The input-classification matrix ``C`` has one row per
supplying industry and one column per customer industry:

```math
C_{k\ell} \in \{0,1,2\}.
```

Here ``2`` marks industry ``k`` as essential to industry ``\ell``. ``1`` marks
a linear input. ``0`` keeps the input in ``W`` and the upstream operator and
contributes to the linear-input denominator. The legacy Boolean vector sets each
row of ``C`` to ``2`` or ``1``.

The one-argument `IndustryInfo(industry_of_firm)` constructor uses ``C = 1``
for every supplier-customer pair, giving the purely linear baseline. The
bundled IHS matrix is opt-in through `ihs_input_classification()`.

The package also uses a capacity-cap vector ``\psi \in [0,1]^N``.

The entry ``\psi_i`` is the exogenous capacity limit for firm ``i``. The firm
can operate at up to a ``\psi_i`` fraction of its normal capacity. Use this vector to
describe a closed firm, a partly constrained sector, or a wider event such as
an energy shortage, sanctions, or a port disruption. ESRI shows how the initial
shock moves through supply chains and affects the wider economy.

In simple language:

- ``\psi_i = 1`` makes full capacity available to firm ``i``.
- ``\psi_i = 0`` closes firm ``i``.
- ``0 < \psi_i < 1`` sets a partial capacity limit for firm ``i``.

A vector of ones represents an unshocked economy.

The ESRI score in Diem et al. is the share of total steady-state output lost
after the full closure of one firm. The ``\psi`` formulation also supports
simultaneous shocks and partial capacity limits.

Define firm output and input totals by

```math
r_i = \sum_{j=1}^N W_{ij},
\qquad
c_j = \sum_{i=1}^N W_{ij}.
```

Let the exogenous scenario be the vector ``\psi \in [0,1]^N``. The package
solves for upstream health ``u \in [0,1]^N`` and downstream health
``d \in [0,1]^N``.

## Upstream operator

The upstream operator propagates customer losses to suppliers.

The upstream operator ``U`` is

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

The downstream operators propagate input shortages to customers.

For each customer ``j``, define the total input weight from supplier industry
``k`` by

```math
E_{jk} = \sum_{i=1}^N W_{ij} \, \mathbf{1}_{g(i)=k}.
```

Define the essential downstream operator ``D^{(e)}`` and the linear-input
downstream operator ``D^{(n)}`` by

```math
D^{(e)}_{ij} =
\begin{cases}
\dfrac{W_{ij}}{E_{j,g(i)}}, & C_{g(i),g(j)} = 2 \text{ and } E_{j,g(i)} > 0, \\
0, & \text{otherwise},
\end{cases}
```

```math
D^{(n)}_{ij} =
\begin{cases}
\dfrac{W_{ij}}{c_j}, & C_{g(i),g(j)} = 1 \text{ and } c_j > 0, \\
0, & \text{otherwise}.
\end{cases}
```

Class-``0`` inputs contribute to ``W`` and the upstream operator. The downstream
operators use class-``2`` and class-``1`` inputs. Class-``0`` weights remain part
of ``c_j`` when class-``1`` inputs are normalized.

## Supplier rationing factor

Given ``d^{(t)}``, define the current industry mass

```math
M_k^{(t)} = \sum_{m:g(m)=k} r_m d_m^{(t)}.
```

The package defines the supplier rationing factor

```math
\sigma_i^{(t)} =
\begin{cases}
0, & r_i = 0, \\
1, & r_i > 0 \text{ and } M_{g(i)}^{(t)} = 0, \\
\min\!\left(\dfrac{r_i}{M_{g(i)}^{(t)}}, 1\right), & r_i > 0 \text{ and } M_{g(i)}^{(t)} > 0.
\end{cases}
```

## Downstream update

For each customer ``j`` and industry ``k``, define the essential shortage term

```math
S_{jk}^{(t)} =
\sum_{i=1}^N
\sigma_i^{(t)}
\left(1 - d_i^{(t)}\right)
D^{(e)}_{ij}
\mathbf{1}_{g(i)=k}.
```

Define the linear-input shortage term

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

It iterates the upstream and downstream recursions in lockstep until

```math
\max\!\left(
\|u^{(t+1)} - u^{(t)}\|_{\infty},
\|d^{(t+1)} - d^{(t)}\|_{\infty}
\right)
< \mathrm{tol},
```

or until `maxiter` iterations have been reached.

Define the final firm health ``f_i`` by

```math
f_i =
\begin{cases}
\min(u_i, d_i), & \texttt{combine} = :\mathrm{min}, \\
u_i, & \texttt{combine} = :\mathrm{upstream}, \\
d_i, & \texttt{combine} = :\mathrm{downstream}.
\end{cases}
```

Let ``w_i`` be the final weights. The package default is ``w_i = r_i``. The
reported scalar is

```math
\mathrm{ESRI} =
\frac{\sum_{i=1}^N w_i (1 - f_i)}{\sum_{i=1}^N w_i}.
```

For ``\sum_i w_i = 0``, the package returns ``0``.

With the default ``w_i = r_i``, this is a share of total output. With custom
`final_weights`, it is a weighted average loss under those weights.

## Relation to the paper

The package follows the ESRI setup from Diem et al., Scientific Reports 12,
7719 (2022). The package supports general shock scenarios through ``\psi``.

For the submitted runtime comparison, see [Performance](performance.md).

## References

Christian Diem, András Borsos, Tobias Reisch, János Kertész, Stefan Thurner.
Quantifying firm-level economic systemic risk from nation-wide supply networks.
Scientific Reports 12, 7719 (2022).
[doi:10.1038/s41598-022-11522-z](https://doi.org/10.1038/s41598-022-11522-z).
