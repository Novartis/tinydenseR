# Graph-Diffused, Density Contrast-Aligned PLS Decomposition (plsD)

Decomposes an expression interaction matrix M.local via NIPALS PLS1
against the density-contrast vector Y. plsD generates candidate features
that drive the density contrast — including population markers (DA),
differentially expressed genes (DE), and their mixture — by maximizing
covariance between graph-smoothed expression and Y.

## Usage

``` r
get.plsD(x, ...)

# S3 method for class 'TDRObj'
get.plsD(
  x,
  .coef.col,
  .model.name = "default",
  .ncomp = NULL,
  .min.prop = 0.005,
  .store.M = FALSE,
  .degree.reg = FALSE,
  .tau.mult = 1,
  .lazy.alpha = 1,
  .YX.interaction = TRUE,
  .loading.method = c("pearson", "ols", "spearman"),
  .residualize = FALSE,
  .verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object after
  [`get.lm()`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md).

- ...:

  Additional arguments passed to methods.

- .coef.col:

  Character: column name in model coefficients to use as density
  contrast Y. Must be a valid column in
  `.tdr.obj$map$lm[[.model.name]]$fit$coefficients`.

- .model.name:

  Character: name of the fitted model to use (default "default").

- .ncomp:

  Integer: number of PLS components. Defaults to
  `ncol(.tdr.obj$pca$embed)`, matching the number of PCs from
  `get.landmarks`.

- .min.prop:

  Numeric: for RNA, minimum proportion of landmarks where a gene must be
  detected (\>0) to be included. Default 0.005 (0.5 percent).

- .store.M:

  Logical: if TRUE, store M.local in output. Default FALSE (saves
  memory; can be large for RNA).

- .degree.reg:

  Logical: if TRUE, apply degree regularization by adding tau \* I
  (self-loops) to SNN before row-normalizing. Default FALSE.

- .tau.mult:

  Numeric: multiplier for tau when .degree.reg = TRUE. tau = .tau.mult
  \* mean(degree). Default 1.

- .lazy.alpha:

  Numeric in (0, 1\]: mixing parameter for lazy random walk. P_lazy =
  (1 - alpha) \* I + alpha \* P. Default 1 (standard walk). Values \< 1
  reduce the spatial propagation range of expression scores on the
  graph, making the decomposition more local. This is particularly
  useful when the density contrast is dominated by a single biologically
  extreme population: reducing `.lazy.alpha` limits how far the
  structural score counterweight (see Details) propagates to distant
  manifold regions. Try 0.5 as a first step when strong unexpected
  balancing is observed. Values \< 1 also incidentally damp oscillations
  on sparse/irregular graphs.

- .YX.interaction:

  Logical: if TRUE (default), construct M.local = P diag(Y) Xc
  (Y-weighted interaction). Loadings capture candidate features driving
  the density contrast through both differential abundance (population
  markers) and differential expression. If FALSE, construct M.local = P
  Xc (graph-smoothed expression only; Y appears only on the response
  side). The FALSE mode is less comprehensive — it misses DA markers
  unless their expression independently covaries with Y. Note: Y is
  centered before forming diag(Y), so the interaction weights landmarks
  by deviation from the mean coefficient, not from zero. Landmarks with
  below-average (but still positive) raw coefficients receive a negative
  weight, contributing to the structural opposite pole in the
  decomposition.

- .loading.method:

  Character: method for computing feature loadings. `"pearson"` computes
  Pearson correlation between component scores and centered expression
  (scale-invariant, fast sparse-BLAS path). `"ols"` computes OLS
  regression coefficients of centered Xc on component scores (same
  sparse-BLAS path; preserves magnitude information but is
  scale-dependent, so high-variance features rank higher). `"spearman"`
  computes Spearman rank correlation via a sparse-aware implementation
  (robust to outliers and nonlinearity; slower, O(nnz \* log(m_g) \*
  K)).

- .residualize:

  Logical: if TRUE, project out nuisance covariates from the expression
  matrix before decomposition. Default FALSE. Nuisance columns are
  identified automatically: when contrasts are present, design columns
  that do not participate in any contrast are nuisance; when no
  contrasts are used, all design columns except `.coef.col` are
  nuisance. If
  [`get.lm()`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md)
  was called with `.block`, the blocking variable is additionally
  included in the nuisance design (dummy-coded). The design is expanded
  from sample level to landmark level via `@config$key`. Sparsity of the
  expression matrix is preserved via implicit operators (the
  residualized matrix is never materialized). Not supported with
  `.loading.method = "spearman"`.

- .verbose:

  Logical: print progress messages? Default TRUE.

## Value

The modified `.tdr.obj` with results stored in
`.tdr.obj$plsD[[.coef.col]]`:

- scores:

  Matrix (landmarks x K): PLS scores (oriented so positive = aligned
  with Y)

- feature.weights:

  Matrix (genes x K): PLS feature weight vectors w (unit-norm)

- x.loadings:

  Matrix (genes x K): PLS deflation loadings p

- y.loadings:

  Numeric vector (K): PLS Y-loadings q (scalar per component)

- raw.loadings:

  Matrix (genes x K): raw feature loadings per component (Pearson r, OLS
  regression, or Spearman rank correlation, depending on
  `.loading.method`); positive = upregulated with Y, negative =
  downregulated

- loadings:

  Matrix (genes x K): concordance-filtered feature loadings, defined as
  `raw.loadings * concordance.weights` for Pearson/OLS. For
  `.loading.method = "spearman"`, `concordance.weights` are `NA` and
  `loadings = raw.loadings`. Positive = upregulated with Y, negative =
  downregulated.

- concordance.weights:

  Matrix (genes/markers x K): per-gene soft concordance weight \\c\_{jk}
  \in \[0,1\]\\ for each component. For gene \\j\\ and component \\k\\,
  \\c\_{jk}\\ is the fraction of the OLS loading numerator
  (\\t\_{k}^\top X\_{c,j}\\) attributable to landmarks where the PLS
  score and the raw density contrast coefficient agree in sign, weighted
  by the posterior probability that each landmark's coefficient is
  genuinely nonzero and concordant: \\s\_{k,i} = \Phi(Y\_{0,i} \cdot
  \operatorname{sign}(t\_{k,c,i}) / \hat\sigma_i)\\, where
  \\\hat\sigma_i = \sqrt{s^{2,\text{post}}\_i} \cdot
  \text{stdev.unscaled}\_{i,c}\\ is the limma eBayes posterior standard
  deviation. \\c\_{jk} = 1\\: loading arises entirely from concordant,
  statistically confident landmarks; \\c\_{jk} = 0\\: loading arises
  entirely from discordant or high-uncertainty landmarks (structural
  counterweight); \\c\_{jk} = 0.5\\: neutral (near-zero loading or
  near-zero/highly-uncertain coefficient). `NA` for
  `.loading.method = "spearman"` (rank-based numerators are not
  additively decomposable). Compare `raw.loadings` to `loadings` to
  assess attenuation of structural balancing effects.

- Y.alignment:

  Numeric vector (K): \|cor(Y, score_k)\| per component

- smoothness:

  Numeric vector (K): Laplacian smoothness per score vector

- Y:

  Numeric vector: the density contrast used, centered to mean zero (Y -
  mean(Y)). Note: this differs from the raw coefficient by a constant
  shift. Use the raw coefficient from the model fit for scientific
  interpretation of absolute density changes.

- Y.mean:

  Numeric: mean of the density contrast used (mean(Y))

- M.local:

  Matrix (optional): the interaction matrix (if .store.M = TRUE)

- params:

  List: parameters used, including: `params$residualize` (logical:
  whether residualization was applied) and `params$nuisance.cols`
  (character: names of nuisance design columns, NULL if not
  residualized)

## Details

plsD is an **interpretive decomposition**, not a formal hypothesis test.
It answers: "What features — through expression patterns, population
identity, or both — explain the density contrast?" Loadings reflect
combined signal, like: markers of depleted populations carry negative
loadings (high expression where density is low), DE genes carry signed
loadings tracking their direction of change, and mixed DA/DE features
appear alongside both.

When `.YX.interaction = TRUE` (default), the interaction term diag(Y)
bakes the density contrast into the data matrix: M.local = P \\ it
captures DA markers, DE genes, and their interplay in a single
decomposition. The interaction amplifies expression signal in landmarks
with strong density contrast, making the method sensitive even to subtle
DE in the absence of DA.

When `.YX.interaction = FALSE`, M.local = P \\ centered expression
without Y-weighting). Y appears only on the response side, so loadings
reflect only features whose graph-smoothed expression independently
covaries with Y. This mode is less comprehensive — it misses DA markers
unless their expression happens to correlate with Y through smoothed
space. In the presence of strongly bimodal density contrasts and/or,
datasets with a small number of extreme landmarks, setting
`.YX.interaction = FALSE` can help distinguish features whose expression
genuinely covaries with Y from those that are strong but purely
geometric counterweights.

The method:

1.  Extracts the density contrast vector Y (coefficients from `get.lm`)

2.  Prepares centered expression matrix Xc (size-factor normalized +
    log2 for RNA; raw for cyto; then centered)

3.  Builds random-walk normalized graph P from SNN (with optional degree
    regularization and lazy walk)

4.  Constructs M.local: if `.YX.interaction = TRUE`, M.local = P \\
    M.local = P \\ Optionally, residualize Xc against nuisance
    covariates before constructing M.local.

5.  Runs NIPALS PLS1: iteratively finds feature weights w maximizing
    cov(M.local w, Y), with deflation of both M.local and Y

6.  Applies sign convention: positive scores = aligned with Y

7.  Computes feature loadings (Pearson r, OLS regression, or Spearman
    rank correlation)

8.  Computes Laplacian smoothness (Sk) for each component

**Diagnostic metrics (Ak, Sk):**

**Ak (Y-alignment):** Measures how strongly component scores covary with
density contrast Y. Ak is high by construction for the leading
components: NIPALS PLS1 maximizes covariance with Y at every deflation
step. Components are ordered by extraction (plsD1 = highest covariance
with Y).

**Sk (graph smoothness):** Derived from the normalized Laplacian applied
to score vectors. High Sk indicates large-scale, graph-smooth structure.

**About Y appearing on both sides (when .YX.interaction = TRUE):** Y
appears in both the data matrix (via diag(Y) in M.local) and the PLS
objective (maximize covariance with Y). This is intentional: it gives
plsD comprehensive sensitivity to features driving density changes. The
design is not circular because the expression matrix Xc sits between:
the product is large only when features exist whose Y-weighted,
graph-smoothed expression genuinely covaries with Y. With permuted Y, Ak
collapses. Setting `.YX.interaction = FALSE` removes Y from the data
matrix entirely, providing a diagnostic comparator.

**Nuisance residualization (.residualize = TRUE):** When the
experimental design includes nuisance covariates (e.g. batch, sex, age),
the density contrast Y from `get.lm` is already adjusted for these.
However, the expression matrix X is not. Setting `.residualize = TRUE`
projects X onto the orthogonal complement of the nuisance design, so
that both sides of the PLS objective are free of nuisance effects. This
improves interpretational symmetry: loadings reflect features whose
expression covaries with the *adjusted* density contrast, with
nuisance-driven expression variation removed. The residualization uses
the Frisch-Waugh-Lovell (FWL) decomposition: X is replaced by (I - H_Z)
X, where H_Z is the hat matrix of the nuisance design expanded to
landmark level via `@config$key`. Sparsity is preserved via implicit
operators. When no nuisance covariates exist (e.g. a simple two-group
design), `.residualize = TRUE` is equivalent to `.residualize = FALSE`.

**Blocking variable handling:** When
[`get.lm()`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md)
was called with a `.block` argument (invoking
[`limma::duplicateCorrelation`](https://rdrr.io/pkg/limma/man/dupcor.html)),
the blocking factor does NOT appear in the design matrix — it enters
only through the GLS covariance structure for the density contrast Y. If
`.residualize = TRUE`, the blocking variable is automatically detected
from `fit$block` and appended to the nuisance design Z as dummy-coded
columns. This ensures that block-level mean expression shifts (e.g.
donor-specific expression baselines) are removed from X, preventing them
from contaminating PLS directions via spurious covariance with the
block-correlated residual structure in Y.

**Structural score constraints and the balancing effect:** The score
vectors computed by plsD are structurally mean-zero across landmarks.
This is a mathematical consequence of column-centering `M.local`: for
any weight vector `w`, the score `t = M.local w` satisfies `sum(t) = 0`.
As a result, a dominant positive pole (landmarks with large positive
scores) must be balanced by compensatory negative scores somewhere on
the manifold.

In datasets where a small number of biologically extreme landmarks
simultaneously have large \|Y_c\| and unusual expression profiles (large
row-wise expression norm), these landmarks can dominate the first PLS
component. The compensatory negative region may then appear coherent and
strong, mimicking a genuine opposing biological program.

**Distinguishing genuine signal from structural counterweight**: the
`plotPlsD` scatter panel colors landmarks by their *raw* (uncentered)
density contrast coefficient. If landmarks with large negative scores
show near-zero or positive raw coefficients, they are geometric
counterweights rather than genuinely depleted populations. Conversely
(depending on contrast direction and component), large-magnitude
positive-score landmarks with near-zero raw Y may also reflect
structural balance.

**Pre-analysis diagnostic**: if unexpected strong balancing is observed,
compute the row-wise interaction norm distribution before running plsD:


    Yc <- coef.mat[, .coef.col]; Yc <- Yc - mean(Yc)
    # Xc: size-factor-normalized, log2-transformed, column-centered expression
    row_norms <- sqrt(rowSums((Yc * as.matrix(Xc))^2))
    quantile(row_norms)

A heavily right-skewed distribution (max \>\> Q75, or top 10 landmarks
holding \>20% of Frobenius norm squared) indicates dominant leverage;
reduce `.lazy.alpha` or winsorize `Yc` before proceeding.

## Note

plsD is designed as an exploratory tool to help interpret density
changes in terms of the features driving them. It does not provide
gene-level p-values or formal multiple testing correction. For rigorous
differential expression testing with p-values following field standards,
use
[`get.pbDE`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md)
(pseudobulk DE via edgeR/limma, both design and marker modes). plsD
complements these methods by providing a multivariate, graph-aware
decomposition that captures joint DA/DE patterns not visible in
gene-by-gene tests.

## See also

[`get.lm`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md)
(required predecessor),
[`plotPlsD`](https://opensource.nibr.com/tinydenseR/reference/plotPlsD.md)
(visualization),
[`plotPlsDHeatmap`](https://opensource.nibr.com/tinydenseR/reference/plotPlsDHeatmap.md)
(heatmap)

## Examples

``` r
if (FALSE) { # \dontrun{
# After fitting linear model
lm.obj <- get.lm(lm.obj, .design = design)

# Run plsD for "Infection" coefficient
lm.obj <- get.plsD(lm.obj, .coef.col = "Infection")

# Access results
lm.obj$plsD$Infection$scores[, "plsD1"]      # PLS scores (Y-aligned)
lm.obj$plsD$Infection$loadings[, "plsD1"]      # concordance-filtered gene loadings
lm.obj$plsD$Infection$raw.loadings[, "plsD1"]  # raw gene loadings before concordance filtering

# Diagnostic table
data.frame(
  component = colnames(lm.obj$plsD$Infection$scores),
  Ak = lm.obj$plsD$Infection$Y.alignment,
  Sk = lm.obj$plsD$Infection$smoothness
)
} # }
```
