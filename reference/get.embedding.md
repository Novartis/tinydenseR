# Compute Sample Embedding from Partial Fitted Values

Computes embeddings for sample-level visualization. When called without
supervised arguments, computes unsupervised embeddings (PCA and
diffusion map trajectory) on landmark densities stored in
`.tdr.obj@density$log.norm`. When called with supervised arguments
(`.contrast.of.interest` or `.red.model`), additionally computes a
partial-effect PCA (pePC) that isolates variation attributable to a
specific effect. Exact for OLS; if duplicateCorrelation/blocking is
used, the decomposition is approximate.

## Usage

``` r
get.embedding(x, ...)

# S3 method for class 'TDRObj'
get.embedding(
  x,
  .full.model = "default",
  .term.of.interest = NULL,
  .red.model = NULL,
  .contrast.of.interest = NULL,
  .n.eigs = 20,
  .n.pcs = 20,
  .ret.trajectory = FALSE,
  .traj.dist.metric = "cosine",
  .seed = 123,
  .verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A
  [`TDRObj`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md),
  Seurat, SingleCellExperiment, or HDF5AnnData (anndataR) object
  processed through
  [`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md).
  Contains `@density$log.norm` (log2-transformed densities) used for
  unsupervised embeddings. Statistical model fits should be stored in
  `@results$lm` via
  [`get.lm()`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md).

- ...:

  Additional arguments passed to methods.

- .full.model:

  Character string naming the full model in `.tdr.obj$map$lm` (default
  "default"). Required only when computing supervised embeddings (pePC).

- .term.of.interest:

  Character string: the covariate/term being isolated when using the
  nested model method (`.red.model`). Must match a column name in
  `.tdr.obj$metadata`. Ignored when using `.contrast.of.interest`.

- .red.model:

  Optional character string naming the reduced model in
  `.tdr.obj$map$lm` (without the effect of interest). If provided,
  computes \\\Delta\hat{Y} = \hat{Y}\_{full} - \hat{Y}\_{red}\\
  directly. Make sure to construct nested models by "dropping terms" (so
  reduced is a strict subset of full)!

- .contrast.of.interest:

  Optional: Character string naming the contrast to extract. Must match
  a column name in the full model's contrasts. When specified, uses the
  Frisch-Waugh-Lovell (FWL) theorem to compute the true partial fitted
  component.

- .n.eigs:

  Integer: number of eigenvectors for diffusion map (default 20).

- .n.pcs:

  Integer: number of PCs for unsupervised PCA and diffusion map (default
  20).

- .ret.trajectory:

  Logical: whether to compute diffusion map trajectory embedding
  (default FALSE).

- .traj.dist.metric:

  Character: distance metric for diffusion map (default "cosine").

- .seed:

  Integer: random seed for reproducibility of sample-level PCA (via
  [`irlba::prcomp_irlba`](https://rdrr.io/pkg/irlba/man/prcomp_irlba.html))
  and diffusion map computation. Default 123.

- .verbose:

  Logical: print progress messages? Default TRUE.

## Value

Modified `.tdr.obj` with embeddings stored in `.tdr.obj$map$embedding`:

- `pca`:

  Unsupervised PCA on log-transformed landmark densities:

  - `coord`: Sample coordinates (samples x PCs)

  - `rotation`: Loadings (landmarks x PCs)

  - `center`, `scale`, `sdev`: Centering, scaling, and standard
    deviations

- `traj`:

  Diffusion map trajectory embedding:

  - `coord`: Sample coordinates (samples x DCs)

  - `eigenvalues`, `eigenvectors`: Raw diffusion map components

  - `sdev`: sqrt(abs(eigenvalues))

- `pePC`:

  Supervised partial-effect embeddings (when pePC args provided):

  - `<contrast>` or `<term>`: Named list containing:

    - `coord`: Nuisance-residualized projection (samples x pePCs)

    - `x`: PCA scores on partial fitted values

    - `rotation`: Landmark loadings from partial-effect PCA

    - `effect.resid.cor`: Per-PC correlation diagnostic

    - `delta.Yhat`: Partial fitted values matrix

    - `method`: "fwl_contrast" or "nested_models"

## Details

**Unsupervised Embeddings (always computed):**

**Method 1: FWL-based contrast extraction** (when
`.contrast.of.interest` is provided)

For a cell-means model with nuisance covariates, uses the
Frisch-Waugh-Lovell theorem:

1.  Compute contrast regressor: \\x_c = X\_{group} \cdot c\\

2.  Residualize against nuisance: \\x\_{c,\perp} = (I - H_Z) x_c\\

3.  Residualize Y against nuisance: \\Y\_\perp = (I - H_Z) Y\\

4.  Compute partial coefficient: \\\gamma_g = (x\_{c,\perp}'
    Y\_{\perp,g}) / (x\_{c,\perp}' x\_{c,\perp})\\

5.  Form \\\Delta\hat{Y} = \gamma_g \otimes x\_{c,\perp}\\

**Method 2: Nested model comparison** (when `.red.model` is provided)

Direct computation: \\\Delta\hat{Y} = \hat{Y}\_{full} - \hat{Y}\_{red}\\

This method is useful when you want to isolate the effect of one or more
covariates without using contrasts (e.g., comparing
`~ disease + sex + age` vs `~ sex + age` to extract the disease effect).

**Embedding via PCA:**

The embedding is computed in two steps:

1.  Learn PCA basis \\V\\ from \\\Delta\hat{Y}^\top\\ (partial fitted
    values):
    `pca <- prcomp(t(delta.Yhat), center = TRUE, scale. = FALSE, rank. = r)`

2.  Project nuisance-residualized samples onto this basis:
    `coord <- t(Y - Yhat.red - pca$center) %*% pca$rotation`

Here \\E\_{red} = Y - \hat{Y}\_{red}\\ are the nuisance-residualized
features (what remains after removing nuisance effects). The loadings
\\V\\ encode the linear subspace (often rank-1 for a single contrast)
that maximizes variance in the partial-effect matrix
\\\Delta\hat{Y}^\top\\.

Geometrically, this evaluates how the residualized samples align with
the effect direction(s). For rank-1 effects (single contrast or
two-level factor), scores increase when residuals align with the unique
partial-effect axis. For rank \> 1 effects, each pePC axis captures a
variance-maximizing direction within the multi-dimensional effect
subspace (see below).

### Multi-level terms (rank \> 1 effects)

When the dropped term has rank \> 1 (e.g., a factor with more than two
levels), the nested model method (`.red.model`) produces
\\\Delta\hat{Y}\\ with rank \\r = p\_{full} - p\_{red} \> 1\\. PCA then
learns \\r\\ orthogonal axes that maximize variance within this effect
subspace.

**Key property:** these axes are rotations that maximize explained
variance. They do *not* correspond to individual level-vs-baseline
comparisons. For example, with a three-level factor (Baseline, D1, D7),
the nested model pePC produces two axes (pePC1, pePC2) that span the
same subspace as the D1-vs-Baseline and D7-vs-Baseline effects, but
pePC1 might capture "shared activation" while pePC2 captures the
difference between D1 and D7.

**The whole-term embedding is correct and parameterization-invariant:**
\\\Delta\hat{Y}\\ depends only on the column spaces of the full and
reduced design matrices, not on the choice of factor coding (treatment,
sum, Helmert) or reference level. The pePC subspace faithfully
represents the complete effect of the dropped term.

**For per-level interpretability, use explicit contrasts instead.** The
FWL-based method (`.contrast.of.interest`) produces rank-1 embeddings
that isolate specific comparisons (e.g., D1 vs Baseline, D7 vs
Baseline). This requires fitting a cell-means model with
[`limma::makeContrasts`](https://rdrr.io/pkg/limma/man/makeContrasts.html).
See Example 3 below.

**Practical guidance:**

- Two-level factors or single contrasts: both methods are equivalent and
  produce a single, interpretable pePC axis.

- Multi-level factors, whole-term view: use nested models. Useful for
  omnibus visualization ("does this factor matter at all?") and for
  downstream analyses that operate on the full effect subspace.

- Multi-level factors, per-level view: use explicit contrasts. Each call
  to `get.embedding()` with a different `.contrast.of.interest` produces
  a separate rank-1 embedding with a clear biological interpretation.

## See also

[`get.lm`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md)
for model fitting,
[`plotSampleEmbedding`](https://opensource.nibr.com/tinydenseR/reference/plotSampleEmbedding.md)
for visualizing the embedding

## Examples

``` r
if (FALSE) { # \dontrun{
# ===========================================================================
# Example 0: Unsupervised embeddings only (before get.lm)
# ===========================================================================

# After get.map, compute unsupervised embeddings
lm.obj <- get.embedding(.tdr.obj = lm.obj)

# Access unsupervised embeddings
lm.obj$map$embedding$pca$coord   # sample PCA coordinates
lm.obj$map$embedding$traj$coord  # sample diffusion map coordinates

# Plot unsupervised embeddings
plotSampleEmbedding(lm.obj, .embedding = "pca", .color.by = "Condition")
plotSampleEmbedding(lm.obj, .embedding = "traj", .color.by = "Timepoint")

# ===========================================================================
# Example 1: Nested model comparison (no contrasts)
# ===========================================================================

# Fit full model with disease effect (stored as "full")
.design <- model.matrix(object = ~ disease + sex + age,
                        data = lm.obj$metadata)
lm.obj <- get.lm(.tdr.obj = lm.obj, .design = .design, 
                    .model.name = "full")

# Fit reduced model without disease (stored as "reduced")
.red.design <- model.matrix(object = ~ sex + age,
                            data = lm.obj$metadata)
lm.obj <- get.lm(.tdr.obj = lm.obj, .design = .red.design, 
                    .model.name = "reduced")

# Extract disease effect embedding (all stored in lm.obj)
lm.obj <- get.embedding(
    .tdr.obj = lm.obj,
    .full.model = "full",
    .term.of.interest = "disease",
    .red.model = "reduced"
)

# Plot the embedding (colored by disease)
plotSampleEmbedding(lm.obj, .embedding = "pePC", .sup.embed.slot = "disease")

# ===========================================================================
# Example 2: FWL-based contrast extraction (cell-means model)
# ===========================================================================

# Fit cell-means model with contrasts
design <- model.matrix(~ 0 + Group + Batch + Age, data = lm.obj$metadata)
contrasts <- limma::makeContrasts(
    TrtVsCtrl = GroupTrt - GroupCtrl,
    KO_TrtVsCtrl = GroupKO_Trt - GroupKO_Ctrl,
    levels = design
)
lm.obj <- get.lm(lm.obj, .design = design, .contrasts = contrasts)

# Extract each contrast embedding
lm.obj <- get.embedding(lm.obj, .contrast.of.interest = "TrtVsCtrl")
lm.obj <- get.embedding(lm.obj, .contrast.of.interest = "KO_TrtVsCtrl")

# Plot the embedding (specify .color.by for metadata column)
plotSampleEmbedding(lm.obj, .embedding = "pePC",
                    .sup.embed.slot = "TrtVsCtrl", .color.by = "Group")

# Access the embedding components (all in lm.obj$map$embedding)
lm.obj$map$embedding$pca$coord              # sample PCA coordinates
lm.obj$map$embedding$traj$coord             # sample trajectory coordinates
lm.obj$map$embedding$pePC$TrtVsCtrl$coord   # sample pePC coordinates
lm.obj$map$embedding$pePC$TrtVsCtrl$rotation # landmark loadings
lm.obj$map$embedding$pePC$TrtVsCtrl$delta.Yhat # partial fitted values

# Multiple embeddings stored together
names(lm.obj$map$embedding)       # c("pca", "traj", "pePC")
names(lm.obj$map$embedding$pePC)  # c("TrtVsCtrl", "KO_TrtVsCtrl")

# ===========================================================================
# Example 3: Multi-level factor — per-level embeddings via explicit contrasts
# ===========================================================================

# Suppose Timepoint has three levels: Baseline, D1, D7.
#
# Option A (nested models): whole-term embedding with rank 2.
# pePC1 and pePC2 are variance-maximizing rotations within the
# Timepoint effect subspace — they do NOT correspond to individual
# level-vs-baseline comparisons.

full.design  <- model.matrix(~ Timepoint + Batch, data = lm.obj$metadata)
red.design   <- model.matrix(~ Batch, data = lm.obj$metadata)
lm.obj <- get.lm(lm.obj, .design = full.design,  .model.name = "tp_full")
lm.obj <- get.lm(lm.obj, .design = red.design,   .model.name = "tp_red")

lm.obj <- get.embedding(lm.obj,
    .full.model = "tp_full", .red.model = "tp_red",
    .term.of.interest = "Timepoint")
# -> produces pePC1, pePC2: omnibus Timepoint effect

# Option B (explicit contrasts): one rank-1 embedding per comparison.
# Each pePC axis directly captures a specific level-vs-baseline effect.

cm.design <- model.matrix(~ 0 + Timepoint + Batch, data = lm.obj$metadata)
tp.contrasts <- limma::makeContrasts(
    D1vsBaseline = TimepointD1  - TimepointBaseline,
    D7vsBaseline = TimepointD7  - TimepointBaseline,
    levels = cm.design
)
lm.obj <- get.lm(lm.obj, .design = cm.design,
                  .contrasts = tp.contrasts,
                  .model.name = "tp_contrasts")

lm.obj <- get.embedding(lm.obj,
    .full.model = "tp_contrasts",
    .contrast.of.interest = "D1vsBaseline")
lm.obj <- get.embedding(lm.obj,
    .full.model = "tp_contrasts",
    .contrast.of.interest = "D7vsBaseline")

# Each embedding is rank 1 with a clear biological meaning:
plotSampleEmbedding(lm.obj, .embedding = "pePC",
    .sup.embed.slot = "D1vsBaseline", .color.by = "Timepoint")
plotSampleEmbedding(lm.obj, .embedding = "pePC",
    .sup.embed.slot = "D7vsBaseline", .color.by = "Timepoint")
} # }
```
