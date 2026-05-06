# Quick Start with tinydenseR

## Introduction

This vignette provides a minimal end-to-end example of a `tinydenseR`
analysis using simulated trajectory data shipped with the package. You
will learn how to:

1.  Run the core TDR pipeline (`RunTDR`)
2.  Fit a linear model across conditions (`get.lm`)
3.  Compute sample-level embeddings (`get.embedding`)
4.  Visualise unsupervised and supervised sample embeddings
5.  Perform plsD to find features driving the density contrast
    (`get.plsD`)

## Setup

``` r

library(tinydenseR)
library(ggplot2)
```

## Load example data

The package includes a simulated trajectory dataset stored as a
`SingleCellExperiment` object.

``` r

data(sim_trajectory_tdr, package = "tinydenseR")

sim_trajectory <- sim_trajectory_tdr$sce
rm(sim_trajectory_tdr)
```

## Run the TDR pipeline

`RunTDR` performs the density estimation and feature construction steps.
Here we specify the sample variable, assay type, and the number of
highly variable genes to use.

``` r

sim_trajectory <-
  tinydenseR::RunTDR(
    x = sim_trajectory,
    .sample.var = "Sample",
    .assay.type = "RNA",
    .nHVG = 500,
    .verbose = FALSE,
    .seed = 123 # for reproducibility
  )
```

## Differential density analysis

### Define the design and contrasts

We set up a design matrix and contrast to compare Condition B vs.
Condition A using `limma`-style syntax.

``` r

.design <-

  model.matrix(~ 0 + Condition,
               data = tinydenseR::GetTDR(sim_trajectory)@metadata)

.contrasts <-
  limma::makeContrasts(
    ConditionB - ConditionA,
    levels = .design
  )
```

### Fit the linear model

``` r

sim_trajectory <-
  tinydenseR::get.lm(
    x = sim_trajectory,
    .design = .design,
    .contrasts = .contrasts,
    .verbose = FALSE
  )
```

### Compute sample embeddings

``` r

sim_trajectory <-
  tinydenseR::get.embedding(
    x = sim_trajectory,
    .contrast.of.interest = "ConditionB - ConditionA",
    .verbose = FALSE
  )
```

## Visualisation

### Unsupervised sample PCA

``` r

tinydenseR::plotSampleEmbedding(
    x = sim_trajectory,
    .embedding = "pca",
    .color.by = "Condition",
    .cat.feature.color = tinydenseR::Color.Palette[1, c(1, 2)],
    .panel.size = 1.5,
    .point.size = 3
  ) +
  ggplot2::labs(title = "PCA") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
```

### Supervised partial-effect PC

``` r

tinydenseR::plotSampleEmbedding(
    x = sim_trajectory,
    .embedding = "pePC",
    .sup.embed.slot = "ConditionB - ConditionA",
    .x.by = "Condition",
    .color.by = "Condition",
    .cat.feature.color = tinydenseR::Color.Palette[1, c(1, 2)],
    .panel.size = 1.5,
    .point.size = 3
  ) +
  ggplot2::labs(title = "partial-effect PC") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
```

## plsD

### Run plsD

``` r

sim_trajectory <-
  tinydenseR::get.plsD(
    x = sim_trajectory,
    .coef.col = "ConditionB - ConditionA",
    .verbose = FALSE
  )
```

### Density PCA coloured by fold-change

``` r

tinydenseR::plotPCA(
    x = sim_trajectory,
    .feature = tinydenseR::GetTDR(sim_trajectory)@results$lm$default$fit$coefficients[, "ConditionB - ConditionA"],
    .panel.size = 1.5,
    .point.size = 1,
    .color.label = "density\nlog2(+0.5)FC",
    .midpoint = 0
  )
```

### plsD dimension plot

``` r

tinydenseR::plotPlsD(
    x = sim_trajectory,
    .coef.col = "ConditionB - ConditionA",
    .plsD.dim = 1,
    .embed = "pca",
    .panel.size = 1.5
  )[[1]] +
  ggplot2::theme(legend.position = "right")
```

### plsD heatmap

``` r

tinydenseR::plotPlsDHeatmap(
  x = sim_trajectory,
  .coef.col = "ConditionB - ConditionA",
  .plsD.dim = 1,
  .order.by = "plsD.dim",
  .panel.height = 3,
  .feature.font.size = 4
)
```

## Session info

``` r

sessionInfo()
```
