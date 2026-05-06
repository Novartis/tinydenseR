# Package index

## Core Pipeline

Entry points for running tinydenseR and managing TDRObj objects

- [`tinydenseR`](https://opensource.nibr.com/tinydenseR/reference/tinydenseR-package.md)
  [`tinydenseR-package`](https://opensource.nibr.com/tinydenseR/reference/tinydenseR-package.md)
  : tinydenseR: Linking Cell-To-Cell Variation to Sample-to-Sample
  Variation
- [`RunTDR()`](https://opensource.nibr.com/tinydenseR/reference/RunTDR.md)
  : Run the full tinydenseR pipeline
- [`GetTDR()`](https://opensource.nibr.com/tinydenseR/reference/GetTDR.md)
  : Extract a TDRObj from a container object
- [`SetTDR()`](https://opensource.nibr.com/tinydenseR/reference/SetTDR.md)
  : Store a TDRObj inside a container object
- [`` `$`( ``*`<TDRObj>`*`)`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  [`` `$<-`( ``*`<TDRObj>`*`)`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  [`names(`*`<TDRObj>`*`)`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  [`show(`*`<TDRObj>`*`)`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  : TDRObj: S4 class for tinydenseR analysis objects
- [`TDRObj()`](https://opensource.nibr.com/tinydenseR/reference/TDRObj.md)
  : Construct a TDRObj
- [`setup.tdr.obj()`](https://opensource.nibr.com/tinydenseR/reference/setup.tdr.obj.md)
  [`setup.lm.obj()`](https://opensource.nibr.com/tinydenseR/reference/setup.tdr.obj.md)
  : Initialize tinydenseR object for landmark-based analysis
- [`is.TDRObj()`](https://opensource.nibr.com/tinydenseR/reference/is.TDRObj.md)
  : Check if an object is a TDRObj

## Differential Density Analysis

Linear modeling of landmark-level density across conditions

- [`get.lm()`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md)
  : Differential Density Testing
- [`get.dea()`](https://opensource.nibr.com/tinydenseR/reference/get.dea.md)
  **\[deprecated\]** : Pseudobulk Differential Expression Analysis
  (Deprecated)
- [`get.density()`](https://opensource.nibr.com/tinydenseR/reference/get.density.md)
  : Access fuzzy density matrices from a TDRObj

## Pseudobulk Differential Expression

Manifold-informed pseudobulk DE in design and marker modes

- [`get.pbDE()`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md)
  : Pseudobulk Differential Expression Analysis
- [`get.markerDE()`](https://opensource.nibr.com/tinydenseR/reference/get.markerDE.md)
  **\[deprecated\]** : Marker Gene/Protein Identification (Deprecated)
- [`get.marker()`](https://opensource.nibr.com/tinydenseR/reference/get.marker.md)
  : Deprecated: Use get.pbDE(.mode = "marker") instead

## plsD Decomposition

Graph-diffused partial least squares decomposition

- [`get.plsD()`](https://opensource.nibr.com/tinydenseR/reference/get.plsD.md)
  : Graph-Diffused, Density Contrast-Aligned PLS Decomposition (plsD)

## Sample Embedding

PCA, diffusion-map, and partial-effect PC embeddings

- [`get.embedding()`](https://opensource.nibr.com/tinydenseR/reference/get.embedding.md)
  : Compute Sample Embedding from Partial Fitted Values

## Cell Typing and Clustering

Cell type annotation, clustering, and resolution management

- [`celltyping()`](https://opensource.nibr.com/tinydenseR/reference/celltyping.md)
  : Manually assign cell type labels to clusters
- [`set_active_celltyping()`](https://opensource.nibr.com/tinydenseR/reference/set_active_celltyping.md)
  : Set active celltyping solution
- [`list_celltyping_solutions()`](https://opensource.nibr.com/tinydenseR/reference/list_celltyping_solutions.md)
  : List stored celltyping solutions
- [`import_cell_annotations()`](https://opensource.nibr.com/tinydenseR/reference/import_cell_annotations.md)
  : Import multiple cell-level annotation columns as landmark-level
  celltyping solutions
- [`recluster()`](https://opensource.nibr.com/tinydenseR/reference/recluster.md)
  : Recluster landmarks
- [`set_active_clustering()`](https://opensource.nibr.com/tinydenseR/reference/set_active_clustering.md)
  : Set active clustering solution
- [`leiden.cluster()`](https://opensource.nibr.com/tinydenseR/reference/leiden.cluster.md)
  : Leiden clustering with straggler absorption
- [`lm.cluster()`](https://opensource.nibr.com/tinydenseR/reference/lm.cluster.md)
  : Leiden clustering of landmarks

## Accessors and Utilities

Extract data from containers and helper functions

- [`get.cells()`](https://opensource.nibr.com/tinydenseR/reference/get.cells.md)
  : Create .cells Object with Automatic Format Detection
- [`get.cells.SCE()`](https://opensource.nibr.com/tinydenseR/reference/get.cells.SCE.md)
  : Create .cells Object from SingleCellExperiment Object
- [`get.cells.Seurat()`](https://opensource.nibr.com/tinydenseR/reference/get.cells.Seurat.md)
  : Create .cells Object from Seurat v4 Object
- [`get.cells.Seurat5()`](https://opensource.nibr.com/tinydenseR/reference/get.cells.Seurat5.md)
  : Create .cells Object from Seurat v5 Object
- [`get.cells.list.mat()`](https://opensource.nibr.com/tinydenseR/reference/get.cells.list.mat.md)
  : Create .cells Object from Count Matrices
- [`get.cellmap()`](https://opensource.nibr.com/tinydenseR/reference/get.cellmap.md)
  : Access per-cell map data for a single sample
- [`get.features()`](https://opensource.nibr.com/tinydenseR/reference/get.features.md)
  : Graph-based feature discovery for landmarks
- [`get.graph()`](https://opensource.nibr.com/tinydenseR/reference/get.graph.md)
  : Graph embedding of landmarks
- [`get.landmarks()`](https://opensource.nibr.com/tinydenseR/reference/get.landmarks.md)
  : Identify landmarks via leverage score sampling
- [`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
  : Mapping cells to landmarks
- [`get.meta()`](https://opensource.nibr.com/tinydenseR/reference/get.meta.md)
  : Extract Sample-Level Metadata with Automatic Format Detection
- [`get.meta.HDF5AnnData()`](https://opensource.nibr.com/tinydenseR/reference/get.meta.HDF5AnnData.md)
  : Extract Sample-Level Metadata from HDF5AnnData Object
- [`get.meta.SCE()`](https://opensource.nibr.com/tinydenseR/reference/get.meta.SCE.md)
  : Extract Sample-Level Metadata from SingleCellExperiment Object
- [`get.meta.Seurat()`](https://opensource.nibr.com/tinydenseR/reference/get.meta.Seurat.md)
  : Extract Sample-Level Metadata from Seurat v4 Object
- [`get.meta.Seurat5()`](https://opensource.nibr.com/tinydenseR/reference/get.meta.Seurat5.md)
  : Extract Sample-Level Metadata from Seurat v5 Object
- [`get.adj.matrix()`](https://opensource.nibr.com/tinydenseR/reference/get.adj.matrix.md)
  : Sparse matrix representation of nearest neighbors
- [`get.subset()`](https://opensource.nibr.com/tinydenseR/reference/get.subset.md)
  : Create a hierarchical subset TDRObj
- [`goi.summary()`](https://opensource.nibr.com/tinydenseR/reference/goi.summary.md)
  : Summarize gene expression patterns for genes of interest
- [`is.hpc()`](https://opensource.nibr.com/tinydenseR/reference/is.hpc.md)
  : Detect if running on High-Performance Computing cluster (internal)
- [`fast.jaccard.r()`](https://opensource.nibr.com/tinydenseR/reference/fast.jaccard.r.md)
  : Shared nearest neighbors via fast Jaccard index calculation
- [`elbow.sec.deriv()`](https://opensource.nibr.com/tinydenseR/reference/elbow.sec.deriv.md)
  : Find elbow point using second derivative method (internal)

## Visualization

Plotting functions for landmarks, samples, and results

- [`plotUMAP()`](https://opensource.nibr.com/tinydenseR/reference/plotUMAP.md)
  : Plot UMAP Plot UMAP
- [`plotPCA()`](https://opensource.nibr.com/tinydenseR/reference/plotPCA.md)
  : Plot PCA
- [`plotDensity()`](https://opensource.nibr.com/tinydenseR/reference/plotDensity.md)
  : Plot Density
- [`plotDEA()`](https://opensource.nibr.com/tinydenseR/reference/plotDEA.md)
  **\[deprecated\]** : Plot Differential Expression Analysis Results
  (Deprecated)
- [`plotBeeswarm()`](https://opensource.nibr.com/tinydenseR/reference/plotBeeswarm.md)
  : Bee Swarm Plot of Density Estimate Change
- [`plotHeatmap()`](https://opensource.nibr.com/tinydenseR/reference/plotHeatmap.md)
  : Plot Mean Expression Heatmap
- [`plotMarkerDE()`](https://opensource.nibr.com/tinydenseR/reference/plotMarkerDE.md)
  : Plot Marker DE Results
- [`plotPbDE()`](https://opensource.nibr.com/tinydenseR/reference/plotPbDE.md)
  : Plot Pseudobulk Differential Expression Results
- [`plotPlsD()`](https://opensource.nibr.com/tinydenseR/reference/plotPlsD.md)
  : Plot plsD Scores (Diagnostic and Component Views)
- [`plotPlsDHeatmap()`](https://opensource.nibr.com/tinydenseR/reference/plotPlsDHeatmap.md)
  : Plot plsD Expression Heatmap
- [`plotSampleEmbedding()`](https://opensource.nibr.com/tinydenseR/reference/plotSampleEmbedding.md)
  : Plot Sample Embedding from get.embedding
- [`plotSamplePCA()`](https://opensource.nibr.com/tinydenseR/reference/plotSamplePCA.md)
  **\[deprecated\]** : Plot Sample PCA (Deprecated)
- [`plotTradPerc()`](https://opensource.nibr.com/tinydenseR/reference/plotTradPerc.md)
  : Plot Traditional Percentages
- [`plotTradStats()`](https://opensource.nibr.com/tinydenseR/reference/plotTradStats.md)
  : Plot Traditional Statistics
- [`plot2Markers()`](https://opensource.nibr.com/tinydenseR/reference/plot2Markers.md)
  : Bidimensional Hexbin Plot for Marker Expression
- [`scatterPlot()`](https://opensource.nibr.com/tinydenseR/reference/scatterPlot.md)
  : Scatter Plot with Feature Coloring

## Data and Simulation

Bundled datasets and simulation functions

- [`sim_trajectory_tdr`](https://opensource.nibr.com/tinydenseR/reference/sim_trajectory_tdr.md)
  : Simulated scRNA-seq trajectory with condition-dependent differential
  abundance
- [`simulate_DA_data()`](https://opensource.nibr.com/tinydenseR/reference/simulate_DA_data.md)
  : Simulate differential abundance (DA) flow cytometry data
- [`simulate_DE_data()`](https://opensource.nibr.com/tinydenseR/reference/simulate_DE_data.md)
  : Simulate differential expression (DE) flow cytometry data
- [`fetch_trajectory_data()`](https://opensource.nibr.com/tinydenseR/reference/fetch_trajectory_data.md)
  : Fetch trajectory simulation dataset from miloR
- [`Color.Palette`](https://opensource.nibr.com/tinydenseR/reference/Color.Palette.md)
  : Color Palette

## Conversion

Convert TDRObj to other Bioconductor containers

- [`as.SummarizedExperiment()`](https://opensource.nibr.com/tinydenseR/reference/as.SummarizedExperiment.md)
  : Convert an object to SummarizedExperiment
- [`as.SummarizedExperiment(`*`<TDRObj>`*`)`](https://opensource.nibr.com/tinydenseR/reference/as.SummarizedExperiment.TDRObj.md)
  : Convert a TDRObj to SummarizedExperiment

## Cache Management

On-disk cache for large per-sample mapping data

- [`tdr_cache_info()`](https://opensource.nibr.com/tinydenseR/reference/tdr_cache_info.md)
  : Print a human-readable summary of the on-disk cache state
- [`tdr_cache_validate()`](https://opensource.nibr.com/tinydenseR/reference/tdr_cache_validate.md)
  : Validate that all on-disk cache files are intact
- [`tdr_cache_cleanup()`](https://opensource.nibr.com/tinydenseR/reference/tdr_cache_cleanup.md)
  : Remove all cached files for a tinydenseR object

## Internal

Internal helper functions (not part of the public API)

- [`Color.Palette`](https://opensource.nibr.com/tinydenseR/reference/Color.Palette.md)
  : Color Palette
- [`GetTDR()`](https://opensource.nibr.com/tinydenseR/reference/GetTDR.md)
  : Extract a TDRObj from a container object
- [`RunTDR()`](https://opensource.nibr.com/tinydenseR/reference/RunTDR.md)
  : Run the full tinydenseR pipeline
- [`SetTDR()`](https://opensource.nibr.com/tinydenseR/reference/SetTDR.md)
  : Store a TDRObj inside a container object
- [`` `$`( ``*`<TDRObj>`*`)`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  [`` `$<-`( ``*`<TDRObj>`*`)`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  [`names(`*`<TDRObj>`*`)`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  [`show(`*`<TDRObj>`*`)`](https://opensource.nibr.com/tinydenseR/reference/TDRObj-class.md)
  : TDRObj: S4 class for tinydenseR analysis objects
- [`TDRObj()`](https://opensource.nibr.com/tinydenseR/reference/TDRObj.md)
  : Construct a TDRObj
- [`as.SummarizedExperiment()`](https://opensource.nibr.com/tinydenseR/reference/as.SummarizedExperiment.md)
  : Convert an object to SummarizedExperiment
- [`as.SummarizedExperiment(`*`<TDRObj>`*`)`](https://opensource.nibr.com/tinydenseR/reference/as.SummarizedExperiment.TDRObj.md)
  : Convert a TDRObj to SummarizedExperiment
- [`celltyping()`](https://opensource.nibr.com/tinydenseR/reference/celltyping.md)
  : Manually assign cell type labels to clusters
- [`fetch_trajectory_data()`](https://opensource.nibr.com/tinydenseR/reference/fetch_trajectory_data.md)
  : Fetch trajectory simulation dataset from miloR
- [`get.cellmap()`](https://opensource.nibr.com/tinydenseR/reference/get.cellmap.md)
  : Access per-cell map data for a single sample
- [`get.cells()`](https://opensource.nibr.com/tinydenseR/reference/get.cells.md)
  : Create .cells Object with Automatic Format Detection
- [`get.cells.SCE()`](https://opensource.nibr.com/tinydenseR/reference/get.cells.SCE.md)
  : Create .cells Object from SingleCellExperiment Object
- [`get.cells.Seurat()`](https://opensource.nibr.com/tinydenseR/reference/get.cells.Seurat.md)
  : Create .cells Object from Seurat v4 Object
- [`get.cells.Seurat5()`](https://opensource.nibr.com/tinydenseR/reference/get.cells.Seurat5.md)
  : Create .cells Object from Seurat v5 Object
- [`get.cells.list.mat()`](https://opensource.nibr.com/tinydenseR/reference/get.cells.list.mat.md)
  : Create .cells Object from Count Matrices
- [`get.density()`](https://opensource.nibr.com/tinydenseR/reference/get.density.md)
  : Access fuzzy density matrices from a TDRObj
- [`get.embedding()`](https://opensource.nibr.com/tinydenseR/reference/get.embedding.md)
  : Compute Sample Embedding from Partial Fitted Values
- [`get.features()`](https://opensource.nibr.com/tinydenseR/reference/get.features.md)
  : Graph-based feature discovery for landmarks
- [`get.graph()`](https://opensource.nibr.com/tinydenseR/reference/get.graph.md)
  : Graph embedding of landmarks
- [`get.landmarks()`](https://opensource.nibr.com/tinydenseR/reference/get.landmarks.md)
  : Identify landmarks via leverage score sampling
- [`get.lm()`](https://opensource.nibr.com/tinydenseR/reference/get.lm.md)
  : Differential Density Testing
- [`get.map()`](https://opensource.nibr.com/tinydenseR/reference/get.map.md)
  : Mapping cells to landmarks
- [`get.marker()`](https://opensource.nibr.com/tinydenseR/reference/get.marker.md)
  : Deprecated: Use get.pbDE(.mode = "marker") instead
- [`get.meta.HDF5AnnData()`](https://opensource.nibr.com/tinydenseR/reference/get.meta.HDF5AnnData.md)
  : Extract Sample-Level Metadata from HDF5AnnData Object
- [`get.meta()`](https://opensource.nibr.com/tinydenseR/reference/get.meta.md)
  : Extract Sample-Level Metadata with Automatic Format Detection
- [`get.meta.SCE()`](https://opensource.nibr.com/tinydenseR/reference/get.meta.SCE.md)
  : Extract Sample-Level Metadata from SingleCellExperiment Object
- [`get.meta.Seurat()`](https://opensource.nibr.com/tinydenseR/reference/get.meta.Seurat.md)
  : Extract Sample-Level Metadata from Seurat v4 Object
- [`get.meta.Seurat5()`](https://opensource.nibr.com/tinydenseR/reference/get.meta.Seurat5.md)
  : Extract Sample-Level Metadata from Seurat v5 Object
- [`get.pbDE()`](https://opensource.nibr.com/tinydenseR/reference/get.pbDE.md)
  : Pseudobulk Differential Expression Analysis
- [`get.plsD()`](https://opensource.nibr.com/tinydenseR/reference/get.plsD.md)
  : Graph-Diffused, Density Contrast-Aligned PLS Decomposition (plsD)
- [`get.subset()`](https://opensource.nibr.com/tinydenseR/reference/get.subset.md)
  : Create a hierarchical subset TDRObj
- [`goi.summary()`](https://opensource.nibr.com/tinydenseR/reference/goi.summary.md)
  : Summarize gene expression patterns for genes of interest
- [`import_cell_annotations()`](https://opensource.nibr.com/tinydenseR/reference/import_cell_annotations.md)
  : Import multiple cell-level annotation columns as landmark-level
  celltyping solutions
- [`is.TDRObj()`](https://opensource.nibr.com/tinydenseR/reference/is.TDRObj.md)
  : Check if an object is a TDRObj
- [`leiden.cluster()`](https://opensource.nibr.com/tinydenseR/reference/leiden.cluster.md)
  : Leiden clustering with straggler absorption
- [`list_celltyping_solutions()`](https://opensource.nibr.com/tinydenseR/reference/list_celltyping_solutions.md)
  : List stored celltyping solutions
- [`lm.cluster()`](https://opensource.nibr.com/tinydenseR/reference/lm.cluster.md)
  : Leiden clustering of landmarks
- [`plot2Markers()`](https://opensource.nibr.com/tinydenseR/reference/plot2Markers.md)
  : Bidimensional Hexbin Plot for Marker Expression
- [`plotBeeswarm()`](https://opensource.nibr.com/tinydenseR/reference/plotBeeswarm.md)
  : Bee Swarm Plot of Density Estimate Change
- [`plotDensity()`](https://opensource.nibr.com/tinydenseR/reference/plotDensity.md)
  : Plot Density
- [`plotHeatmap()`](https://opensource.nibr.com/tinydenseR/reference/plotHeatmap.md)
  : Plot Mean Expression Heatmap
- [`plotMarkerDE()`](https://opensource.nibr.com/tinydenseR/reference/plotMarkerDE.md)
  : Plot Marker DE Results
- [`plotPCA()`](https://opensource.nibr.com/tinydenseR/reference/plotPCA.md)
  : Plot PCA
- [`plotPbDE()`](https://opensource.nibr.com/tinydenseR/reference/plotPbDE.md)
  : Plot Pseudobulk Differential Expression Results
- [`plotPlsD()`](https://opensource.nibr.com/tinydenseR/reference/plotPlsD.md)
  : Plot plsD Scores (Diagnostic and Component Views)
- [`plotPlsDHeatmap()`](https://opensource.nibr.com/tinydenseR/reference/plotPlsDHeatmap.md)
  : Plot plsD Expression Heatmap
- [`plotSampleEmbedding()`](https://opensource.nibr.com/tinydenseR/reference/plotSampleEmbedding.md)
  : Plot Sample Embedding from get.embedding
- [`plotSamplePCA()`](https://opensource.nibr.com/tinydenseR/reference/plotSamplePCA.md)
  **\[deprecated\]** : Plot Sample PCA (Deprecated)
- [`plotTradPerc()`](https://opensource.nibr.com/tinydenseR/reference/plotTradPerc.md)
  : Plot Traditional Percentages
- [`plotTradStats()`](https://opensource.nibr.com/tinydenseR/reference/plotTradStats.md)
  : Plot Traditional Statistics
- [`plotUMAP()`](https://opensource.nibr.com/tinydenseR/reference/plotUMAP.md)
  : Plot UMAP Plot UMAP
- [`recluster()`](https://opensource.nibr.com/tinydenseR/reference/recluster.md)
  : Recluster landmarks
- [`scatterPlot()`](https://opensource.nibr.com/tinydenseR/reference/scatterPlot.md)
  : Scatter Plot with Feature Coloring
- [`set_active_celltyping()`](https://opensource.nibr.com/tinydenseR/reference/set_active_celltyping.md)
  : Set active celltyping solution
- [`set_active_clustering()`](https://opensource.nibr.com/tinydenseR/reference/set_active_clustering.md)
  : Set active clustering solution
- [`setup.tdr.obj()`](https://opensource.nibr.com/tinydenseR/reference/setup.tdr.obj.md)
  [`setup.lm.obj()`](https://opensource.nibr.com/tinydenseR/reference/setup.tdr.obj.md)
  : Initialize tinydenseR object for landmark-based analysis
- [`sim_trajectory_tdr`](https://opensource.nibr.com/tinydenseR/reference/sim_trajectory_tdr.md)
  : Simulated scRNA-seq trajectory with condition-dependent differential
  abundance
- [`simulate_DA_data()`](https://opensource.nibr.com/tinydenseR/reference/simulate_DA_data.md)
  : Simulate differential abundance (DA) flow cytometry data
- [`simulate_DE_data()`](https://opensource.nibr.com/tinydenseR/reference/simulate_DE_data.md)
  : Simulate differential expression (DE) flow cytometry data
- [`tdr_cache_cleanup()`](https://opensource.nibr.com/tinydenseR/reference/tdr_cache_cleanup.md)
  : Remove all cached files for a tinydenseR object
- [`tdr_cache_info()`](https://opensource.nibr.com/tinydenseR/reference/tdr_cache_info.md)
  : Print a human-readable summary of the on-disk cache state
- [`tdr_cache_validate()`](https://opensource.nibr.com/tinydenseR/reference/tdr_cache_validate.md)
  : Validate that all on-disk cache files are intact
