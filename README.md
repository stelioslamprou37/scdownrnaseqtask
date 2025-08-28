# scDown-based Single-Cell Downstream Analysis Task

## Project Overview
This submission evaluates the **scDown** R package for integrated downstream analysis of single-cell RNA-seq data. It enables rapid and reproducible analyses of cell-type proportions, cell–cell communication, pseudotime trajectories, and RNA velocity using long-form outputs from Seurat or Scanpy.

## Data Source
The scDown pipeline is designed for *post-processed* single-cell transcriptomic data (e.g., Seurat objects or h5ad files), not raw FASTQ reads. The data for this task were obtained following the access instructions described in the scDown publication (Sun et al., 2025).

## Pipeline Steps
scDown automates the following analyses:

1. **Cell-Type Proportion Analysis** (*scProportionTest*)  
2. **Cell–Cell Communication Networks** (*CellChat*)  
3. **Pseudotime and Trajectory Inference** (*Monocle3*)  
4. **RNA Velocity Modeling** (*scVelo via velociraptor and PAGA*)  

## How to Run
1. Place processed scRNA-seq data (Seurat `.rds` or `.h5ad`) into the `/data` directory.
2. Place the scDown pipeline script or config in the `/workflow` directory.
3. Run scDown via R or command line as documented in the GitHub repository.
4. Answer the task questions using outputs stored in `outputs/` (e.g., cell type counts, trajectory plots).

## References
- **Paper**: *scDown: A Pipeline for Single-Cell RNA‑Seq Downstream Analysis* (Sun et al., 2025)  
- **DOI**: 10.3390/ijms26115297  
- **GitHub**: https://github.com/BCH-RC/scDown
