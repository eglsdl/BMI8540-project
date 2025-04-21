# Open Chromatin Peak Enrichment Analysis
My project for course BMI8540

# Background
Chromatin is formed by tightly wrapping the DNA around nucleosomes which causes these parts of the sequence to be transcriptionally inactive. To obtain information about open chromatin, ATAC-seq method is used - DNA is fragmented using Tn5 transposase which cuts it and inserts sequencing adapters into open chromatin regions resulting in small DNA fragments. CUTAC-seq is a similar method that fragments open chromatin surrounding a specific protein target, e.g. RNA Polymerase II, a transcription factor of interest, etc. Using open chromatin peak enrichment analysis we can analyse transcriptional events and identify transcriptionally active regions.

# Objectives and outcomes
1) **Download and prepare data:** Using bash, download ATAC-seq and/or CUTAC-seq data and pre-process it for further analysis (quality control, mapping). **Outcome:** A bash script that is used to download and preprocess the data.

2) **Visualize accessible chromatin:** Peaks of interest will be provided in BED format. For each pair of sequencing data file (as processed in previous objective) and peaks of interest file (BED), deepTools computeMatrix is used to generate matrix with both peak and signal information. Then, using a custom python script, visualize peak enrichment data in the form of a heatmap. **Outcome:** Python script used to plot the heatmaps.

3) **Store information about each heatmap in a database:** Generate a script which stores information about each heatmap in the SQL database. **Outcome:** SQL script used to store information about generated heatmaps in the database.

# Run information
The pipeline is set up to run on HCC Swan server where the used software packages and Slurm Workload Manager is available.

# References
Brahma, S., Henikoff, S. The BAF chromatin remodeler synergizes with RNA polymerase II and transcription factors to evict nucleosomes. Nat Genet 56, 100–111 (2024). https://doi.org/10.1038/s41588-023-01603-8

Grandi, F.C., Modi, H., Kampman, L. et al. Chromatin accessibility profiling by ATAC-seq. Nat Protoc 17, 1518–1552 (2022). https://doi.org/10.1038/s41596-022-00692-9

Meers, M.P., Tenenbaum, D. & Henikoff, S. Peak calling by Sparse Enrichment Analysis for CUT&RUN chromatin profiling. Epigenetics & Chromatin 12, 42 (2019). https://doi.org/10.1186/s13072-019-0287-4
