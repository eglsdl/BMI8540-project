# BMI8540-project
My project for course BMI8540.

# Objectives and outcomes:
1) **Download data:** Using bash, download pre-mapped and pre-processed CUTAC-seq data (bigWig file format). I will use data from GEO accession GSE224292. **Outcome:** Bash script used to download the data.

2) **Visualize accessible chromatin:** Extract peak information from bigWig file and store it in a BED file. Use deepTools computeMatrix to generate matrix with both peak and signal information. Then, using a custom python script, visualize peak enrichment data in the form of a heatmap. **Outcome:** Python script used to plot the heatmaps.

3) **Store information about peaks with highest signal in a database:** Generate a script which stores information about 5% of peaks with highest signal in the SQL database (Note: the percentage of peaks might change). **Outcome:** SQL script used to store peak information in the database.
