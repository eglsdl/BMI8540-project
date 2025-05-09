CREATE TABLE my_heatmap_records (
    srr_id VARCHAR(11),
    total_reads INT,
    trimmed_reads INT,
    aligned_reads INT,
    alignment_rate VARCHAR(6),
    peak_fn VARCHAR(100),
    plot_fn VARCHAR(100),
    date VARCHAR(10),
    PRIMARY KEY(srr_id)
);
INSERT INTO my_heatmap_records VALUES("SRR23310241", 6676989, 6517213, 5220393, "80.10%", "outputs/SRR23310242.peaks.stringent.bed", "outputs/plots/SRR23310241.png", "2025-05-09");
INSERT INTO my_heatmap_records VALUES("SRR23310242", 9170828, 9156682, 8551013, "93.39%", "outputs/SRR23310242.peaks.stringent.bed", "outputs/plots/SRR23310242.png", "2025-05-09");
