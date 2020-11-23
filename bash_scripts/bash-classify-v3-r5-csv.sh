#!bin/bash
./classify_editing_homology_v3_csv.x editing_homology_Cirri_b_R ../../Bla_annot-FINAL_v4_names_tr.fa ../../Bla_annot-FINAL_v4_names.gtf -r 5 -p classified_editing_r5 2>> classifyR5_CSV.log
./classify_editing_homology_v3_csv.x editing_homology_Embr8h_b_ ../../Bla_annot-FINAL_v4_names_tr.fa ../../Bla_annot-FINAL_v4_names.gtf -r 5 -p classified_editing_r5 2>> classifyR5_CSV.log
./classify_editing_homology_v3_csv.x editing_homology_Embr10h_b_ ../../Bla_annot-FINAL_v4_names_tr.fa ../../Bla_annot-FINAL_v4_names.gtf -r 5 -p classified_editing_r5 2>> classifyR5_CSV.log
./classify_editing_homology_v3_csv.x editing_homology_Embr15h_b_ ../../Bla_annot-FINAL_v4_names_tr.fa ../../Bla_annot-FINAL_v4_names.gtf -r 5 -p classified_editing_r5 2>> classifyR5_CSV.log
./classify_editing_homology_v3_csv.x editing_homology_Embr36h_b_ ../../Bla_annot-FINAL_v4_names_tr.fa ../../Bla_annot-FINAL_v4_names.gtf -r 5 -p classified_editing_r5 2>> classifyR5_CSV.log
./classify_editing_homology_v3_csv.x editing_homology_Epidermis_b_R ../../Bla_annot-FINAL_v4_names_tr.fa ../../Bla_annot-FINAL_v4_names.gtf -r 5 -p classified_editing_r5 2>> classifyR5_CSV.log
./classify_editing_homology_v3_csv.x editing_homology_FemGonads_b_R ../../Bla_annot-FINAL_v4_names_tr.fa ../../Bla_annot-FINAL_v4_names.gtf -r 5 -p classified_editing_r5 2>> classifyR5_CSV.log
./classify_editing_homology_v3_csv.x editing_homology_Gills_b_R ../../Bla_annot-FINAL_v4_names_tr.fa ../../Bla_annot-FINAL_v4_names.gtf -r 5 -p classified_editing_r5 2>> classifyR5_CSV.log
./classify_editing_homology_v3_csv.x editing_homology_Gut_b_ ../../Bla_annot-FINAL_v4_names_tr.fa ../../Bla_annot-FINAL_v4_names.gtf -r 5 -p classified_editing_r5 2>> classifyR5_CSV.log
./classify_editing_homology_v3_csv.x editing_homology_Hepatic_b_ ../../Bla_annot-FINAL_v4_names_tr.fa ../../Bla_annot-FINAL_v4_names.gtf -r 5 -p classified_editing_r5 2>> classifyR5_CSV.log
./classify_editing_homology_v3_csv.x editing_homology_MaleGonads_b_R ../../Bla_annot-FINAL_v4_names_tr.fa ../../Bla_annot-FINAL_v4_names.gtf -r 5 -p classified_editing_r5 2>> classifyR5_CSV.log
./classify_editing_homology_v3_csv.x editing_homology_Muscle_b_R ../../Bla_annot-FINAL_v4_names_tr.fa ../../Bla_annot-FINAL_v4_names.gtf -r 5 -p classified_editing_r5 2>> classifyR5_CSV.log
./classify_editing_homology_v3_csv.x editing_homology_NeuralTube_b_ ../../Bla_annot-FINAL_v4_names_tr.fa ../../Bla_annot-FINAL_v4_names.gtf -r 5 -p classified_editing_r5 2>> classifyR5_CSV.log
./DNA_reads_editing_filter_csv.x classified_editing_r5/Cirri/editing_homology_Cirri.csv /data/regvolution/mzawi/reads/ERR1419085_variants_filtered_coding > classified_editing_r5/Cirri/editing_homology_Cirri_vs_DNA_reads.csv
./DNA_reads_editing_filter_csv.x classified_editing_r5/Embr8h/editing_homology_Embr8h.csv /data/regvolution/mzawi/reads/ERR1419085_variants_filtered_coding > classified_editing_r5/Embr8h/editing_homology_Embr8h_vs_DNA_reads.csv
./DNA_reads_editing_filter_csv.x classified_editing_r5/Embr10h/editing_homology_Embr10h.csv /data/regvolution/mzawi/reads/ERR1419085_variants_filtered_coding > classified_editing_r5/Embr10h/editing_homology_Embr10h_vs_DNA_reads.csv
./DNA_reads_editing_filter_csv.x classified_editing_r5/Embr15h/editing_homology_Embr15h.csv /data/regvolution/mzawi/reads/ERR1419085_variants_filtered_coding > classified_editing_r5/Embr15h/editing_homology_Embr15h_vs_DNA_reads.csv
./DNA_reads_editing_filter_csv.x classified_editing_r5/Embr36h/editing_homology_Embr36h.csv /data/regvolution/mzawi/reads/ERR1419085_variants_filtered_coding > classified_editing_r5/Embr36h/editing_homology_Embr36h_vs_DNA_reads.csv
./DNA_reads_editing_filter_csv.x classified_editing_r5/Epidermis/editing_homology_Epidermis.csv /data/regvolution/mzawi/reads/ERR1419085_variants_filtered_coding > classified_editing_r5/Epidermis/editing_homology_Epidermis_vs_DNA_reads.csv
./DNA_reads_editing_filter_csv.x classified_editing_r5/FemGonads/editing_homology_FemGonads.csv /data/regvolution/mzawi/reads/ERR1419085_variants_filtered_coding > classified_editing_r5/FemGonads/editing_homology_FemGonads_vs_DNA_reads.csv
./DNA_reads_editing_filter_csv.x classified_editing_r5/Gills/editing_homology_Gills.csv /data/regvolution/mzawi/reads/ERR1419085_variants_filtered_coding > classified_editing_r5/Gills/editing_homology_Gills_vs_DNA_reads.csv
./DNA_reads_editing_filter_csv.x classified_editing_r5/Gut/editing_homology_Gut.csv /data/regvolution/mzawi/reads/ERR1419085_variants_filtered_coding > classified_editing_r5/Gut/editing_homology_Gut_vs_DNA_reads.csv
./DNA_reads_editing_filter_csv.x classified_editing_r5/Hepatic/editing_homology_Hepatic.csv /data/regvolution/mzawi/reads/ERR1419085_variants_filtered_coding > classified_editing_r5/Hepatic/editing_homology_Hepatic_vs_DNA_reads.csv
./DNA_reads_editing_filter_csv.x classified_editing_r5/MaleGonads/editing_homology_MaleGonads.csv /data/regvolution/mzawi/reads/ERR1419085_variants_filtered_coding > classified_editing_r5/MaleGonads/editing_homology_MaleGonads_vs_DNA_reads.csv
./DNA_reads_editing_filter_csv.x classified_editing_r5/Muscle/editing_homology_Muscle.csv /data/regvolution/mzawi/reads/ERR1419085_variants_filtered_coding > classified_editing_r5/Muscle/editing_homology_Muscle_vs_DNA_reads.csv
./DNA_reads_editing_filter_csv.x classified_editing_r5/NeuralTube/editing_homology_NeuralTube.csv /data/regvolution/mzawi/reads/ERR1419085_variants_filtered_coding > classified_editing_r5/NeuralTube/editing_homology_NeuralTube_vs_DNA_reads.csv

