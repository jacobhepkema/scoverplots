#!/bin/bash
tar -xzvf data/humankidney/20230123_kidney_all_loo_scores_df.csv.tar.gz
tar -xzvf data/humankidney/20230125_kidney_repr_aligned_motif_patterns_in_genes.csv.tar.gz
tar -xzvf data/tm/20221111_tm_all_loo_scores_df.csv.tar.gz
tar -xzvf data/tm/20230126_tm_repr_aligned_motif_patterns_in_genes.csv.tar.gz
tar -xzvf data/humanbrain/20230110_trevino_100neighbours_all_loo_scores_df.csv.tar.gz

mv 20230123_kidney_all_loo_scores_df.csv data/humankidney/
mv 20230125_kidney_repr_aligned_motif_patterns_in_genes.csv data/humankidney/
mv 20221111_tm_all_loo_scores_df.csv data/tm/
mv 20230126_tm_repr_aligned_motif_patterns_in_genes.csv data/tm/
mv 20230110_trevino_100neighbours_all_loo_scores_df.csv data/humanbrain/

