# scoverplots
Reproduce plots for scover manuscript. Output will be generated in `output` dir. For the scover repository with model code, please see [https://github.com/jacobhepkema/scover](https://github.com/jacobhepkema/scover)

## Human kidney

First, make sure to untar the following files in `data/humankidney`:

```
20230123_kidney_all_loo_scores_df.csv.tar.gz
20230125_kidney_repr_aligned_motif_patterns_in_genes.csv.tar.gz
```

```bash
# Work in progress to be updated
Rscript humankidney.R
```

## Tabula Muris

First, make sure to untar the following files in `data/tm`:
```
20221111_tm_all_loo_scores_df.csv.tar.gz
20230126_tm_repr_aligned_motif_patterns_in_genes.csv.tar.gz
```

```bash
# Work in progress to be updated
Rscript tabulamuris.R
```

## Human brain
```bash
# Work in progress, will be added soon
Rscript humanbrain.R
```
