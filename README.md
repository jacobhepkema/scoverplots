# scoverplots
Reproduce plots for scover manuscript. Output will be generated in `output` dir. For the scover repository with model code, please see [https://github.com/jacobhepkema/scover](https://github.com/jacobhepkema/scover)

## Setup

This repository uses Git LFS for file storage. See [this page](https://git-lfs.com/) for more information.

First, make sure Git LFS is installed, and after cloning the repo, enter and run:
```
git lfs pull
```
This will download the large files associated with the repository.

To untar the large files, run the following command:
```
bash untar_files.sh
```

Now, you can the figure creation using the following commands.

## Human kidney

```bash
Rscript humankidney.R
```

## Tabula Muris

```bash
Rscript tabulamuris.R
```

## Human brain

```bash
Rscript humanbrain.R
```
