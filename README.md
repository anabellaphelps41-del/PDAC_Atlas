# PDAC Atlas

Rebuild pipeline for an integrated human PDAC single-cell atlas (R / Seurat on
the Pitt CRC SLURM cluster).

The code is in [`code/`](code/) — see [`code/README.md`](code/README.md) for
the steps, how to run them, and the decisions they encode (which samples are
included, how gene symbols are harmonized, doublet detection, storage).

Data and objects are not in this repository. They live on the cluster at
`/ix/tbruno/anp795/PDAC_Atlas/` (the source atlas is on Zenodo record 21795216,
fetched with `code/fetch_atlas.sh`).
