# scripts to generate the data in extdata folder

## Zebrafish annotations: fish.rds

The `fish.rds` is GRanges object with 22802 ranges and gene_name as metadata.
The names of GRanges object is the ensembl gene IDs.
The data were created by `fish_homologs.R` at the date Oct 11, 2024.

## Zebrafish homologs annotations: homologs.rds

The `homologs.rds` is a list of homologous genes of 
'human', 'house mouse', 'Japanese medaka', 'turquoise killifish',
'three-spined stickleback' for zebrafish.
It is a list of GRanges objects with gene_name and homolog_ensembl_gene_ids as
metadata. The names of GRAnges object is the zebrafish ensembl gene IDs.
The corresponding homolog ensembl gene IDs for each species are saved as
homolog_ensembl_gene_ids in metadata.
The data were created by `fish_homologs.R` at the date Oct 11, 2024.

## orthoFinder sample data: 

The files in the subfolder orthofinder is the toy data to demo for 
ortholog group assignments.
It was created by the script `orthofinder.sh`, and
`createOrthogroups_grange.obj.rds.R`.

## speices name mapping data: sysdata.rda

The `sysdata.rda` in R folder saved the species names used for 
`guessSpecies` function. It was created by the script `createNameTable.R` at
the date Jul 11, 2025.
