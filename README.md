# Batch effect correction in Baseline experiments from Expression Atlas

The pipeline to correct batch effect between Baseline experiments is described in details below. The functions are written in the file `pipeline.R`. They are now available and documented in the [package `eigenangles`](https://github.com/gheager/eigenangles).

## Installing the dependencies (first use)
To run the following functions, you will need some packages that can be installed using these commands in R :
```r
install.packages(c('magrittr','stringr','purrr','ggplot2','igraph','gtools'))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('sva')
BiocManager::install('RUVSeq')
BiocManager::install('batchelor')
```

## Loading data in R
If the data files are already on your computer, you can use this step. If you want to download data from Expression Atlas, skip this part and go directly to "Downloading data from Expression Atlas".

- Put the data files in `.Rdata` format in a directory containing only them. (The R objects contained in those files must be of class `SimpleList` with a `$rnaseq` slot which is a `SummarizedExperiment` object containing the data.)
- Load all the experiments in a list using the function `load_experiments` : 

```r
experiments <- load_experiments('directory_path')
```

## Downloading data from Expression Atlas
If you want to use data from Expression Atlas that can be downloaded in `.Rdata` format, you can use the function `download_experiments_from_ExpressionAtlas` in this way :

```r
experiments <- download_experiments_from_ExpressionAtlas('E-MTAB-3718','E-MTAB-3725','E-GEOD-44366','E-GEOD-74747')
```

The experiments IDs can be found on the bottom left corner of your screen when you browse Expression Atlas experiments and point an experiment.

![](Expression%20Atlas%20screenshot.png)

This downloads the experiments in a new directory called "experiments" in your working directory and loads all the experiments in R within a list, using `load_experiments` function.

After having loaded the experiments, you get a list of `SummarizedExperiment` object :
```
> experiments
$GEOD44366
class: RangedSummarizedExperiment 
dim: 53465 6 
metadata(4): pipeline filtering mapping quantification
assays(1): counts
rownames(53465): ENSMUSG00000000001 ENSMUSG00000000003 ... ENSMUSG00000115849 ENSMUSG00000115850
rowData names(0):
colnames(6): SRR729296 SRR729297 ... SRR729300 SRR729301
colData names(7): AtlasAssayGroup developmental_stage ... strain technical_replicate_group

$GEOD74747
class: RangedSummarizedExperiment 
dim: 45513 9 
metadata(4): pipeline filtering mapping quantification
assays(1): counts
rownames(45513): ENSMUSG00000000001 ENSMUSG00000000003 ... ENSMUSG00000106670 ENSMUSG00000106671
rowData names(0):
colnames(9): SRR2927735 SRR2927736 ... SRR2927742 SRR2927743
colData names(7): AtlasAssayGroup age ... sex strain

$MTAB3718
class: RangedSummarizedExperiment 
dim: 45513 20 
metadata(4): pipeline filtering mapping quantification
assays(1): counts
rownames(45513): ENSMUSG00000000001 ENSMUSG00000000003 ... ENSMUSG00000106670 ENSMUSG00000106671
rowData names(0):
colnames(20): SRR306757 SRR306758 ... SRR306775 SRR306776
colData names(8): AtlasAssayGroup organism ... biosource_provider technical_replicate_group

$MTAB3725
class: RangedSummarizedExperiment 
dim: 45513 6 
metadata(4): pipeline filtering mapping quantification
assays(1): counts
rownames(45513): ENSMUSG00000000001 ENSMUSG00000000003 ... ENSMUSG00000106670 ENSMUSG00000106671
rowData names(0):
colnames(6): SRR579545 SRR579546 ... SRR579549 SRR579550
colData names(6): AtlasAssayGroup organism ... sex strain
```

## Removing the isolated experiments
To correct batch effect, one needs to take the biological characteristics of the samples into account (organism part in our example). If no sample of an experiment shares biological characteristics with samples from other batches, it is not possible to correct batch effect with these batches since one cannot distinguish the biological difference from the artifact. The function `remove_isolated_experiments` removes the isolated experiments and plots graphs of intersections between the experiments before and after removal.

```r
experiments %<>% remove_isolated_experiments('organism_part')
```

WARNING : this function only removes the isolated experiments. Although it is still possible that two or more unconnected groups of experiments remain, within which the experiments are connected. In this case, batch effect correction is not possible neither and one has to choose a group of experiments manually.

The two following plots are displayed by the function. The first one shows the graph of intersections of all the experiments before the removal of isolated ones. The second shows the same graph after their removal.
![](https://github.com/gheager/Pipeline-for-batch-effect-correction-in-Baseline-experiments/blob/master/graph%20before%20removal.png)
![](https://github.com/gheager/Pipeline-for-batch-effect-correction-in-Baseline-experiments/blob/master/graph%20after%20removal.png)

## Merging experiments in a single dataset
The function `merge_experiments` merges all the experiments in the list in a single `SummarizedExperiment` object and doesn't perform any correction. This function has two additional arguments `log` and `filter` (respectively set to `TRUE` and `FALSE` by default).
- The `log` argument determines whether to perform log transformation on the data (recommended).
- The `filter` argument determines whether to filter genes for which all the samples of a batch have zero-counts. Set it to `TRUE` if you have issues in running ComBat at the next step.

```r
experiments %<>% merge_experiments
#OR
experiments %<>% merge_experiments(log=TRUE,filter=FALSE)
#or any other settings
```

`experiments` is now an only `SummarizedExperiment` object containing the information about batches both in its `@metadata` and `@colData` slots :
```
> experiments
class: SummarizedExperiment 
dim: 45513 35 
metadata(1): batch
assays(1): log_counts
rownames(45513): ENSMUSG00000000001 ENSMUSG00000000003 ... ENSMUSG00000106670 ENSMUSG00000106671
rowData names(0):
colnames(35): SRR2927735 SRR2927736 ... SRR579549 SRR579550
colData names(11): AtlasAssayGroup age ... technical_replicate_group batch
```

## Correcting batch effect
This is the most important step. The function `correct_batch_effect` takes three or four arguments (depending on the choice of the correction algorithm).
- `experiment` which is the merged dataset obtained at previous step
- `model` which is a R formula mentioning the biological factor to take into account during correction (still organism part in our example).
- `method` which can be `'ComBat'`, `'RUV'` or `'MNN'` that use respectively the functions `ComBat` from package `sva`, `RUVs` from `RUVSeq` and `mnnCorrect` from `batchelor`
- `k` which is an argument for `RUV` and `MNN`. In `RUV`, it denotes the "number of factors of unwanted variation to remove" while in `MNN`, it denotes the number of nearest neighbours to consider (kNN) (see documentation of these packages for more information).

```r
experiments %<>% correct_batch_effect(model = ~organism_part, method = 'ComBat')
```