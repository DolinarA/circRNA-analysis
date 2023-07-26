# Analysis of circRNA expression

This analysis is written for Agilent Human circRNA array V2 (8x15K) microarray. It uses a dataset [GSE101684](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101684) as an example dataset for the analysis. This dataset from Jiangmin Zhao, Le Li, Qian Wang, Hongxiu Han, Qing Zhan, and Ming Xu analyses the circRNA expression profiles in early stage lung adenocarcinoma ([article](https://doi.org/10.1159/000485953)). Some other publicly available datasets for this platform are available [here](https://https.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL21825).

Based on your selected dataset, you need to change two variables, `my_colors` in row 24 and `factors` (from row 42 onward). Raw data should be in folder `RawData` and phenotypic data should be in file `samples_phenoData.txt`.

In the folder `Annotation`, you can find a compiled annotation file. Several sources of annotations were combined to produce a more complete annotation file. Sources can be found in `Annotation files V2.txt`.
