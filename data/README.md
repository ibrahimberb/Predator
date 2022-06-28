<br/>
<h1 align="center">Predator Data</h1>
<br/>

The _data_ folder contains all the datasets to reproduce the results discussed 
in the article.

## Input

### IntAct Mutations Influencing Interactions 
IntAct data contains annotations of experimental evidence where mutations have been 
shown to affect a protein interaction. The data is downloaded from 
[here](https://www.ebi.ac.uk/intact/download/datasets#mutations) (released version 
July, 2020).

The dataset after filtration of human mutations and as a result of the retrieval 
of features from [ELASPIC web server](http://elaspic.kimlab.org) can be found 
in _training_data_M1.txt_.

### Cancer Somatic Mutations (Single Nucleotide Variations)
The cancer mutation datasets (Single Nucleotide Variations) are obtained from TCGA 
project. These datasets are downloaded from [GDC](https://portal.gdc.cancer.gov) 
using R/Bioconductor 
_[TCGAbiolinks](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html)_ 
package on August 11th, 2020.

The mutation datasets include these six cohorts:
1. Breast invasive carcinoma (BRCA)
2. Colon adenocarcinoma (COAD)
3. Esophageal carcinoma (ESCA)
4. Glioblastoma multiforme (GBM)
5. Head-Neck squamous cell carcinoma (HNSC)
6. Ovarian serous cystadenocarcinoma (OV)

### ELASPIC Merged Results

We retrieve the features for mutations present in cancer somatic mutation datasets 
using [ELASPIC web server](http://elaspic.kimlab.org). For each cohort, there are 
two type of mutations: _core_ and _interface_.

### CGC Genes
[The Cancer Gene Census (CGC)](https://pubmed.ncbi.nlm.nih.gov/30371878/) genes 
are reference sets of known cancer genes.

## Output

### Patient Interaction Datasets

Patient interaction datasets are generated Excel files containing information number of disrupted
and nondisrupted interactions across patients for each gene.

### Predictions Datasets
The Predator's predictions for each _(protein, mutation, interactor)_ triplet can be found
in the _predictions_soft_<date>.csv file across each TCGA cohort. The _Prediction_ is the binary classification results, being either 0 (nondisruptive) 
or 1 (disruptive). The _Median_Probability_ is the probability of being class 1 with a
threshold of 0.50. Any triplets above 0.50 are assigned 1 and those with a probability 
lower than or equal to 0.50 are assigned 0.

The _<TCGA\>\_preliminary_data_cgc\_<date\>.xlsx_ files are Excel files generated as a 
result of the analysis. In it, several counts such as the number of entries and 
number patients are listed for each proteins (and genes).

### References

* [Capturing variation impact on molecular interactions in the IMEx Consortium mutations data set](https://www.nature.com/articles/s41467-018-07709-6), _Nature Communications_, 2019

* [ELASPIC web-server: proteome-wide structure-based prediction of mutation effects on protein stability and binding affinity](https://academic.oup.com/bioinformatics/article/32/10/1589/1743335?login=true), _Bioinformatics_, 2016

* [Comprehensive genomic characterization defines human glioblastoma genes and core pathways](https://www.nature.com/articles/nature07385), _Nature_, 2008

* [TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data](https://academic.oup.com/nar/article/44/8/e71/2465925?login=true), Nucleic Acids Research, 2016
