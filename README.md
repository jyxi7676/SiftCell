# SiftCell
SiftCell is a robust framework to identify and filter cell-missing droplets from scRNA-seq data. The SiftCell R package includes three tools: *SiftCell-Shuffle* to visualize cell-containing and cell-missing droplets in manifold space with randomization, *SiftCell-Boost* to classify between two types of droplets, and *SiftCell-Mix* to quantify the contribution of ambient RNAs for each droplet.  For details, please refer to [SiftCell paper](http://xxx)

---

## Getting Started
STtools is a R package, and it is recommended to run in R with version >=4.1.1.See Installtion for required software tools to run STTools. 
A [tutorial](https://colab.research.google.com/drive/1ebV42lohhkTWGjkHHjUs1Ma92IXAxx2t) is shown on brain nuclei data.
Explanation of separate main functions are shown as following: 
- **SiftCell-Shuffle**
	- It works with arbitrary digital expression matrix to visualize the distribution of potentially cell-missing barcoded droplets in a manifold space using randomization. 
 	- It takes working directory "workingdir" as input parameters and write shuffleDGE with barcodes.tsv, features.tsv, and matrix.mtx in the working directory. 
	```ruby
	SiftCellShuffle(workingdir)
	```
- **SiftCell-Boost**
	- It leverages results from *SiftCell-Shuffle* and applies a gradient boosting classification algorithm [XGBoost](https://www.kdd.org/kdd2016/papers/files/rfp0697-chenAemb.pdf) by assigning randomized droplets as negative labels (representing ambient RNAs) and droplets confidently predicted to contain cells as positive labels using an overdispersion test.
	- It takes working directory as input, and output a txt file with cell-containing droplets.
	- It allows the user to set flag genes to to avoid including unintended cell types (i.e. platelets), see pacakge documentation for details.
	- It also allows a manually curated version of training labels by leveraging visualizations from SiftCell-Shuffle. See example format of [manual version of training labels](./examples/manual_labels.csv)  manual version of training labels.
	```ruby
	SiftCellBoost(workingdir)
	```
- **SiftCell-Bayes**
	- It is a model based method that estimate the contribution of ambient RNAs per droplet
	- It takes input of external cell type information either from *SiftCell-Boost* (needs to modify func to take SiftCellBoost result) or from external sources. Please refer to the example format of [cell type](./examples/celltype.csv) files.
	- The format of the cell type file should have two columns, the first column is barcodes, the 2nd column contains the cell type info.Estimated proportion coeffcient of the Dirichlet Multinomial Mixture model is returned. 
	- It outputs the estimated proportion coefficient of each cell type (first k columns for k various cell types)and the ambient cell types(last column).
	```ruby
	celltype = read.csv("/path_to_celltype.csv",row.names=1)
	prop = SiftCellBayes(workingdir,celltype)
	```

## Installation
You also need to install the following software tools 
- R version >= 4.1.1
- Python version >= 2.7.0

To install STtools:

```ruby
install.packages("devtools")
library(devtools)
install_github("jyxi7676/SiftCell")
library(SiftCell)
 ```
---

## Contact
- Jingyue Xi <jyxi@umich.edu>
- Hyun Min Kang <hmnkang@umich.edu>

---

## License & copyright
Jingyue Xi, Hyun Min Kang, University of Michigan
