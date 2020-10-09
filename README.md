# SiftCell
**Version 1.0.0**
A robust framework to filter cell-free droplets in scRNA-seq data. For details, please refer to the paper:xxx

---

## Getting Started
Please follow the following steps
- **Prerequisites**
	- Latest version of R is required; Is' better to install the latest version of python, but not required.
- **Installation**
  	- Install the released version of SiftCell from github
  	- 

	```ruby
	install.packages("devtools")
 	library(devtools)
 	install_github("jyxi7676/SiftCell")
	library(SiftCell)
	 ```



---

## Usage
An example is given to show how to run SiftCell framework
- **SiftCell-Shuffle**
 	- example files can be found in , please download this data and put the file in the working directory with folder name "DGE"
	- The function *SiftCellShuffle* takes working directory "workingdir" as input parameters and write shuffleDGE in the working directory. And write barcodes for cell-containing droplets in the working directory
	```ruby
	SiftCellShuffle(workingdir)
	```
- **SiftCell-Boost**
	- After running *SifCell-Shuffle*, two folders(DGE and shuffleDGE) will exist in the working directory
	- As for the example we take default parameters. Please refer to the documentation for details about meaning of the paramters.
	```ruby
	SiftCellBoost(workingdir)
	```
- **SiftCell-Bayes**
	- The function *SiftCellBayes* needs input of external cell type information; The format of the cell type file should have two columns, the first column is barcodes, the 2nd column contains the cell type info.Estimated proportion coeffcient of the Dirichlet Multinomial Mixture model is returned. The example cell type file can be found in xxx. 
	```ruby
	celltype = read.csv("/path_to_celltype.csv",row.names=1)
	prop = SiftCellBayes(workingdir,celltype)
	```
---



## Contact
- Jingyue Xi <jyxi@umich.edu>
- Hyun Min Kang <hmnkang@umich.edu>

---

## License & copyright
Jingyue Xi, Hyun Min Kang, University of Michigan
