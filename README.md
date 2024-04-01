## PHISGAE

Code and Datasets for "Host specificity-guided graph autoencoder for predicting phage-host interaction"

### Developers

Zhen Xiao (xiaozh@mails.ccnu.edu.cn) and Xingpeng Jiang (xpjiang@mail.ccnu.edu.cn) from School of Computer Science, Central China Normal University.

### Datasets

- data/PHI_all.csv is the dataset with known phage-host pairs between bacteriophages and hosts from different taxa.
- data/phage_ncbi.txt is used to download the gene sequence and protein sequence information of all bacteriophages (fasta file and gb file).
- data/DNA_f is a set of files with features encoded by DNA sequences corresponding to all phages.
- data/Protein_f is a set of files with features encoded by protein sequences corresponding to all phages.
- data/dnapro_f is a set of files with features encoded by DNA sequences and protein sequences corresponding to all phages.
- data/phage_host.csv is an interaction matrix between bacteriophages and hosts(species).


### Environment Requirement
The code has been tested running under Python 3.7 and R 4.2.0. The required packages are as follows:

* numpy == 1.21.5
* pandas == 1.3.5
* tensorflow == 1.15.0
* keras == 2.3.1
* torch == 1.4.0+cpu
* [ExPASy](https://web.expasy.org/translate/)

### Usage
#### 1 Prepare the dataset

```
git clone https://github.com//PHISGAE
cd PHISGAE/code
python dna_f.py   
python pro_f.py   ####compute features derived from DNA and protein sequences
```

Users can use their own data to train prediction models. 

1. Download all phage DNA and protein sequences from [NCBI batchentrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez?), named **phage.fasta** and **phage.gb** into the *data* folder.
2. Run `python dna_f.py` and `python pro_f.py` to compute the features derived from DNA and protein sequences. **Note** If your protein sequences in one **.gb** file, you need to run `python protein_processing.py` to split them into different **.txt** files. For phages missing protein sequence information, you can use the [ExPASy](https://web.expasy.org/translate/) tool to translate their gene sequence into protein sequence.
3. Use `python sim.py` and module *## Other similarity calculation method ##* in `data_processing.R` to calculate the different similarity matrix between phages with `data/dna_pro_f/phage_dna_pro.csv`.
4. Use module *KNN graph* in `data_processing.R` to obtain refined phage-phage connection network.
5. Use `python hostI.py` to generate a represention of host-host connection network.


### 2 Predicting phage-host interactions

Users need to transform the known phage-host pairs into different interaction matrices by running Run `python RAWA.py`, then  using `python main.py` to predict phage-host interactions.

**Example**

You can use the files in *example* folder as the input and run `python main.py` to predict phage-host interactions.


### Contact

Please feel free to contact us if you need any help.
