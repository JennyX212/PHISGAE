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
- data/phage_host.csv is a dataset representing the interaction matrix between bacteriophages and hosts.
- data/phage_sim.csv is a dataset representing the refined phage-phage connection network.
- data/host_I.csv is is a dataset representing the host-host connection network.


### Environment Requirement
The code has been tested running under Python 3.7 and R 4.2.0. The required packages are as follows:

* numpy == 1.21.5
* pandas == 1.3.5
* tensorflow == 1.15.0 
* keras == 2.3.1
* torch == 1.4.0+cpu
* [ExPASy](https://web.expasy.org/translate/)

### Usage

Users can use their own data to train prediction models. If you want to train PHISGAE with your own dataset, you need to go through the following two processes: **1 Prepare the dataset** and **2 Predicting phage-host interactions**. In the meantime, we have provided the data that has been processed and you can implement PHISGAE directly using procedure: **2 Predicting phage-host interactions**.
 
### 1 Prepare the dataset

```
git clone https://github.com/JennyX212/PHISGAE
cd PHISGAE/code
python dna_f.py   
python pro_f.py   ####compute features derived from DNA and protein sequences
```

1. Download all phage DNA and protein sequences from [NCBI batchentrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez?), named **phage.fasta** and **phage.gb** into the *data* folder.
2. Run `python dna_f.py` and `python pro_f.py` to compute the features derived from DNA and protein sequences. **Note** If your protein sequences in one **.gb** file, you need to run `python protein_processing.py` to split them into different **.txt** files. For phages missing protein sequence information, you can use the [ExPASy](https://web.expasy.org/translate/) tool to translate their gene sequence into protein sequence.
3. Use `python sim.py` and module *## Other similarity calculation method ##* in `data_processing.R` to calculate the different similarity matrix between phages with `data/dna_pro_f/phage_dna_pro.csv`.
4. Use module *KNN graph* in `data_processing.R` to obtain refined phage-phage connection network.
5. Use `python hostI.py` to generate a represention of host-host connection network.
6. Convert the known phage-host pairs into different interaction matrices by running `python RAWA.py`


### 2 Predicting phage-host interactions

* tensorflow == 1.15.0

```
git clone https://github.com/JennyX212/PHISGAE
cd PHISGAE/code
python main.py  --input = "input_path" --output = "output_path"	 #### train model and get the prediction result
```

You just need to run `python main.py` and type `input =../data/train/species/` to train PHISGAE and get the result of PHISGAE for predicting phage-host interactions at species level. You can save the result by typing `output ="output_path"`, where "output_path" is the file path of the outcome. Besides, if you want the prediction results at different taxonomic level, you can directly modify the input path when typing `input = ` and change "species" to the other taxonomic level.

If you want to run PHISGAE with your own dataset, simply enter the path to the data set in the input prompt.
### Contact

Please feel free to contact us if you need raw data or any other help.

