# AgingMicrobiome
This is the CHARGE collaboration to detect an aging signature in the human stool microbiome. 

American cohorts can use the Qiita environment to perform the preprocessing steps. Preprocessed data is expected to be called

For European cohorts, this is due to the GDPR forbidden; hence, preprocessing should be performed locally. 

## Stool microbiome data harmonization


Data harmonization will take place using GreenGenes2. For both 16S- and shotgun-data this is done in the Qiime2 environment. We first need to install that.

```{bash}
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
conda update conda
conda install wget
wget https://data.qiime2.org/distro/core/qiime2-2023.7-py38-linux-conda.yml
conda env create -n qiime2-2023.7 --file qiime2-2023.7-py38-linux-conda.yml
wget http://ftp.microbio.me/greengenes_release/2022.10/2022.10.taxonomy.asv.nwk.qza
conda init bash
pip install q2-greengenes2
```

### Preprocessing 16S-data
Make a file `manifest` with a column `sample-id` and a column `absolute-filepath` with the sample ids and the paths to the **first reads** of each sample. 
Download the `PreProcessing16S.zip` file. 

Run `Preprocessing16S.sh` 

### Preprocessing Shotgun-data
Make a file `manifest` with a column `sample-id` and a column `absolute-filepath` with the sample ids and the paths to each sample's **alligned reads**. 

Download the `PreProcessingShotgun.zip` file. 

Install Woltka to be able to run the Web of Life
```{bash}
pip install woltka
```
You also have to download the Web of Life: [http://ftp.microbio.me/pub/wol2/genomes/] (file named all.fna). Place this file into a directory wol2 within the directory where your files are located. 

Run `PreprocessingShotgun.sh` 

### GreenGenes2
According to the GreenGenes2 harmonization, there are 3 options; choose the one suiting for your data:
- 16S V4 data: Download, place in the same folder as preprocessed data, and run `closed_reference16SV4.sbatch` and `taxonomic_table16SV4.sbatch`
- 16S data not V4: Download `get_repset.py`, place in the same folder as preprocessed data, and run `closed_reference16SNonV4.sbatch` and `taxonomic_table16SNonV4.sbatch`
- Shotgun data: Download, place in the same folder as preprocessed data, and run `closed_referenceShotgun.sbatch` and `taxonomic_tableShotgun.sbatch` 
