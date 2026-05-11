# AgingMicrobiome
This is the CHARGE collaboration to detect an aging signature in the human stool microbiome. 

*Information for the microbiombiom R-package are in the microbiombiom subdirectory*

American cohorts can use the Qiita environment to perform the preprocessing steps. Preprocessed data is expected to be called `feature.table.gg2-2022.10.qza` and the corresponding taxonomy file `df.gg2.taxonomy.qza`.

For European cohorts, this is due to the GDPR forbidden; hence, preprocessing should be performed locally. 

## Stool microbiome data harmonization


Data harmonization will take place using GreenGenes2. For both 16S- and shotgun-data this is done in the Qiime2 environment. We first need to install that.

```{bash}
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
PATH=$PATH:~/miniconda3/bin/
conda update conda
conda install wget
wget https://data.qiime2.org/distro/core/qiime2-2023.7-py38-linux-conda.yml
conda env create -n qiime2-2023.7 --file qiime2-2023.7-py38-linux-conda.yml
wget http://ftp.microbio.me/greengenes_release/2022.10/2022.10.taxonomy.asv.nwk.qza
conda init bash
source activate qiime2-2023.7
pip install q2-greengenes2
```

### Preprocessing 16S-data
Make a file `manifest` with a column `sample-id` and a column `absolute-filepath` with the sample ids and the paths to the **QC'd forward reads** of each sample. 
Use the `16S/PreProcessing16S` directory. 

Run `Preprocessing16S.sh` 

### Preprocessing Shotgun-data
Make a file `manifest` with a column `sample-id` and a column `absolute-filepath` with the sample ids and the paths to each sample's **QC'd forward reads**. 

Use the `shotgun/PreProcessingShotgun` directory. 

Install Woltka to be able to run the Web of Life
```{bash}
pip install woltka
```
You also have to download the Web of Life: [http://ftp.microbio.me/pub/wol2/genomes/] (file named all.fna). Place this file into a directory wol2 within the directory where your files are located. 

This can be done by the following code

```{bash}
mkdir -p wol2
cd wol2
curl -L -O https://ftp.microbio.me/pub/wol2/genomes/all.fna.xz
```

Run `PreprocessingShotgun.sh` 

### GreenGenes2
Following the GreenGenes2 harmonization, there are three options available. Please select the one that best suits your data.
- 16S V4 data: 
	Download, place in the same folder as preprocessed data, and run `closed_reference16SV4.sbatch` and `taxonomic_table16S.sbatch`
- 16S data not V4: 
	Get full backbone of Web of Life
	 ```{bash}
	wget http://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.fna.qza
	``` 
	Download `get_repset.py`, place in the same folder as preprocessed data, and run `closed_reference16SNonV4.sbatch` and `taxonomic_table16S.sbatch`
- Shotgun data: 
	Download, place in the same folder as preprocessed data, and run `closed_referenceShotgun.sbatch` and `taxonomic_tableShotgun.sbatch`

## Associations with aging phenotypes
Now using the previously created GreenGenes2 called data `feature.table.gg2-2022.10.qza` and taxonomy file `df.gg2.taxonomy.qza` we will run the analyses with the aging phenotypes.

### Metadata
The Python-script is under the assumption that metadata is a tab-separate .txt document. The first column needs to be sampleid, it is very important that the values in `sampleid` column overlap with the sampleid in the feature-table .qza file.

Run the following commands to make sure all the necessary packages are installed in the qiime2-2023.7 conda environment 
```
source activate qiime2-2023.7
conda install -c conda-forge click scikit-bio pandas numpy biom-format statsmodels lifelines patsy
```

Change the file paths in changenames.sbatch to start the analyses; you can here also add study-specific covariates that should be considered a factor in the factor line. For instance, if you need the study site to be treated as a factor add it here.

Make sure variables are coded according to the analysis plan:
* sex: "men" and "women";
* ppump: proton-pump inhibitor usage, if participant indicates use 1, not 0;
* metfor: metformin usage, if participant indicates use 1, not 0;
* statin: stating usage, if participant indicates use 1, not 0;
* bmi: numeric, as many digits as available;
* race: as factor, "white" will be reference
* dietscore: numeric, as many digits as available (Nettleton: DOI: 10.1093/aje/kws297)
* fiber: numeric, as many digits as available

### Models
Change models in `mods_agingmicrobiome.txt` to let them include study-specific variables such as Site.

### Outcomes
Remove outcomes your cohort does not participate in in `outs_agingmicrobiome.txt` 

To clarify:
* age: chronological age
* fi: frailty index
* mortality: all-cause mortality; If mortality is the outcome of interest the script also expects the variables `age` (chronological age) and `studytime` (time between sampling and censoring) to be present in the meta-data

### Subset
Remove from the `subsets_agingmicrobiome.txt` file subsets your cohort does not take part in.

To clarify:
* all: all participants
* men: men
* women: women
* age_1: participants with a chronological age >= 18 & chronological age < 40
* age_2: participants with a chronological age >= 40 & chronological age < 50
* age_3: participants with a chronological age >= 50 & chronological age < 60
* age_4: participants with a chronological age >= 60 & chronological age < 70
* age_5: participants with a chronological age >= 70

### IMPORTANT
Please change your cohort name in `cohort.txt` 

