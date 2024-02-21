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
Download the `PreProcessing16S.zip` file. 

Run `Preprocessing16S.sh` 

### Preprocessing Shotgun-data
Make a file `manifest` with a column `sample-id` and a column `absolute-filepath` with the sample ids and the paths to each sample's **QC'd forward reads**. 

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

## Associations with aging phenotypes
Now using the previously created GreenGenes2 called data `feature.table.gg2-2022.10.qza` and taxonomy file `df.gg2.taxonomy.qza` we will run the analyses with the aging phenotypes.

### Metadata
The R-script is under the assumption that metadata is a separate .txt document named metafile `metadata_agingmicrobiome.txt` 
Make sure variables are coded according to the analysis plan:
* sex: "men" and "women";
* ppump: proton-pump inhibitor usage, if participant indicates use 1, not 0;
* metfor: metformin usage, if participant indicates use 1, not 0;
* statin: stating usage, if participant indicates use 1, not 0;
* bmi: numeric, as many digits as available;
* race: as factor, "white" will be reference
* dietscore: numeric, as many digits as available
* fiber: numeric, as many digits as available

### Models
Change models in `mods_agingmicrobiome.txt` to let them include study-specific variables such as Site.

### Outcomes
Remove outcomes your cohort does not participate in in `outs_agingmicrobiome.txt` 

To clarify:
* age: chronological age
* fi: frailty index
* cont: continuous frailty phenotype
* mortality: all-cause mortality

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

### Send back
Please include a README Word file with the following study descriptive information:
General	
* Named individuals with contact details who should be considered in future publications, including: analysts and PIs (if possible indicate their role in your study, so that we can direct any queries to the appropriate person)
*	Brief description of your study.
*	Acknowledgements for your study, including funding sources

For analytic sample
*	Sample size
*	N of women (%)
*	N participants living in the same household as other study participants (%)
*	Mean age (years ± SD)
*	N PPI users (%)
*	N participants with missing information on PPI use (%)
*	N metformin users (%)
*	N participants with missing information on metformin use (%)
*	N lipid lowering medication users (%)
*	N participants with missing information on metformin use (%)
*	Mean BMI (± SD)
*	N participants with missing information on BMI
*	N per ethnic subgroup (%)
*	N participants with missing information on ethnicity
*	Mean fiber intake (± SD)
*	N participants with missing information on fiber intake
*	Mean diet score (± SD)
*	N participants with missing information on diet score

For microbiome data
*	Sample collection.
*	Microbial DNA isolation
*	DNA sequencing

For phenotypes
*	Chronological age:
  *	Minimal age
  *	Maximum age
  *	Median age
  *	IQR age
*	Frailty index
  *	List of variables included in frailty index
  *	Mean frailty index (± SD)
  *	Minimal frailty index
  *	Maximum frailty index
  *	Median frailty index
  *	IQR frailty index
*	Continuous frailty phenotype
  *	Mean continuous frailty phenotype (± SD)
  *	Minimal continuous frailty phenotype
  *	Maximum continuous frailty phenotype
  *	Median continuous frailty phenotype
  *	IQR continuous frailty phenotype 
*	All-cause mortality
  *	Number of deaths (%)
  *	Number of people lost to follow-up (%)
  *	Mean follow-up time (± SD)
  *	Minimal follow-up time
  *	Maximum follow-up time
  *	Median follow-up time
  *	IQR follow-up time

**Send the document with study participants with descriptives and your used scripts back to l.m.kuiper[at]erasmusmc.nl**

