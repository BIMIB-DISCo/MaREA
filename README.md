# Load Dataset (browsing the file system)

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/dataset_loaded.png?raw=true "Title")

Use **Select dataset** to browse your files and select a dataset in csv or tsv format. This file have to contain one column with gene ID write in one of the following format:
- Hugo Symbol (*IDH1*)
- Hugo ID (*HGNC:5382*)
- Ensemble ID (*ENSG00000138413*)
- Entrez ID (*3417*)

This column have to be specified with the **Gene IDs col** drop-down listist.

**Gene Identifier** allow to specify which kind of gene identifier is used in the dataset.

The **First patient col** drop-down list is used to specify the first patient (or sample) in the dataset. Following columns will be considered others patients as well till the end of file.


# Process Dataset
Specify type of ID for genes and name of column for IDs and first patients
![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/addDataset.png?raw=true "Title")

Push **Add dataset** to import the selected dataset in the app.

Little summary of dataset imported is created during the process and is displayed in *Datasets* text box.

# Load and Process second (optional) dataset

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/second_dataset.png?raw=true "Title")

With this process any number of datasets could be imported.

# Load Model
![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/metabolicModel.png?raw=true "Title")

A metabolic model *COBRA toolbox* compliant have to be chosen to perform the further analysis. MaREA already provides *Recon 2.2* a genome wide model (ref) and *HMRcore* a core model manually curated by our research group. 
With *Custom* another model can be used, but user have to specify in which format model.genes are written. (This Gene ID can be the same of the dataset or different, in the latter case MaREA automatically convert model genes ID in dataset ID.

In **Model** text box a summary of the metabolic model chosen will appear. 

# Compute RAS

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/computeRAS.png?raw=true "Title")


![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/advancedOptions.png?raw=true "Title")


![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/RASready.png?raw=true "Title")


Select at least one entry from *Dataset* list box and MaREA will compute RAS (reaction activity score) for each reactions in metabolic model and for each patient or sample in dataset. 

**Advanced Option** 
 - Turn on **Parallele mode** will enable faster computation if Parallele Toolbox is installed.
 
 *Following option are not available yet*
 - **Resolve (A and *nan*) rules**: Rules could consider two genes in *and*, one of them don't have information about it's transcript level (*nan* value).
	- with *A* selected MaREA will not consider *nan* gene so the corresponding result will be just transcript level of gene **A**
	- with *nan* selected the result of the rule will be *nan*
	**N.B.** if transcript level of **A** is 0 the result will always be 0
	
 - **Convert log scaled to decimal values**
	If *yes* is selected MaREA will consider dataset values (t) log2 scaled so before perform analysis 2^t will apply.

 - **Threshold 0 transcript**
	Define a threshold: each value below it will replaced by 0


# Compare datasets
![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/compareDatasets.png?raw=true "Title")

**Export Data** button allow to save dataset with previously performed computation in a *.mat* file. (Select at least one entry from **Datasets with RAS** list box.


***Enrichment* section
First of all select at least two entry from *Dataset with RAS* list box.

Select options from drop-down list
 - **Mode**
	- *Pairwise*: each dataset will compared again each other.
	- *One vs other*: each dataset will compared again mean values of others datasets.

 - **Hypotesis test** chose one test to compute p-value for difference between transcript distributions. Test available are: *Kolmogorov-Smirnov*, *Kruskal-Wallis*, *t-student*, *One way anova* or *Do not compute* for skip this analysis.
 
 - **Fold Change** chose if or not to compute the ratio of each reaction value between dataset.
 
 Use **Compare dataset** to create .txt file with the statistical analysis output.
	
# Color map

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/colorMap.png?raw=true "Title")
![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/map_enriched.png?raw=true "Title")

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/mapDownloaded.png?raw=true "Title")

If HMRcore was chosen as metabolic model MaREA allows to color a metabolic map to better visualize enrichment results. 

Chose a *pV - Thr*: reactions with pValue under this threshold will draw with gray dashed line. If no test was used skip this point.

Chose a *FC - Thr*: reactions with Fold Change under this threshold will gray colored. This value must be a coefficient between 0 and 1. (0.2 = 20%)

# Get clusters

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/getClusters.png?raw=true "Title")
![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/clustersDownloaded.png?raw=true "Title")

MaREA allows to stratify dataset(s) with K means cluster method.

Use **Number of clusters** to chose a range of K values. Each value correspond to a  number of groups dataset will be divide. So chose a range from 2 to 5 will divide same dataset in 2,3,4 and 5 groups. 

Use **scale variables** to normalize dataset values.
 - *divide row by max*: divide each row by their maximum values.
 - *divide col by max*: divide each col by their maximum values.
 - *none*: do nothing
 
Use **Data** to chose if stratify dataset based on metabolic transcript or RAS score

**Number of replicates**: perform analysis many times and return the best clusterization.

**Clusterize** start the analysis, many files will be created, we suggest to choose an empty directory.


