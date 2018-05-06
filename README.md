# Requirements	
MaREA requires MATLAB and the COBRA Toolbox.

# Run MaREA
To run this MATLAB tool simply add the <*MaREA-master*> directory in your Matlab path and digit MaREA in MATLAB Command Window.

# Load Dataset (browsing the file system)

Use **Select dataset** (tab Dataset) to browse your files and select a dataset in csv or tsv format. This file must contain one column with the genes ID in one of the following formats:
- Hugo Symbol (*IDH1*)
- Hugo ID (*HGNC:5382*)
- Ensemble ID (*ENSG00000138413*)
- Entrez ID (*3417*)

Once the dataset is correctly loaded the message "dataset loaded" appears (see screenshot below).

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/dataset_loaded.png?raw=true "Title")


# Process Dataset
Once the dataset is loaded, the user must specify the name of the column reporting the genes ID from the **Gene IDs col** drop-down list.

**Gene Identifier** allows to specify which kind of gene identifier is used in the dataset.

The **First patient col** drop-down list is used to specify the first patient (or sample) in the dataset. Following columns till the end of file will be considered as patients/samples.


Push **Add dataset** to import the selected dataset in the app.
A Summary of  imported datasets will be created and displayed in the *Datasets* text box (see screensnot below).

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/addDataset.png?raw=true "Title")


# Load and Process a second (optional) dataset

By repeating the process above any number of datasets can be imported.

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/second_dataset.png?raw=true "Title")


# Load Model
A ( *COBRA toolbox* compliant - extension .mat or .xml) metabolic model must be chosen to perform further analyses. MaREA already provides the *Recon 2.2*  genome wide model and the *HMRcore* core model (manually curated by our group). 
By selecting *Custom*, a different model can be imported. In this case, the user must specify which identificative for genes the model vector model.genes uses. This Gene ID can either be the same used for the dataset or a different one, in the latter case MaREA automatically converts model genes IDs into dataset IDs.

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/metabolicModel.png?raw=true "Title")

Once the model is loaded, a summary of the metabolic model will appear in **Model**. 

# Compute RAS
Select at least one entry from *Dataset* list box and MaREA will compute the RAS (Reaction Activity Score) for each reaction in the metabolic model and for each patient or sample in the dataset. 

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/computeRAS.png?raw=true "Title")

Wait until RASs are computed.

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/RASready.png?raw=true "Title")


Slide the botton **Advanced Option** to display more options (see screenshot below).

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/prepareData.png?raw=true "Title")

 - Turn **Parallele mode** on to enable faster computation, provided that Parallele Toolbox is installed.
 
 - **Resolve (A and *nan*) rules**: Rules could consider two genes in *and*, one of them don't have information about it's transcript level (*nan* value).
	- with *A* selected MaREA will not consider *nan* gene so the corresponding result will be just transcript level of gene **A**
	- with *nan* selected the result of the rule will be *nan*
	**N.B.** if transcript level of **A** is 0 the result will always be 0
	
 - **Convert log scaled to decimal values**
	If *yes* is selected MaREA will consider dataset values (t) log2 scaled so before perform analysis 2^t will apply.

 - **Threshold 0 transcript** (feature still to be implemented)
	Define a threshold: each value below it will be replaced by 0


# Compare datasets
![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/compareDatasets.png?raw=true "Title")

**Export Data** button allows to save the dataset with previously performed computation in a *.mat* file. Select at least one entry from the **Datasets with RAS** list box.


***Enrichment* section
First of all select at least two entries from the *Dataset with RAS* list box.

Choose the following options from drop-down list:
 - **Mode**
	- *Pairwise*: each dataset will be compared again each other.
	- *One vs other*: each dataset will compared again mean values of other datasets.

 - **Hypotesis test** choose one test to compute p-value for the difference between transcript distributions. Test available are: *Kolmogorov-Smirnov*, *Kruskal-Wallis*, *t-student*, *One way anova* or *Do not compute* for skip this analysis.
 
 - **Fold Change** choose wether to compute the ratio of each reaction value between datasets.
 
 Use **Compare dataset** to create a .txt file with the statistical analysis output.
	
# Color map

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/colorMap.png?raw=true "Title")
![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/map_enriched.png?raw=true "Title")

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/mapDownloaded.png?raw=true "Title")

Provided that HMRcore was chosen as metabolic model, MaREA allows to color a metabolic map to better visualize enrichment results. 

Chose a *pV - Thr*: reactions with pValue under this threshold will draw with gray dashed line. If no test was used skip this point.

Chose a *FC - Thr*: reactions with Fold Change under this threshold will gray colored. This value must be a coefficient between 0 and 1. (0.2 = 20%)

# Get clusters

![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/getClusters.png?raw=true "Title")
![Alt text](https://raw.githubusercontent.com/BIMIB-DISCo/MaREA/master/Images/clustersDownloaded.png?raw=true "Title")

MaREA allows to stratify dataset(s) with K-means cluster method.

Use **Number of clusters** to choose a range of K values. Each value corresponds to a number of clustes in which the dataset will be split. Chosing a range from 2 to 5 will divide the same dataset in 2,3,4 and 5 groups. 

Use **scale variables** to normalize dataset values.
 - *divide row by max*: divide each row by their maximum values.
 - *divide col by max*: divide each col by their maximum values.
 - *none*: do nothing
 
Use **Data** to chose whether stratify the dataset based on metabolic transcript or on RAS.

**Number of replicates**: perform analysis many times and return the best clusterization.

**Clusterize** start the analysis, many files will be created, we suggest to choose an empty directory.


