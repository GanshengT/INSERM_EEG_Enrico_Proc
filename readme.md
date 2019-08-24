This is the repository for processing EEG data using MNE python package <br>
Readme file is under construction

<a href="https://www.inserm.fr/">
    <img src="https://upload.wikimedia.org/wikipedia/fr/c/cd/Inserm.svg" alt="INSERM logo" title="INSERM" align="right" height="60" />
</a>

Spectral Analysis and preprocessing of EEG signals
======================
As both sample sizes and EEG channel densities increase, traditional processing approaches like manual data rejection are becoming unsustainable. Different from the classical EEG signals analysis, here we provide a convenient wrapper of [mne](https://martinos.org/mne/stable/index.html) functions to do the semi-automatic preprocessing and spectral analysis. The examples and scripts are correspond to this [specific experiment design](https://www.researchgate.net/publication/326542100_Differential_effects_of_non-dual_and_focused_attention_meditations_on_the_formation_of_automatic_perceptual_habits_in_expert_practitioners).



## Table of content

- [Prerequest](#prerequest)
    - [Environment](#environment)
    - [Database](#database)

- [Workflow](#workflow)
    - [Preprocessing](#preprocesing)
    - [Statistical analysis](#composer)

- [Scripts explanation and examples](#scripts explanation and examples)
    - [ready-to-run](#ready-to-run)
    - [Go to the import view](#go-to-the-import-view)
    - [Import the uploaded page tree file](#import-the-uploaded-page-tree-file)
-[Troubleshooting]
- [Authors](#authors)
- [License](#license)
- [Extentions](#links)

## Prerequest
### Environment
Language we used for programming is python, python 3.5 or 3.6 is required, MNE-Python is the essential tool to preprocess the data and perform sensor-level analysis.  
[How to install mne]  
To perform mathematical easily, we use numpy and scipy

### Database
The database are the EEG recording files using eeglab format. Each file stands for one state, three in total. (Please refer to the experiment design). Additionally, each subject performs two sessions, thus, each subjects has 6 raw-data files. The recording contains 64 EEG channels and 2 EOG channels (VEOG and HEOG)


## Workflow
The workflow consists of two parts. 
### Preprocessing
Electroenchephalography (EEG) recordings have a high degree of artifact contamination, so the first part is artifact rejection which includes following steps:
* visual inspection: This steps helps to define the bad electrodes and do not take them into conside
ration in the following steps, bad eletrodes will be interpolated after preprocessing.
* filter - 1-100Hz passband filter and 50 Hz notch filter (France) *Notice: this happens after extracting signal in practice but it might have edge effet, for practical purpose, we fix sampling rate as 512 for all the recording files*
* signal extraction and annotation engineering:</br>
 The raw data is cut from event 254 to event 255, then we recode the events as follows:  
    * events code: state + condition + session:
        1. state: 1:VD 2:FA 3:OP  
        2. condition: 1:baseline 2:safe 3:threat  
        3. session: 1:session1 2:session2  
    At the end, we concatenate three recording files for one sessions into one epoching file that we call *full_epochs*.
* ASR - [artifact subspace reconstruction](https://www.ncbi.nlm.nih.gov/pubmed/30440615)  
    * selecting baseline signal and using yule walker to amplify artifact components
    * cut blinks based on VEOG channel and remove periods with abnormally high-power content from continuous data
    * calculate the rejecting threshold based on cleaned baseline signal (please refer tp ASR documentation)
    * apply yule walker on *full_epochs* to get *epochs4detect*, reconstruct *full_epochs* window by window with respect to the correlation between *epochs4detect* and *full_epochs*.
    
* *full_epochs* concatenation and apply [ICA](https://www.sciencedirect.com/science/article/pii/S0893608000000265) on *full_epochs* to exclude the residus of artifact components (especially blinks and saccades that are not focus by ASR because the yule walker that we used aims to amplify high-frequency artifacts)

* Visually exclude ICA components and run [Autoreject](https://autoreject.github.io/index.html) with local threshold initially, while the rejecting rate is higher than 10%, we choose a more tolerant global threshold. The datas after above procedures are called *precleaned_epochs*, one for each subject.  


### Statistical analysis
As for the statistical analysis, afk processing, we use R language to fit mix-effeted model to perform tests which includes:**to be verified by Arnaud**
* extract psd and band power from *precleaned_epochs*. ->panda dataframe or csv
* visualize data distribution and normalise data and eliminate outliners
* define variable and test the difference within groups, states etc, and the correlation between variables.
* post-hoc tests  
we use mne-python to perform spatio-clustering permutation test
* information extraction:
    * one matrix for paired test, two matrix for unpaired test.

### Scripts explanation and examples
This chapter explains most of the methods in folder script_tan, those not being described are in developpement version or of small importance.
## ready-to-run (in calculation machine)
* methods .py in preprocessing
    * utils_preProcessingWorkflowJuly05.py
    this file consists of basic wrap-up function of mne-python
        * *autorej_rate* takes epochs after autoreject as argument, and return the percentage of the epoch that have been rejected by autoreject.
        * *get_epochs_ASR_clean* takes *subject id* and *session number* as arguments, and it returns cleaned epochs after artifact subspace reconstruction. To run this method, one has to well define:
        ```python
 raw_data_path = 'your/path'
 montage_fname = 'your/path'
 preProc_ica_path = 'path for storing ica mixing matrix in the format of .fif'
 report_path = 'your/path'
```


## Authors
* [**Gansheng Tan**](https://ganshengt.github.io/) - *Initial work* 
* **Arnaud Poublan-couzardot** - *Statistical analysis* 



## License
The scripts are open-source but intended to be modified by the member of Dycog Team. The Dataset (suitable policy)

## Extentions

* [Compared to HAPPE](https://www.frontiersin.org/articles/10.3389/fnins.2018.00097/full)
* [automatic ICA components rejection](https://www.ncbi.nlm.nih.gov/pubmed/21810266)

## Acknowledgments

* Advice given by Manu MABY and Françoise LECAIGNARD has been a great help in
* I am particularly grateful for the assistance given Arnaud Poublan-couzardot and the comments by Antoine LUTZ
