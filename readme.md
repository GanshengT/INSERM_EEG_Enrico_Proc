This is the repository for processing EEG data using MNE python package <br>
Readme file is under construction

<a href="https://www.inserm.fr/">
    <img src="https://upload.wikimedia.org/wikipedia/fr/c/cd/Inserm.svg" alt="INSERM logo" title="INSERM" align="right" height="60" />
</a>

Spectral Analysis and preprocessing of EEG signals
======================
As both sample sizes and EEG channel densities increase, traditional processing approaches like manual data rejection are becoming unsustainable. Different from the classical EEG signals analysis, here we provide a convenient wrapper of [mne](https://martinos.org/mne/stable/index.html) functions to do the semi-automatic preprocessing and spectral analysis. The examples and scripts are correspond to this [specific experiment design](https://www.researchgate.net/publication/326542100_Differential_effects_of_non-dual_and_focused_attention_meditations_on_the_formation_of_automatic_perceptual_habits_in_expert_practitioners).



## Table of content

- [Workflow](#workflow)
    - [Preprocessing](#typo3-extension-repository)
    - [Statistical analysis](#composer)
- [Prerequest](#typo3-setup)
    - [Environment](#extension)
    - [Database](#database)
- [Scripts explanation and examples](#page-setup)
    - [Upload the page tree file](#upload-the-page-tree-file)
    - [Go to the import view](#go-to-the-import-view)
    - [Import the uploaded page tree file](#import-the-uploaded-page-tree-file)
-[Troubleshooting]
- [Authors](#authors)
- [License](#license)
- [Extentions](#links)

## Workflow
The workflow consists of two parts. 
### Preprocessing
Electroenchephalography (EEG) recordings have a high degree of artifact contamination, so the first part is artifact rejection which includes following steps:
* visual inspection: This steps helps to define the bad electrodes and do not take them into conside
ration in the following steps, bad eletrodes will be interpolated after preprocessing.
* filter - 1-100Hz passband filter and 50 Hz notch filter (France) *Notice: this happens after extracting signal in practice but it might have edge effet*
* signal extraction and annotation engineering:</br>
 The raw data is cut from event 254 to event 255, then we recode the events as follows:
    *events code: state + condition + session:
        1. state: 1:VD 2:FA 3:OP\n",
        2. condition: 1:baseline 2:safe 3:threat\n",
        3. session: 1:session1 2:session2\n",





### Statistical analysis



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
