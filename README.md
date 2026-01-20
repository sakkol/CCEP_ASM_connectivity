# README

# CCEP\_ASM\_connectivity



## General Info

This folder contains the data analysis and figure generation MATLAB functions and scripts for the Akkol et al. with the title "Anti-Seizure Medications Alter Functional and Effective Connectivity as Measured with Intracranial Electroencephalography". Scripts that will generate the figures and analyses results are named accordingly. The neural data and electrode localization must be obtained from the first or senior authors directly. Save the neural data and electrode localization into the "CCEP\_ASM\_connectivity" folder and keep the same folder structure. This is needed for the analyses to run accurately (also see Requirements below).



## Requirements:

* Main dependency of the analyses is the [Fieldtrip toolbox](https://www.fieldtriptoolbox.org/).
* For the folder structure and naming convention, [CoreAdminCodes](https://github.com/sakkol/CoreAdminCodes) repository is used.
* For the colormaps, [master\_ColorMaps](https://github.com/sakkol/master_ColorMaps) repository of the first author needs to be used.



## File information

We have included scripts and functions necessary to obtain the analyses and figures for the article. Two main scripts (starting with "Pipeline") include the analyses steps for both CCEPs and resting state functional connectivity. Scripts that were named "FigureXX" produce the important parts of the figures (then further processed in InkScape). All the functions are used by different scripts/functions.

