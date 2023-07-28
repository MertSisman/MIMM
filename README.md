# MIMM
Microstructure Informed Myelin Mapping

This code is the implementation of the Microstructure Informed Myelin Mapping (MIMM) algorithm. It provides the myelin volume fraction (MVF) maps from multi gradient echo (mGRE) MRI data.

The code requires mGRE magnitude data, QSM map, and weighting factor lambda_chi for the basic version. For orientation informed approach, the FA map and fiber orientation (theta) map also needs to be provided.

An example dataset, corresponding output maps, and dictionaries can be downloaded [here](https://zenodo.org/record/8193673).

If you want to process your own data, the suggested toolboxes for QSM reconstruction is [MEDI Toolbox](https://pre.weill.cornell.edu/mri/pages/qsm.html), and for DTI reconstruction [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSL).
