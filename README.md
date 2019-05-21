##codes to measure hcal depth distribution of jets##

This code is tested on CMSSW_10_5_0 with added depth segmentation
Modify `jetemet_tree.py` for desired variables

###Produce flat Ntuples###
Make a list of files stored in files.txt
run: `python jetmet_analyzer inputFiles=files.txt`
This will read all files from `files.txt` and produce an flat root file called `jetmetNtuples.root`

###Make jet resolution plots###
run `python jet_resolution.py`
This will read from `jetmetNtuples.root` and make plots in the `result/` directory

###Make MET resolution plots###
For ZMM samples only!
run `python met_resolution.py`
This will read from `jetmetNtuples.root` and make plots in the `result/` directory

###Plot Hcal depth distributions for jets###
`python plot_depth.py`
Also read from `jetmetNtuples.root`

###Plot Hcal depth distributions for pf particles###
`python plot_pf_distributions_in_jets.py inputFiles=files.txt`
This need to read from MINIAOD files
