# gCaMP-Analysis
*Analysis code for processing in vivo calcium-imaging data collected in the Haptics Lab.
A collection of analysis Scripts, Functions, and GUIs primarilly written in MATLAB.
  -mduhain*

## CONTENTS

### 1. 'dataParser.m'
  - The main script to preform the initial creation of an im2p class for each experimmental session.
  - Combines session metadata with gfp stacks, rfp stacks.
  - designed to be called at root of nested folders organized: ~/m000(mouseNum)/YYYY-MM-DD/Exp000/TactileStim_00000.tif 
  - TODO: Add functionality for vasculature images.
  
### 2. 'im2p.m'
  - The custom class to hold all data relevant to an experiment session.
  
### ~ ~ ~ ~ ~ ~

### 8. 'Main_Analysis.m'
  - The old version of analysis code, a MATLAB script to quantify neuron activity over time, requires manual identification of neuron locations.

### 9. 'quick_2p_analysis.mlapp'
  - The new GUI version of analysis, designed to be used during experiments to provide rapid ientification of neuron tuning preferences.
  - Acceps in '.tif' image stacks, splits into RFP/GFP channels, and attempts to automatically identify RFP-tagged neuron locations.
  - Also accepts in '.mat' output data from (https://github.com/mduhain/Haptics-Lab-GUI/Stimulator-Gui.mlapp).
  - Parses through stim times and frame times, to identify neuron tuning from the gCaMP signal.

### 10. 'imneuron.m'
  - A new, automatic cell tracking function written in MATLAB, accepts '.tif' image stacks and tracks RFP-tagged cell bodies over time
  - More efficient and accurate than the cell tracking in 'quick_2p_analysis.mlapp', will be integrated in the future.



## OTHER USEFUL REPOS
1. imageAnalysis_gui-Master by @cdeister (https://github.com/cdeister/imageAnalysis_gui)
2. 
