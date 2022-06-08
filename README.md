# gCaMP-Analysis
*Analysis code for processing in vivo calcium-imaging data collected in the Haptics Lab.*
*A collection of analysis Scripts, Functions, and GUIs primarilly written in MATLAB.*
*-mduhain*
..
..

## CONTENTS
### 1. Main_Analysis.m 
  - The old version of analysis code, a MATLAB script to quantify neuron activity over time, requires manual identification of neuron locations.

### 2. quick_2p_analysis.mlapp
  - The new GUI version of analysis, designed to be used during experiments to provide rapid ientification of neuron tuning preferences.
  - Acceps in '.tif' image stacks, splits into RFP/GFP channels, and attempts to automatically identify RFP-tagged neuron locations.
  - Also accepts in '.mat' output data from (https://github.com/mduhain/Haptics-Lab-GUI/Stimulator-Gui.mlapp).
  - Parses through stim times and frame times, to identify neuron tuning from the gCaMP signal.

### 3. imneuron.m
  - A new, automatic cell tracking function written in MATLAB, accepts '.tif' image stacks and tracks RFP-tagged cell bodies over time
  - More efficient and accurate than the cell tracking in 'quick_2p_analysis.mlapp', will be integrated in the future.



##OTHER USEFUL REPOS
1. imageAnalysis_gui-Master by @cdeister (https://github.com/cdeister/imageAnalysis_gui)
2. 
