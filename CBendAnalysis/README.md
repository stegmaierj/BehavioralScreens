# CBendAnalysis

# Prerequisites
To be able to run the CBendAnalysis script, make sure to have one of the latest MATLAB versions installed including the Image Processing Toolbox. Note that compressed video files may not be processible if the codices are not installed properly on your system.

# Using the CBendAnalysis Scripts
The main script to start an analysis is named `BatchAnalysis.m`. If you run the script it will ask for an input folder and expects a folder containing `*.avi` video files.

To speedup the processing if only quantitative information is required, you may want to disable the writing of video files in the main script `BatchAnalysis.m` by setting in line 36 to `generateResultVideo = false;`.
