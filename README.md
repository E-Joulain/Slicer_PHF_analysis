# Slicer_PHF_analysis
This extension is designed for the analysis of proximal humerus fractures (PHF) using 3D Slicer. It provides tools to place anatomical landmarks, compute key parameters, and simulate virtual reconstructions.

## Installation

If not already installed, add the StudyPHF extension to your 3D Slicer environment via the Extension Manager. Fisrt, download the folder "PHF_analysis" which contains the module. Then from the extension manager click on "Install from file...", select the downloaded folder. 

## Getting Started
### Step 1: Load DICOM Data

Download the DICOM images you plan to analyze. This may include:\
  •	Pre-operative (preop) images \
  •	Post-operative (postop) images \
  •	Single or multiple patients \
Import the DICOM data into 3D Slicer. Be careful to add/annotate an indicative name for each volume and patient. It will be helpful when choosing the input volume for the analysis.

### Step 2: Set Up Patient Information
Enter a patient name or ID that will be used when saving the data. It helps you identify the analysis results, especially if you are analyzing multiple patients.

Select the X-ray volume (called a volume in Slicer) from the dropdown menu labeled in red.

### Step 3: Define Image Side
Specify whether the shoulder in the image is left or right.

### Step 4: Begin the Analysis
Choose the appropriate workflow: \
o	Use the green buttons for preoperative analysis \
o	Use the blue buttons for postoperative analysis \
o	Or one after the other if you are performing both analysis. You will need to change the selected volume when changing from pre to post (the red button).

#### Landmark Placement
1.	Click Place Fiducials to begin landmark annotation. \
2.	Place all 13 anatomical fiducial points as described in the Analysis Table. \
3.	Once done, click Compute Parameters to calculate key anatomical metrics. \
To visualize the computed data, click the Results button (or similar if named differently in the UI). \

#### Managing Fiducials
To re-place landmarks, click Reset Fiducials, then Place Fiducials again. When satisfied with your placement, click Export Fiducials to save the point coordinates.
All saved data (coordinates, results, etc.) will be stored in the extension's Data folder. The full path will appear in the Python console


## Optional: Virtual Humerus Analysis
You can simulate a virtual replacement of the humeral head by clicking Virtual Humerus Analysis, located below the blue  postop buttons.
