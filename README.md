# AccuSleep

## DESCRIPTION

AccuSleep is a set of graphical user interfaces for scoring rodent
sleep using EEG and EMG recordings.

## Installation instructions:

Save the MATLAB package somewhere on your computer, then add it
to your MATLAB path. You can do this in the MATLAB "Current Folder"
window by right-clicking the AccuSleep folder, clicking "Add to Path"
--> "Selected Folders and Subfolders", then running the command
`savepath`
in the Command Window.

To get started, run `AccuSleep_GUI` and click the Help button, or run
`doc AccuSleep_instructions`
for a full explanation of these functions and the types of input
they require.

**`AccuSleep_GUI`** provides a convenient interface with all of the functions
you need, but if you want to batch process many recordings, you can
call the required functions yourself.

## Functions

- **`AccuSleep_GUI`**. A user interface for labeling sleep states, either
    manually or automatically
- **`AccuSleep_viewer`**. A user interface for manually labeling sleep states
- **`AccuSleep_classify`**. Automatically labels sleep states using a
    pre-trained neural network
- **`AccuSleep_train`**. Trains a neural network for labeling sleep states
- **`createCalibrationData`**. Generates a file that is required for automatic
    sleep state labeling for a given combination of subject and
    recording equipment

## Requirements:
- MATLAB version 2016b or later
- Statistics and Machine Learning Toolbox
