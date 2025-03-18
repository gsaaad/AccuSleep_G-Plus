# AER-Accusleep

The AER Lab has developed a customized version of Accusleep for efficient sleep-scoring.

## Updates - **03/12/2025**

### Added Features:

- Introduced lab-standard preferences for the number of bins/epochs to display [15 epochs].
- Enhanced efficiency and consistency in scoring by computing and plotting the previous and forward 1-minute EEG and EMG.
- Computed and plotted the PSD for up to 9 epochs, displaying essential frequencies (0-20 Hz).
- Color-coded PSD plots in bins (0-4, 5-9, 10-15, 16-20 Hz).
- Enabled auto-scroll by default.
- Implemented up to 6 scoring labels (REM, WAKE, NREM, REM*, WAKE*, NREM\*) with color coding for easy identification.
- Redesigned the GUI for improved usability.
- Opened Viewer on full screen – this is how it should be used.
- Added auto-saving functionality for every 450 keystrokes (~30 minutes of data).
- Compute average Delta, Theta and Theta Ratio, displaying it on the GUI.

## Acknowledgements

We extend our gratitude to [George Saad](https://github.com/gsaaad) for customizing lab-specific preferences and developing a tool for widespread use in the lab.

Special thanks to [Zeke Barger](https://github.com/zekebarger) for providing the original Accusleep final product, and to [Franz Weber](https://www.med.upenn.edu/weberlab/) for creating an early version of the manual labeling interface.

## Screenshots

### Primary Interface (AccuSleep_GUI)

![Primary Interface](assets/Accusleep_Gui.jpg)

### Interface for Manual Sleep Scoring (AccuSleep_viewer)

![Manual Sleep Scoring Interface](assets/GS_Accusleep_Viewer.jpg)
