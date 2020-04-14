# BioXAS_Toolbox
Set of tools to convert/process Acquaman data files. Beamline-independent, works with files recorded at IDEAS, BioXAS, VESPERS. More templates can be added to tje acquamanDataFormats.py.
All created in python 3.7, not tested in other incarnations of python3, but might work right away.

1. Acquaman CDF converter (acquaman_exporter_sqlite3_DarkCorrection.py)
Simple script to export Acquaman raw data files to the HDF5 format. Dark currents corrected columns are added automatically as well as most of the metadata - motor positions, scan regions and so on.

Dependencies:
spacepy, h5py

Script uses the CDFlib wrapper included in spacepy, it's necessary to install both the spacepy package and the library dlls.
Recommended way to install spacepy is pip:

pip install spacepy

CDF library is available at https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf37_1/windows/
Just download and run the installer for your platform (cdf37_1_0-setup-64.exe for 64-bit Windows)

h5py is a python binding for HDF5, once again use pip:

pip install h5py

Usage:

Check out the inline comments in the script, it's simple. In principle, one should define the input and output folders and scan name (the same as automatically exported file).
Leaving the percent symbol in the scanNameFiler will make the script export all spectra (may take some time and gigabyte of disk space).
Specify the input directory with raw Acquaman data, usually it's something like AcquamanMainData\users\PROJECT\userData, check if it contains the userdata.db file.
If the output folder doesn't exist, the script will create it.


2. Acquaman HDF Explorer, or just ACHE (AcquamanHDFExplorer.py)
GUI-enabled tool for fast analysis and visualisation of the fluorescent spectra previously exported to the HDF5 format. Allows you to inspect the full MCA data in the individual channes, create/modify the ROIs while interactively monitoring the quality of the resulting EXAFS.
Uses the live EXAFS extractor created by Vadim Murzin (PETRA III P64) and modified by Konstantin Klementiev (MAX-IV Balder). See https://github.com/vadmu/EXAFS_Monitor
This tool is in development, use it at your own risk, but don't hesitate to report bugs and request features.

Use pip to install the dependencies:
pip install numpy scipy h5py matplotlib pyqtgraph

Usage:

Run the script, the use the directory tree on the left to navigate to your HDF5 files. Double click the file to open it in the viewer (patience, these files are GB-sized). Once the structure is recognized and the data is loaded, the pixel selection table will appear. Slider on top allows to view the MCA maps for individual pixels (channel 0 corresponds to the sum along all enabled channels). 
All ROIs created in Acquaman for the scan will be present in the ROI table. If no ROI was defined for the scan some arbitrary settings will be used. Red horizontal line will denote the emission peak position, yellow vertical line will be drawn at the absorption edge.
Use the matplotlib toolbar to zoom in and out. Default mouse click action (if no tool is enabled in the toolbar) is ROI selection. For better precision use the ROI panel to define the bounds.

EXAFS monitor (opens in separate window) will interactively show the counts (I0, I1, Fluo_in_ROI), mu(E) and even chi(k). Sum over all enabled pixels is used for the fluo curve, therefore you see right away what this or that pixel does with your EXAFS.

Once all bad pixels are disabled and the best ROI is identified, it's time to export the result into ASCII. To do that just set the file name and press the button in the Export panel (export status available in the terminal, let me know if you need a visible confirmation).


