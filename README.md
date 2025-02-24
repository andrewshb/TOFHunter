# TOFHunter
![TOC](https://github.com/user-attachments/assets/cbda75ad-4107-4a24-b352-0e0b365a5146)

TOFHunter was designed to screen ICP-TOF-MS data from TOFWERK icpTOF instruments. The programs are written in python and are provided as .py files.

## Installation: 
Users can download and run TOFHunter either by:
Directly downloading a copy of this repository 
or by performing a git clone (https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)

### Windows
1. Open Git Bash.
2. Change the current working directory to the location where you want the cloned directory.
3. Type git clone, and then paste the URL you copied earlier.
git clone https://github.com/andrewshb/TOFHunter

### MAC
1. Open Terminal.
2. Change the current working directory to the location where you want the cloned directory.
3. Type git clone, and then paste the URL you copied earlier.
git clone https://github.com/andrewshb/TOFHunter

### Folder configuration
Within the TOFHunter folder you will find two .py files 
1. TOFHunter.py
2. TOFHunter GUI.py

These can be run in the users prefered python environment; however, they were designed for use with the Spyder IDE which can be installed using Anaconda (anaconda.org).
Both files should be kept in the same working directory as included the "Reference Sheets" and "assets" folders to function properly. 
Data files can be saved within this directory or elsewhere.

### TOFHunter.py
The script can be run by clicking "Run File" or (F5)![image](https://github.com/user-attachments/assets/38ea7451-ad14-4429-aeaf-e54243ce2c80)
 in Spyder. 
**An .h5 file path must be entered on line 24. **
Other features can be toggled on lines 23 - 37. Note, modifying lines 38 onward may effect TOFHunter operation.
![image](https://github.com/user-attachments/assets/6c7b6b7a-41e0-412b-b3ec-098ab9303ec3)


### TOFHunter GUI.py
The script can be run by clicking "Run File" or (F5)![image](https://github.com/user-attachments/assets/5378682e-5393-48fa-a217-f70c19b3c6d0)
 in Spyder. 
**The Spyder graphics setting must be set to tkinter **
To modify this click the "preferences" icon ![image](https://github.com/user-attachments/assets/be9c0168-6d2b-46db-bb96-2244f53a5c7e)
Within the preferences window, select the "IPython console" option in the left pane. 
Next select the "Graphics" tab and change the Backend selection to "Tkinter". 
Lastly, select apply to modify settings. This should only need to be completed once unless Spyder is reinstalled. 
See the image below with annotated directions.
![image](https://github.com/user-attachments/assets/001b9bff-faea-447a-9f53-b7827930fc4f)


A GUI window will appear (see below) and files and features can be modified within the GUI. 
Note, GUI buttons may need to be clicked several times in the current version. The button will appear depressed when the click is registered. 
Default feature values can be toggled on lines 30-45. Note, modifying lines 48 onward may effect TOFHunter GUI operation.
In some cases, closing the GUI window will not terminate the script from running. If this happens you should click the 'interupt kernel' ![image](https://github.com/user-attachments/assets/4fcbfaa5-8460-4216-ade6-d0db519c9671) button in Spyder's console.

![image](https://github.com/user-attachments/assets/b6f954fb-d621-410a-871b-f4da2db7db0d)


_Additional details are provided within each script._

**Happy Hunting!**

## Reference:
For more information on TOFHunter and example applications of its use see the associated publication. We ask that if TOFHunter is used in your study that the below article be cited. 

Hunter B. Andrews, Lyndsey Hendriks, Sawyer B. Irvine, Daniel R. Dunlap, Benjamin T. Manard. 
"TOFHunter – Unlocking Rapid Untargeted Screening by ICP-TOF-MS," _Journal of Analytical Atomic Spectrometry_, 2025.   

## Required Python Libraries:
If using Spyder IDE, it is recommended to use "conda install package" for all the below. Note, alive_progress can only be installed via pip. Mixing conda and pip may cause issues with Spyder. If you experience issues with Spyder, uninstall and reinstall. Then consider a new environment be made to use this program.
For more information on creating a secondary environment for Spyder see the linked video at timestamp 6:18 (https://www.youtube.com/watch?v=Ul79ihg41Rs). 

pathlib (https://docs.python.org/3/library/pathlib.html)

tkinter (https://docs.python.org/3/library/tkinter.html)

matplotlib (https://matplotlib.org/)

pandas (https://pandas.pydata.org/)

scipy (https://scipy.org/)

numpy (https://numpy.org/)

h5py (https://docs.h5py.org/en/stable/)

alive_progress (https://github.com/rsalmei/alive-progress)

re (https://docs.python.org/3/library/re.html)

scikit-learn (https://scikit-learn.org/stable/)

openpyxl (https://openpyxl.readthedocs.io/en/stable/)
