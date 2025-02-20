# TOFHunter
![TOC](https://github.com/user-attachments/assets/cbda75ad-4107-4a24-b352-0e0b365a5146)

TOFHunter was designed to screen ICP-TOF-MS data from TOFWERK icpTOF instruments. The programs are written in python and are provided as .py files.

## For use: 
Download and unzip the folder from the most current released version. 
Within this folder you will find two .py folders - TOFHunter.py and TOFHunter GUI.py
These can be run in the users prefered python environment; however, they were designed for use with the Spyder IDE (/https://www.spyder-ide.org/).
Both files should be kept in the same directory as included the "Reference Sheets" and "assets" folders to function properly. 
Data files can be saved within this directory or elsewhere.

### TOFHunter.py
The script can be run by clicking "Run File" or (F5)![image](https://github.com/user-attachments/assets/38ea7451-ad14-4429-aeaf-e54243ce2c80)
 in Spyder. 
**An .h5 file path must be entered on line 23. **
Other features can be toggled on lines 22 - 36. Note, modifying lines 38 onward may effect TOFHunter operation.
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
Default feature values can be toggled on lines 35-45. Note, modifying lines 48 onward may effect TOFHunter GUI operation.
![image](https://github.com/user-attachments/assets/b6f954fb-d621-410a-871b-f4da2db7db0d)


_Additional details are provided within each script._

**Happy Hunting!**

## Reference:
For more information on TOFHunter and example applications of its use see the associated reference 

Hunter B. Andrews, Lyndsey Hendriks, Sawyer B. Irvine, Daniel R. Dunlap, Benjamin T. Manard. 
"TOFHunter – Unlocking Rapid Untargeted Screening by ICP-TOF-MS," _Journal of Analytical Atomic Spectrometry_, 2025.   

## Required Python Libraries:
If using Spyder IDE, it is recommended to use "conda install package" for all the below. Note, alive_progress can only be installed via pip. Mixing conda and pip may cause issues with Spyder and it is recommended that a new environment be made to use this.
For more information on creating a secondary environment for spyder see the linked video at timestamp 6:18 (https://www.youtube.com/watch?v=Ul79ihg41Rs). 

pathlib (https://docs.python.org/3/library/pathlib.html)
tkinter (https://docs.python.org/3/library/tkinter.html)
matplotlib (https://matplotlib.org/)
pandas (https://pandas.pydata.org/)
scipy (https://scipy.org/)
numpy (https://numpy.org/)
h5py (https://docs.h5py.org/en/stable/)
os (https://docs.python.org/3/library/os.html)
alive_progress (https://github.com/rsalmei/alive-progress)
re (https://docs.python.org/3/library/re.html)
scikit-learn (https://scikit-learn.org/stable/)
openpyxl (https://openpyxl.readthedocs.io/en/stable/)
