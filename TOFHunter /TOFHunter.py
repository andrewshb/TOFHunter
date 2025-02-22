"""
TOFHunter (Version 1.0.5)

Description of the program can be found at the article below. 
If you use this code in your research, please cite the reference:
Hunter B. Andrews, Lyndsey Hendriks, Sawyer B. Irvine, Daniel R. Dunlap, Benjamin T. Manard. 
"TOFHunter – Unlocking Rapid Untargeted Screening by ICP-TOF-MS," Journal of Analytical Atomic Spectrometry, 2025.                                                        

For proper operation this .py file should be located within the "TOFHunter (vX.X.X)" folder with the "Reference Sheets" folder as well
"""

# all of the below libraries need to be installed for proper operation
import pandas as pd
from scipy.signal import find_peaks
import numpy as np
import h5py as h5py
import matplotlib.pyplot as plt
from alive_progress import alive_bar
import re
from sklearn.decomposition import PCA
from pathlib import Path

# copy and paste the .h5 file pathway, for example: r'/Users/CurrentUser/Desktop/TOFHunter build/my hdf5 file.h5'
filename = r''
# if you would like the peak summary table to be exported as an excel sheet please provide a file path, for example: 'r'/Users/CurrentUser/Desktop/SummaryTableName.xlsx'
# if you do not want the peak summary table to be exported please enter 'None'
export_path='None'
#set as 0.99 for 99% of the variance to be explained by PCA or enter a integer number of components
pca_components = 0.99 

#select number of unique spectra to keep from the IFF or define a percentage of the top IFF frequency to keep spectra if they are found above this frequency
unique_spectra = 0.1 

#define the criteria for peak finding algorithm, note distance is number of data points between peaks and not a direct m/z value
peak_height = 50
prominence = 5
distance = 60

def run_TOFHunter(filename, export_path, pca_components, unique_spectra, peak_height, prominence, distance):
    print('TOFHunter working...')
    
    h5_file = h5py.File(filename)
    if export_path == 'None': export = False

    #############################################################################################################
    # Define functions
    #############################################################################################################
    def find_closest_value_index(array, target_value):
        """
        Returns the index of the closest value to the target value.

        Parameters:
        - array: numpy array of values.
        - target_value: desired value to find closest.

        Returns:
        - closest index: index of the value closest to the target_value within the array.
        """
        
        abs_diff = np.abs(array - target_value)
        closest_index = np.argmin(abs_diff)

        return closest_index

    def round_half_number(number):
        return round(number*2)/2

    def check_double_charge(s, mass):
        """
        Corrects interference/mass list for doubling charged species to reflect their half mass if applicable.

        Parameters:
        - s: A string containing a nuclide name.
        - mass: mass number associated with the nuclide string.

        Returns:
        - s: a modified mass number if the nuclide is a doubly charged species.
        """
        
        if '++' in s:
            number = re.search(r'\d+', s).group()
            divided_number = int(number) / 2
            if divided_number != mass:
                return None
            else:
                return s
        else:
            return s
        
    def IFF(spectra, num_randvectors=10000):
        """
        Selection of the purest spectra in a data set.

        Parameters:
        - spectra: The data set you want to explore with spectra along the rows.
        - num_randvect: Number of random vectors to be generated (should be >= 2000).

        Returns:
        - ind: Indexes of selected spectra (sorted by decreasing order of selection frequency).
        - freq: Selection frequency of each spectrum in ind.
        """

        rows, columns = spectra.shape

        # Mean centering of the data set
        mean = np.mean(spectra, axis=0)
        spectra_centered = spectra - np.tile(mean, (rows, 1))

        # Generation of 'num_randvect' random vectors
        randvector = 2 * np.random.rand(num_randvectors, columns) - 1

        # Initiarowszation of the selection frequency rowsst
        votes = np.zeros(rows)

        # Projection of all the data onto a random vector in the loop
        with alive_bar(num_randvectors) as bar:
            for k in range(num_randvectors):
                temp = np.dot(randvector[k, :], spectra_centered.T)
                index = np.argmax(temp)
                votes[index] += 1
                index = np.argmin(temp)
                votes[index] += 1
                
                bar()
        # Sorting spectra in a decreasing order of frequency selection
        sorted_indices = np.argsort(votes)[::-1]
        freq = votes[sorted_indices]
        indices = sorted_indices

        return indices, freq
    
    #############################################################################################################
    # Extract the file attributes and peak data
    #############################################################################################################

    write = h5_file['/'].attrs['NbrWrites'][0]
    buf = h5_file['/'].attrs['NbrBufs'][0]
    segment = h5_file['/'].attrs['NbrSegments'][0]

    massAxis = np.array(h5_file.get('/FullSpectra/MassAxis'))
    sumSpectrum = np.array(h5_file.get('/FullSpectra/SumSpectrum'))


    temp = h5_file.get('/PeakData/PeakData')
    peakData = np.zeros(np.shape(temp))
    temp.read_direct(peakData)
    peakData = peakData.reshape(-1, peakData.shape[-1])

    peakAxis = np.array(h5_file.get('/PeakData/PeakTable')[:,'mass'])
    peakLabels = np.array(h5_file.get('/PeakData/PeakTable'))
    print('data loaded')
    
    #############################################################################################################
    # Run PCA analysis on peak data
    #############################################################################################################

    #set as 0.99 for 99% of the variance to be explained by PCA 
    pca = PCA(n_components=pca_components)
    #or enter a integer number of components
    # pca = PCA(n_components=10)

    TOFpca = pca.fit_transform(peakData)

    loadings = pca.components_
    explained_variance_ratio = 100*pca.explained_variance_ratio_
    cumulative_variance = np.cumsum(explained_variance_ratio)

    loadingAxis = np.arange(peakAxis[0],peakAxis[-1],0.2)
    int_loadings = np.zeros((len(loadings), len(loadingAxis)))
    for i in range(len(peakAxis)):
        index = find_closest_value_index(loadingAxis, peakAxis[i])
        int_loadings[:,index] = loadings[:,i]


    plt.close(1)
    plt.figure(1,(8,4))
    ax1 = plt.subplot2grid((1,4),(0,0))
    ax2 = plt.subplot2grid((1,4),(0,1),colspan=3)
    ax1.plot(np.arange(1, len(explained_variance_ratio) + 1), explained_variance_ratio, '.-')
    ax1.set_xlabel('Principal Component')
    ax1.set_ylabel('Explained Variance Ratio')

    ax2.plot(loadingAxis,int_loadings.T)
    ax2.set_xlabel('Mass-to-charge Ratio (Th)')
    ax2.set_ylabel('Loading Value')
    plt.legend(np.arange(1, len(explained_variance_ratio) + 1))
    ax2.plot(peakAxis[[0,-1]],[0,0],'k',linewidth=0.3)
    # plt.yscale('log')
    plt.title('Principal Component Analysis')
    plt.tight_layout()
    print('pca complete')

    #############################################################################################################
    #Run IFF analysis on peak data and extract full spectra
    #############################################################################################################

    # unique_spectra= 10 #select number of unique spectra to keep from the IFF
    # or define a percentage of the top IFF frequency to keep spectra if they are found above this frequency

    ind, freq = IFF(peakData, 2000)
    if isinstance(unique_spectra, int) == False:
        threshold = unique_spectra*freq[0]
        unique_spectra = int(len(freq[freq>threshold]))
    top_spectra = list(ind[:unique_spectra])

    Data = np.empty([len(top_spectra), len(massAxis)])
    i = 0

    for i in range(0, len(top_spectra)):
        spectrum = top_spectra[i]
        w = int(np.floor(spectrum/(buf*segment)))
        s = np.remainder(spectrum, (buf*segment))/segment
        if s<1:
            s = int(s*segment)
            b = 0
        else:
            b = int(np.floor(s))
            s = int(round((s-b)*segment))

        Data[i, :] = np.array(h5_file.get('/FullSpectra/TofData')[w, b, s, :])

               
    plt.close(2)
    plt.figure(2,(8,4))
    ax1 = plt.subplot2grid((1,4),(0,0))
    ax2 = plt.subplot2grid((1,4),(0,1),colspan=3)
    ax1.plot(np.arange(1, len(top_spectra)+1), freq[:len(top_spectra)], '.-')
    ax1.set_xlabel('Unique Spectra')
    ax1.set_ylabel('IFF Frequency')

    ax2.plot(massAxis, Data.T)
    ax2.set_xlabel('Mass-to-charge Ratio (Th)')
    ax2.set_ylabel('Intensity (a.u.)')
    plt.legend(np.arange(1, len(top_spectra)+1))
    ax2.plot(peakAxis[[0,-1]],[0,0],'k',linewidth=0.3)
    # plt.yscale('log')
    plt.title('Interesting Feature Finders')
    plt.tight_layout()     
    print('IFF complete')
    
    #############################################################################################################
    # Run peak finder
    #############################################################################################################

    spectra_list = pd.DataFrame()
    peak_list = pd.DataFrame()
    counts_list = pd.DataFrame()

    for i in range(0, unique_spectra):
        peaks_idx, _ = find_peaks(
            Data[i, :].T, height=peak_height, prominence = prominence, distance = distance)
        peaks = pd.DataFrame(massAxis[peaks_idx])
        counts = pd.DataFrame(Data[i, peaks_idx])
        spectra_list = pd.concat([spectra_list, pd.DataFrame(
            np.ones([len(peaks)])*i)], axis=0, ignore_index=True)
        peak_list = pd.concat([peak_list, peaks], axis=0, ignore_index=True)
        counts_list = pd.concat([counts_list, counts], axis=0, ignore_index=True)
    peak_summary = pd.concat([spectra_list+1, peak_list, counts_list], axis=1)

    #example peak finder
    spectrum_number=0
    x = Data[spectrum_number-1, :].T
    peaks, _ = find_peaks(x, height=peak_height, prominence=prominence, distance=distance)
    peak_list = massAxis[peaks]

    plt.close(3)
    plt.figure(3,(8,4))
    plt.plot(massAxis, x)
    plt.plot(massAxis[peaks], x[peaks], "x")
    plt.xlabel('Mass-to-charge Ratio (Th)')
    plt.ylabel('Intensity (a.u.)')
    plt.legend(['Mass Spectrum','Identified Peaks'])
    plt.plot(peakAxis[[0,-1]],[0,0],'k',linewidth=0.3)
    plt.title('Peak Finder Example')
    plt.tight_layout()
    plt.show()

    #############################################################################################################
    # Run peak identification and match with nuclides and interferences
    #############################################################################################################

    NISTAMU = pd.read_excel(
        str(Path.cwd())+r'/Reference Sheets/NISTAMU.xlsx')
    interference = pd.read_excel(
        str(Path.cwd())+r'/Reference Sheets/ICP elements and interferents_revised.xlsx',
        sheet_name='Ion DB', header = 1)

    interference = interference.loc[interference['Ion type'] != 'elemental']

    # This loop searches through the NISTAMU list to concat matched names within an acceptence to each peak
    rounded_peak_summary = round_half_number(peak_summary.iloc[:, 1])
    dfNIST = pd.DataFrame(np.full((1000, len(peak_summary)), 'nan'))
    dfINF = pd.DataFrame(np.full((1000, len(peak_summary)), 'nan'))

    for i in range(len(rounded_peak_summary)):
        temp_peak_match = NISTAMU[NISTAMU.iloc[:,2] == rounded_peak_summary[i]].iloc[:,5]
        temp_INF = interference[interference.iloc[:,0] == rounded_peak_summary[i]].iloc[:,3] #checks if interference is in peak match column
        temp_INF =temp_INF.drop_duplicates()
        dfNIST.iloc[:len(temp_peak_match),i] = temp_peak_match
        dfINF.iloc[:len(temp_INF),i]=temp_INF
        
    dfNIST = dfNIST.T
    dfINF = dfINF.T

    dfNIST['compressed'] = dfNIST.apply(lambda row: row.tolist(),axis=1)
    dfNIST = dfNIST.iloc[:,[-1]]

    dfINF['compressed'] = dfINF.apply(lambda row: row.tolist(),axis=1)
    dfINF = dfINF.iloc[:,[-1]]

    for i in range(len(dfNIST)):
        dfNIST.iloc[i,0] =  [x for x in dfNIST.iloc[i,0] if x != 'nan']
        dfINF.iloc[i,0] = [x for x in dfINF.iloc[i,0] if x != 'nan']

    summary = pd.concat([peak_summary,dfNIST,dfINF],axis=1)
    summary.columns = ['Spectrum', 'm/z','Intensity (a.u.)','Nuclide Matches','Potential Interferences']

    if export == True:
        summary.to_excel(export_path,header=True, index=False)
    print('export complete')

run_TOFHunter(filename, export_path, pca_components, unique_spectra, peak_height, prominence, distance)