""""
In this code, intervals between detected beats are here referred to 
as RRIs (R peak – to R peak – intervals), as is applicable
for electrocardiogram (ECG) data. The actual correct term when using 
photopletysmography (PPG) data is IBIs (Interbeat Intervals).”
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# The module preprocessing can be found at
# https://github.com/Aura-healthcare/hrv-analysis/blob/master/hrvanalysis/preprocessing.py
from preprocessing import interpolate_nan_values


"""
Note: 
This file takes as input the file LL_CumData_RRI.xlsx. 
This file has both the cumulative data, as well as the corresponding 
RR interval (RRI) data.
"""

df = pd.read_excel(r'LL_CumData_RRI.xlsx')

"""
The below comment turns the dataframe df, containing the inputted Excel data, 
into a matrix where 
- The first column has the measurements numbers
- The second column the cumulative time data
- The third column has the RRI data
"""

data = df.to_numpy()

"""
The function 'mad_detection' defines the quotient filter which indicates all RRI's that 
differ more than a multiplicative factor alpha with one of its neigbouring values.

It replaces such NNI's with a NaN-value
"""


def mad_detection(rri,alpha=0.2):
    n = len(rri) # Number of rri's per measurement
    rri_nan = rri.copy().astype('float')
    for i in range(0,n): 
        if (i > 0 and i < n-1): #Let's assume the first and last rri-value are correct
            if ((rri[i]/rri[i-1] < 1 - alpha or rri[i]/rri[i-1] > 1 + alpha) and \
                (rri[i]/rri[i+1] < 1-alpha or rri[i]/rri[i+1] > 1+alpha)):
                    rri_nan[i] = np.nan   
    return rri_nan

"""
The function below replaces all NaN-values by an appropriate amount of 
substitute beats (if beats were missed, NaN-values are added, 
                  if there are too many beats, some are deleted)

Function is designed to detect isolated errors.
It may not detect all series of faulty rri-measurements.
"""

def mad_beats_insert_delete(rri,alpha=0.2):
    #Input is rri values 
    n = len(rri)
    
    #rri_nan gives rri with outlier replaced by np.nan
    rri_nan = mad_detection(rri,alpha)

    #The below parameters are statistics concerning deletion/addition
    intervals_added = 0 #Counts in how many intervals beats added 
    intervals_deleted = 0 # " " " deleted
    intervals_same = 0 #Measures 'ectopic beats', i.e., intervals where number of beats did not change
    total_beats_added = 0 #Count total numbe of beats added over all relevant intervals
    total_beats_deleted = 0 # " " " deleted " " "
    i = 1
    while i < n - 1:     
        mask = np.isnan(rri_nan)
        if mask[i] == True:
            #Start of NaN-interval
            i_start = i

            #Determine end of interval with NaN-values
            if mask[i+1] == True:
                i_end = np.argmin(mask[i+1:]) + i
            else:
                i_end = i
            
            #Total length of interval
            total_nan = np.sum(rri[i_start:i_end+1])
            
            #Determine #beats to insert
            average = np.mean([rri[i_start-1],rri[i_end+1]])
            beats_to_insert = np.round(total_nan/average).astype('int')
            
            #Update relevant statistics parameters
            diff = i_end+1 - i_start
            if beats_to_insert == diff:
                intervals_same += 1
            elif beats_to_insert > diff:
                intervals_added += 1
                total_beats_added += beats_to_insert - diff
            elif beats_to_insert < diff:
                intervals_deleted += 1
                total_beats_deleted += diff - beats_to_insert
            
            #Here we delete the NaN-values in the current interval
            rri_nan = np.delete(rri_nan,np.s_[i_start:i_end+1])
            rri = np.delete(rri,np.s_[i_start:i_end+1])
            
            #Here we add corrected number of NaN-values in current interval
            count = beats_to_insert
            while(count > 0):
                rri_nan = np.insert(rri_nan,i_start,np.nan)
                rri = np.insert(rri,i_start,0)
                count = count - 1
                
            i =  i_start + beats_to_insert + 1
            n = len(rri_nan) #Redefine length of rri if it has been adjusted
        else:
            i = i+1
            
        #Store statistics in dictionary    
        correction_stats = {
            'Intervals with beats deleted' : intervals_deleted,
            'Total number of beats deleted' : total_beats_deleted,
            'Intervals with beats added' : intervals_added,
            'Total number of beats added' : total_beats_added,
            'Number of ectopic beats' : intervals_same
                }
    return rri_nan, correction_stats

"""
This function performs linear interpolation on NaN values
"""
    
def mad_nan_interpolation(rri,alpha=0.2):
    rri_corr = mad_beats_insert_delete(rri,alpha)[0]
    rri_interpolated = \
        interpolate_nan_values(rri_corr, \
                               interpolation_method="linear")
    return np.array(rri_interpolated)

  
"""
Below is function for plotting
"""

def rri_comparison_plot(meas_num,data,alpha):
    """
    Takes as input measurement number and big data set
    """
    mask = data[:,0] == meas_num #Select data corresponding to 'meas_num'

    rri = data[mask,2] #Gives the NNI data of measurement j
    rri = np.delete(rri, 0, axis=0) #Delete first zero in rri sequence
    rri_cum = np.cumsum(rri)

    rri_corr = mad_nan_interpolation(rri,alpha)
    rri_corr_cum = np.cumsum(rri_corr)

    plt.figure()
    plt.plot(rri_cum,rri,label='Original')
    plt.plot(rri_corr_cum,rri_corr,label='Corrected')
    plt.legend()
    plt.title('Measurement %i' % meas_num)
    plt.savefig('Comparison for Measurement_%i.png' % j)
    
    time_diff = rri_cum[-1] - rri_corr_cum[-1]
    return time_diff

"""
Below we correct the data
"""

#Choice of alpha determines strictness of the quotient filter 
alpha = 0.2

total_meas = 1280 # Total number of measurements in input data
stat_count = 0
stats_info = np.zeros((total_meas,5+1)) #The '5' comes from that we compute five statistics in mad_beats_insert_delete

for j in range(1,total_meas+1):
    mask = data[:,0] == j #Select data corresponding to meas. 'j'
    rri = data[mask,2] #Gives the raw NNI/RRI data of measurement j
    rri = np.delete(rri, 0, axis=0) #Delete zero at the beginning
    
    #Compute corrected (cumulative) rri data
    rri_corr = mad_nan_interpolation(rri,alpha)
    rri_corr_cum = np.cumsum(rri_corr)   
    n = len(rri_corr)
    
    #Column with measurement number and new data
    meas_col = j*np.ones((n,1)).astype('int')    
    new_meas_data = np.hstack((meas_col,rri_corr_cum[:,None],rri_corr[:,None]))
    
    if j == 1:
        new_data = new_meas_data
    elif j > 1:
        new_data = np.vstack((new_data,new_meas_data))

    meas_corr_stats = mad_beats_insert_delete(rri,alpha=0.2)[1]
    meas_stats = np.array(list(meas_corr_stats.values()))
    stats_info[j-1,:] = np.hstack((np.array([j]),meas_stats))
    
    rri_comparison_plot(j,data,alpha)
    
"""
Writing corrected data to Excel file
"""

column_names = ['Meas. number', \
                'Corrected cum. data ' , \
                'Corrected RRI data']

df_new= pd.DataFrame(new_data,columns=column_names)


stats_col_names = ['Meas. number', \
                'Number of intervals with beats del.', \
                'Total number of beats deleted', \
                'Number of intervals with beats add.', \
                'Total number of beats added.', \
                'Number of ectopic beats']
 
df_new_stats = pd.DataFrame(stats_info,columns=stats_col_names)

# Create an Excel writer object
with pd.ExcelWriter('LL_CumData_RRI_corrected.xlsx') as writer:
    # Write each dataframe to a different sheet
    df_new.to_excel(writer, sheet_name='CumData_RRI, alpha=%.2f' % alpha, index=False)
    df_new_stats.to_excel(writer, sheet_name='Add-del counts, alpha=%.2f' % alpha, index=False)
    
