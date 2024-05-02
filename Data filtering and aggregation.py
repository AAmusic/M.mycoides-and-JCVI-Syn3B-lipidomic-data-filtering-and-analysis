
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import openpyxl
from pandas import ExcelWriter
from scipy.optimize import curve_fit

#PART 1. CHOOSE THE LIPIDOMIC SUBSET, REMOVE THE STORAGE LIPIDS, FILTER THE DATA TO 98% CUMLULATIVE AVERAGE CUTOFF AND AGGREGATE DATA BY DIFFERENT CATEGORIES (OPTIONAL)

#keywords to operate parts of the dataset. 'strain' contatins two keywords that allow to choose M. mycoides lipidomes ('MM') or Syn3B lipdiomes ('S3'). 'dataset' picks lipidomes that were acquired from cells, equlibrated to grow at different temperatures (= 'adapted'),
#and those that are obtained from cells, subjected to rapid temperature change, followed by 24 hours adaptation to a new temperature (='temporal')
strain = ['MM', 'S3']
dataset = ['adapted', 'temporal']

def pick_data(file_path, keyword1, keyword2):
    sample_data = []
    data  = pd.ExcelFile(file_path)
    for sheet_name in data.sheet_names:
        if keyword1 in sheet_name and keyword2 in sheet_name:
            df = pd.read_excel(file_path, sheet_name = sheet_name)
            sample_data.append(df)
    return sample_data
#The path to call the main data. Upoin downloading the data, insert the customary file path before the file name to call all the lipidomic data. The file name = 'TempAdapt_Manuscript_Data.xlsx'
file_path = '/file path/TempAdapt_Manuscript_Data.xlsx'
#The example code here samples lipdiomic for M. mycoides lipidomes, adapted to different growth temperatures
keyword1 = 'MM'
keyword2 = 'adapted'

#Sample_data is the raw lipdiomic dataset for a chosen bacterial strain and number of growth conditions:
sample_data = pick_data(file_path, keyword1, keyword2)


#Filter the data to 98% cumulative average cutoff threshold:
Filtered = []
for df in sample_data:
    no_storage_lipids = df.query('classes != "CE" & classes != "TAG"')
    d = no_storage_lipids[['R_1', 'R_2', 'R_3']]
    dsum = d.sum(numeric_only = True, axis = 0)
    dn = (d/dsum)*100
    dj = pd.concat([no_storage_lipids[['feature','classes','totallength','totaldb']], dn], axis = 1)
    TAV = dj[['R_1', 'R_2', 'R_3']].mean(axis = 1)
    dj['TAV'] = TAV
    dj = dj.sort_values(by = 'TAV', ascending = False)
    dj['cumTAV'] = dj['TAV'].cumsum()
    filtered = dj.query('cumTAV < 98')
    mainclasses = filtered.loc[filtered['classes'].isin(['Chol', 'PG', 'CL', 'PC', 'SM'])]
#Recalculate the remaining values to 100 mol%. Skip this step if proceed with further data aggregation by category or class
    f = mainclasses[['R_1', 'R_2', 'R_3']]
    fsum = f.sum(numeric_only = True, axis = 0)
    fn = (f/fsum)*100
    final = pd.concat([mainclasses[['feature','classes','totallength','totaldb']], fn], axis = 1)
    Filtered.append(final)

#Filtered data export to excel (at this or any other stage):
    #writer = ExcelWriter("/Volumes/Present/Mycoplasma on FBS dataset/16.xlsx")
    #final.to_excel(writer, sheet_name= 'all') 
    #writer.save()

#To aggregate the data by class (or call a particular species, parameter):
    PG = final.query('classes == "PG"')
    CL = final.query('classes == "CL"')
    PC = final.query('classes == "PC"')
    SM = final.query('classes == "SM"')
#Next, proceed with the above steps to recalculate the species sum to 100% in a class



#PART II. POWER FIT OF THE LIPIDOMES
#First, need to merge all 5 filtered lipidomes for all growth temperatures into one dataset to calculate its total average
Merged = Filtered[0]
for df in Filtered[1:]:
    Merged = pd.merge(Merged, df, on = 'feature', how = 'outer')
symbol = 'R_'
reps = [col for col in Merged.columns if symbol in col]
TAVM = Merged[reps].mean(axis = 1)
Merged['TAVM'] = TAVM
Merged = Merged.sort_values(by = 'TAVM', ascending = False)
Merged.columns = [ 'species', 'classes37', 'length37', 'db37', 'R1_37', 'R2_37', 'R3_37', 'classes33.5', 'length33.5', 'db33.5', 'R1_33.5', 'R2_33.5', 'R3_33.5', 'classes30', 'length30', 'db30', 'R1_30', 'R2_30', 'R3_30', 'classes27', 'length27', 'db27', 'R1_27', 'R2_27', 'R3_27', 'classes25', 'length25', 'db25', 'R1_25', 'R2_25', 'R3_25', 'TAVM']
classes_all = Merged['classes37'].fillna('') + Merged['classes33.5'].fillna('') + Merged['classes30'].fillna('') + Merged['classes27'].fillna('') + Merged['classes25'].fillna('')
classes_all = [s[:2] if s.startswith('PG') else s for s in classes_all]
classes_all = [s[:2] if s.startswith('CL') else s for s in classes_all]
classes_all = [s[:2] if s.startswith('PC') else s for s in classes_all]
classes_all = [s[:2] if s.startswith('SM') else s for s in classes_all]
length_all = Merged['length37'].combine_first(Merged['length33.5']).combine_first(Merged['length30']).combine_first(Merged['length27']).combine_first(Merged['length25'])
db_all = Merged['db37'].combine_first(Merged['db33.5']).combine_first(Merged['db30']).combine_first(Merged['db27']).combine_first(Merged['db25'])
Merged['classes'] = classes_all
Merged['length'] = length_all
Merged['db'] = db_all
Merged.drop(columns = [ 'classes37', 'length37', 'db37', 'classes33.5', 'length33.5', 'db33.5', 'classes30', 'length30', 'db30','classes27', 'length27', 'db27','classes25', 'length25', 'db25'], inplace = True)

#Skip the above step and proceed with power fit right away if wish to obtain the fit for each of the temperatures in the dataset
#Define the power fit function to fit sorted mycoplasma lipidomes:
x = list(range(1, len(Merged['TAVM']) +1))

def f(x,a,b):
    return a*np.power(x,b)
popt, pcov = curve_fit(f,x, Merged['TAVM'])
par = popt[0]*np.power(x, popt[1])
#in the above equation, popt[1] represents the actual value of a power fit, while popt[0] is a scale constant of the power equation




#PART III. CALCULATE LIPIDOMIC REMODELING ACROSS TEMPERATURES
#Keep using the merged dataset to be able to directly subtract the species abundances at different temperatures
LRspecies2537 = (Merged[['R1_37', 'R2_37', 'R3_37']].mean(axis = 1)- Merged[['R1_25', 'R2_25', 'R3_25']].mean(axis = 1)).abs()
SD = np.sqrt(((Merged[['R1_37', 'R2_37', 'R3_37']].var(axis = 1)) + Merged[['R1_25', 'R2_25', 'R3_25']].var(axis = 1))/2)
TotalLR = LRspecies2537.sum(axis = 0)
TotalLRSD =np.sqrt(((Merged[['R1_37', 'R2_37', 'R3_37']].var(axis = 1) + Merged[['R1_25', 'R2_25', 'R3_25']].var(axis = 1))/2).sum(axis = 0)/(len(Merged['TAVM']) +1))
#the above code can simply be used to look at other temperature differences and combined with the above parts to look at subssets of lipidome


#PART IV. CALCULATE AVERAGE ACYL CHAIN FEATURES OF PHOSPHOLIPIDS
#Use length and db (double bonds) columns, provided in the original dataset to calculate the average entry fro length and unsaturation per lipid. For PG, PC and SM, divide total length and unsaturation by 2, for CL - by 4
#Start modifying the merged dataset by excluding cholesterol from the lipidome to only look at phosphplipids. In the similar manner, phospholipids can be further grouped by glycerophospholipids and sphingomyelin
Merged.drop(columns = ['TAVM'], inplace = True)
PLs = Merged.query('species != "Chol"').fillna(0)
symbol2 = 'R'
reps2 = [col for col in PLs.columns if symbol2 in col]
psum = PLs[reps2].sum(axis = 0)
pnew = (PLs[reps2]/psum)*100
AVL = PLs.apply(lambda row: row['length'] / 4 if row['classes'] == 'CL' else row['length'] / 2, axis=1)
AVDB = PLs.apply(lambda row: row['db'] / 4 if row['classes'] == 'CL' else row['db'] / 2, axis=1)
PLsnew = pd.concat([PLs[['species', 'classes', 'length', 'db']], AVL, AVDB, pnew], axis = 1)
#Now that we have average length and unsaturation per species, we can calculate the weighted average of average acyl chain features based on species abundance in phospholipidome:
abundances = PLsnew[reps2]
AverageL = pd.DataFrame()
AverageDB = pd.DataFrame()
for col in abundances:
    AverageL[col + 'AVL'] = (PLsnew[col]*PLsnew[0]/100)
    AverageDB[col + 'AVDB'] = (PLsnew[col]*PLsnew[1]/100)
AVLfinal = AverageL.sum(axis = 0)
AVDBfinal = AverageDB.sum(axis = 0)
#AVLfinal and AVDBfinal are the average acylchain lenght and unsaturation, repsectively, for mycoplasma phospholipidome, calcualted for 3 replicates at 5 growth temperatures.






