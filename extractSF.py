#!/usr/bin/env /usr/bin/python3

import sys, os
import pandas


#inFile_base = sys.argv[1]


all_files = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]
fileList = [f for f in all_files if 'log' in f]
fileList = [item for item in fileList if "summary" not in item]
fileList = [item for item in fileList if ".py" not in item]

tempData = []

header = ['BE', 'J', 'pi', 'T']
isHeaderDone = False
data = []
BigExList = [] # that contain many nuclei
IsoList = []
nOrbital = 0

orbtial = ""
SFData = []

nucleusLeft = ""
nucleusRight = ""

for fileIndex, f in enumerate(fileList):
  print("!" + f)
  try:
    with open(f, "r") as file:

      content = file.read()
      
      file.seek(0)

      #========= if SF file
      if 'Spectroscopic factor' in content:
        print("Processing SF data")
        for line in file:
          columns = line.strip().split()
          if 'fn_load_wave_l' in line: nucleusLeft = (columns[2].split('_'))[0]
          if 'fn_load_wave_r' in line: nucleusRight = (columns[2].split('_'))[0]
          if len(columns) >= 7 and columns[0] == 'orbit':
            orbtial = columns[6] + "-" + columns[7]
          if len(columns) == 8 and columns[0][1] == '(':
            tempData.append(orbtial)
            tempData.append(float(columns[2]))
            tempData.append(float(columns[5]))
            tempData.append(float(columns[7]))
            tempData.append(nucleusLeft)
            tempData.append(nucleusRight)
            SFData.append(tempData)
            tempData = []


      else:
        print("Processing Ex data")

        data = []
        
        IsoList.append( f.split('_')[1])

        tempData = []
        for line in file:
          columns = line.strip().split()

          #==== header
          if( isHeaderDone == False ):
            if len(columns) >=3 and columns[0] == 'p' and columns[1] == 'orbit':
              for index, a in enumerate(columns):
                if index > 1 : header.append('p-'+a)
            if len(columns) >=3 and columns[0] == 'n' and columns[1] == 'orbit':
              for index, a in enumerate(columns):
                if index > 1 : header.append('n-'+a)
              isHeaderDone = True
              nOrbital = (len(header)-4)/2
              
          if isHeaderDone :
            #==== get the energy and orbtial occupancy
            if len(columns) == 9 and columns[1] == '<H>:' :
              tempData.append(float(columns[2]))
              tempData.append(columns[6])
              tempData.append(columns[8])
            if len(columns) == 6 and columns[0] == '<Hcm>:' :
              tempData.append(columns[5])
            if len(columns) == 4 and columns[0] == '<TT>:' :
              tempData.append(columns[3])
            if len(columns) >= nOrbital and columns[0] == '<p' and columns[1] == 'Nj>' :
              for index, a in enumerate(columns): 
                if index > 1 :  tempData.append(float(a))
            if len(columns) >= nOrbital and columns[0] == '<n' and columns[1] == 'Nj>' :
              for index, a in enumerate(columns): 
                if index > 1 :  tempData.append(float(a))
            if len(columns) >=1 and columns[0] == '<Lp>':
              data.append(tempData) 
              tempData = []

        BigExList.append(data)

  except FileNotFoundError:
    print(f"File not found: {f}")

  except IOError:
    print("Error occurred while reading the file.")


pandas.set_option('display.max_columns', None)
pandas.set_option('display.width', None)

print("================================")

SFData = sorted(SFData, key=lambda x: (x[1], x[2]))

dfSF = pandas.DataFrame(SFData, columns=['orb', 'E1', 'E2', 'SF', 'Nu1', 'Nu2'])

# for row in SFData: print(row)

# print(dfSF)

# print(isHeaderDone)
# print(header)

# print( BigExList )

print("================================ no. orbital : %d" % nOrbital)

CondensedIso = list(set(IsoList))
print(CondensedIso)
#=== sort CondensedIso by the number part
import re
def custom_sort_key(s):
    letter_part, number_part = re.match(r'([A-Za-z]+)(\d+)', s).groups()
    return (letter_part, int(number_part))

CondensedIso.sort(key = custom_sort_key)

print(CondensedIso)

#create an empty list of list
CondensedData = [ [] for _ in range(len(CondensedIso))]

for k, iso in enumerate(CondensedIso):
  for index, iso2 in enumerate(IsoList):
    if iso == iso2 :
      CondensedData[k] += BigExList[index]
  
# for k, haha in enumerate(CondensedData):
#   print("=========" + CondensedIso[k])
#   for row in haha: print(row)

DF = pandas.DataFrame()

for k, haha in enumerate(CondensedData):
  data = sorted(haha, key=lambda x: x[0])

  # print(header)
  # for row in data: print(row)
  df =  pandas.DataFrame(data, columns=header)

  #===== calculate excitation energy
  ex = [0]
  isoName = []

  if k == 0 : isoName.append(CondensedIso[k])

  for index, row in enumerate(data):
    if index > 0 :
      isoName.append(CondensedIso[k])
      ex.append(row[0] - data[0][0])

  df.insert(1, 'Ex', ex)
  df.insert(0, 'Iso', isoName)

  DF = pandas.concat([DF, df])


for iso in CondensedIso:
  print("================ " +  iso)
  print(DF[DF['Iso']==iso])

#======================== re-arrange SF data 

#Fine out which iso is the middle, assume only 3 iso

if len(CondensedIso) != 3 : exit()


orbtialList = header[4:]
# print(orbtialList)


df = DF[DF['Iso'] == CondensedIso[1]].sort_values(by=['BE'])
dft = DF[DF['Iso'] == CondensedIso[0]].sort_values(by=['BE'])

removal = dfSF[ dfSF['E1'] == round(df['BE'].iloc[0],3)]
# print(removal)
print("================================ %s remove to %s" % (CondensedIso[1], CondensedIso[0]))

removeSF = []
tempData = []
for index, row in removal.iterrows():
  tempData = []
  BE2 = row['E2']
  tempData.append(BE2)
  haha = dft[round(dft['BE'],3) == BE2]
  tempData.append(haha['Ex'].iloc[0])
  tempData.append(haha['J'].iloc[0])
  tempData.append(haha['pi'].iloc[0])
  tempData.append(haha['T'].iloc[0])
  for orb in orbtialList:
    if row['orb'] == orb :
      tempData.append(row['SF'])
    else:
      tempData.append(0)
  removeSF.append(tempData)  

removeDF = pandas.DataFrame(removeSF, columns= ['BE', 'Ex', 'J', 'pi', 'T'] + orbtialList)

print(removeDF)

if len(CondensedIso) < 3 : exit()
#===== adding reaction

df = DF[DF['Iso'] == CondensedIso[1]].sort_values(by=['BE'])
dft = DF[DF['Iso'] == CondensedIso[2]].sort_values(by=['BE'])


adding = dfSF[dfSF['E2'] == round(df['BE'].iloc[0],3)]
# print(adding)
print("================================ %s adding to %s" % (CondensedIso[1], CondensedIso[2]))

addingSF = []
tempData = []
for index, row in adding.iterrows():
  tempData = []
  BE1 = row['E1']
  tempData.append(BE1)
  haha = dft[round(dft['BE'],3) == BE1]
  tempData.append(haha['Ex'].iloc[0])
  tempData.append(haha['J'].iloc[0])
  tempData.append(haha['pi'].iloc[0])
  tempData.append(haha['T'].iloc[0])
  for orb in orbtialList:
    if row['orb'] == orb :
      tempData.append(row['SF'])
    else:
      tempData.append(0)
  addingSF.append(tempData)  

addingDF = pandas.DataFrame(addingSF, columns= ['BE', 'Ex', 'J', 'pi', 'T'] + orbtialList)

print(addingDF)

#============================ Calcualte the ESPEs base on the transfer reaction

#TODO use real data to see the code is correct

print("====== Ground state eergy:")
groundStateBE = []
for iso in CondensedIso:
  BE = DF[DF['Iso']==iso].sort_values(by='Ex')['BE'].iloc[0]
  groundStateBE.append(BE)
  print("=== %s : %f" %  (iso, BE) )

#experimental C12 adding and removal data,  the algorithm testing
# groundStateBE = [-73.4410, -92.1617, -97.1080]
# # Sample data
# data = {
#     'Ex': [0.000, 3.089, 3.685, 3.854],
#     'J': ['1/2', '1/2', '3/2', '5/2'],
#     'n-0p_1/2': [0.77, 0, 0, 0],
#     'n-0p_3/2': [0.00, 0, 0.14, 0],
#     'n-1s_1/2': [0.00, 0.65, 0.0, 0],
#     'n-0d_5/2': [0.00, 0.00, 0.0, 0.58],
#     'n-0d_3/2': [0.00, 0.00, 0.0, 0.00],
# }

# addingDF = pandas.DataFrame(data)
# # Sample data
# data = {
#     'Ex':      [0.000,   2.0,   4.3,   4.8,  6.34,  6.49,  6.92,  7.53,  8.13,  8.14],
#     'J':       ['3/2', '1/2', '5/2', '3/2', '1/2', '7/2', '5/2', '3/2', '3/2', '5/2'],
#     'n-0p_1/2': [    0,  0.61,     0,     0,     0,     0,     0,     0,     0,     0],
#     'n-0p_3/2': [  2.5,     0,     0,  0.33,     0,     0,     0,     0,  0.0175,   0],
#     'n-1s_1/2': [    0,     0,     0,     0, 0.00075,   0,     0,     0,     0,   0],
#     'n-0d_5/2': [    0,     0,     0,     0,     0,     0,0.0175,     0,     0,   0],
#     'n-0d_3/2': [    0,     0,     0,     0,     0,     0,     0,  0.01,     0,   0],
# }
# removeDF = pandas.DataFrame(data)
# print(addingDF)
# print(removeDF)

for pp, orb in enumerate(orbtialList):
  if pp < 5 : continue
  Ep = addingDF[addingDF[orb] > 0 ]
  Em = removeDF[removeDF[orb] > 0 ]
  if Ep.empty or Em.empty : continue
  print("========= " + orb)
  # print(Ep)
  # print(Em)

  EpAvg = 0
  SumSp = 0
  Gp = 0
  for index, row in Ep.iterrows():
    SumSp += row[orb]
    Gp += row[orb] * (float(row['J'][0])+1)
    EpAvg += row['Ex'] * row[orb]

  EpAvg /= SumSp

  EmAvg = 0
  SumSm = 0
  Gm = 0
  for index, row in Em.iterrows():
    SumSm += row[orb]
    Gm += row[orb]
    EmAvg += row['Ex'] * row[orb]

  EmAvg /= SumSm

  Epp = EpAvg + groundStateBE[2] - groundStateBE[1]
  Emm = -EmAvg - groundStateBE[0] + groundStateBE[1]

  print( " Ep' :  %.3f , Ep'' : %.3f, Gp : %.3f" % (EpAvg, Epp, Gp))
  print( " Em' :  %.3f , Em'' : %.3f, Gm : %.3f" % (EmAvg, Emm, Gm))

  maxOccpancy = int(orb.split('_')[1][0]) + 1
  print( " sum G :  %.3f, 2J+1 : %.0f, quenching : %.3f" % (Gp+Gm, maxOccpancy, (Gp+Gm)/maxOccpancy))

  print(" ESPS %s : %.3f" % (orb, (Epp * Gp + Emm * Gm)/(Gp + Gm)))





