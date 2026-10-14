#!/usr/bin/env python3 
#use the pthon3 in $PATH

#================== The folder containing the log_* files must be given as the
#================== first argument (run with -h for usage).

import argparse, os, sys
import pandas

def build_parser():
  parser = argparse.ArgumentParser(
      description="Parse kshell log_* files: extract the energy lists and the "
                  "spectroscopic factors (C^2*S), build the adding/removal tables, "
                  "and compute the ESPE of each neutron orbital from the transfer data.",
      epilog="Usage notes:\n"
             "  * Any file in the folder whose name contains 'log' (and does not\n"
             "    contain 'summary', and is not a .py file) is used.\n"
             "  * Files containing 'Spectroscopic factor' are treated as transit (tr)\n"
             "    logs; all other log files are treated as single-nucleus energy logs.\n"
             "  * Energy logs are grouped by isotope, taken from the second '_' field\n"
             "    of the file name (log_N15_ysox_m1p.txt -> N15); all partitions of\n"
             "    the same isotope are combined into one spectrum.\n"
             "  * Isotopes are sorted by element name and mass number. The removal\n"
              "    reaction is isotope[1] -> isotope[0] and the adding reaction is\n"
              "    isotope[0] -> isotope[1] (ground state of the donor nucleus).\n"
              "  * A folder with a single nucleus (or without 'Spectroscopic factor'\n"
              "    files) only gets the combined spectrum printed; the transfer and\n"
              "    ESPE sections are skipped.\n"
             "  * The ESPE is referenced to the ground-state energy of the lightest\n"
             "    isotope (negative = bound). An orbital is skipped if it has no\n"
             "    strength on either the adding or the removal side.\n"
             "  * The first 5 orbitals of the model space are assumed to be proton\n"
             "    orbitals and are skipped in the ESPE section.",
      formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("folder", help="folder containing the log files")
  return parser

parser = build_parser()
if len(sys.argv) < 2:
  parser.print_help()
  sys.exit(0)
args = parser.parse_args()

all_files = [f for f in os.listdir(args.folder) if os.path.isfile(os.path.join(args.folder, f))]
fileList = [f for f in all_files if 'log' in f] # pick any file has 'log'
fileList = [item for item in fileList if "summary" not in item] # pick any file does not have 'summary'
fileList = [item for item in fileList if ".py" not in item] # pick any file does not have '.py'

tempData = []

header = ['BE', 'J', 'pi', 'T']
isHeaderDone = False
tempData = []
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
    with open(os.path.join(args.folder, f), "r") as file:

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
            tempData = []
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
        
        isoParts = f.split('_')
        IsoList.append( isoParts[1] if len(isoParts) > 1 else os.path.splitext(f)[0] )

        for line in file:
          columns = line.strip().split()

          #==== header
          if( isHeaderDone == False ):
            if len(columns) >=3 and columns[0] == 'p' and columns[1] == 'orbit':
              for index, a in enumerate(columns):
                if index > 1 : header.append('p-'+a)
              continue
            if len(columns) >=3 and columns[0] == 'n' and columns[1] == 'orbit':
              for index, a in enumerate(columns):
                if index > 1 : header.append('n-'+a)
              isHeaderDone = True
              nOrbital = (len(header)-4)/2
              continue
              
          if isHeaderDone :
            #==== get the energy and orbtial occupancy
            if len(columns) == 9 and columns[1] == '<H>:' :
              tempData.append(float(columns[2]))
              tempData.append(columns[6])
              tempData.append(columns[8])
              continue
            if len(columns) == 6 and columns[0] == '<Hcm>:' :
              tempData.append(columns[5])
              continue
            if len(columns) == 4 and columns[0] == '<TT>:' :
              tempData.append(columns[3])
              continue
            if len(columns) >= nOrbital and columns[0] == '<p' and columns[1] == 'Nj>' :
              for index, a in enumerate(columns): 
                if index > 1 :  tempData.append(float(a))
              continue
            if len(columns) >= nOrbital and columns[0] == '<n' and columns[1] == 'Nj>' :
              for index, a in enumerate(columns): 
                if index > 1 :  tempData.append(float(a))
              continue
            if len(tempData) > 0 and   ('------' in line):
              data.append(tempData) 
              tempData = []
              continue

        BigExList.append(data)
        # print(BigExList)

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
# print(CondensedIso)
#=== sort CondensedIso by the number part
import re
def custom_sort_key(s):
    letter_part, number_part = re.match(r'([A-Za-z]+)(\d+)', s).groups()
    return (letter_part, int(number_part))

CondensedIso.sort(key = custom_sort_key)

#create an empty list of list
CondensedData = [ [] for _ in range(len(CondensedIso))]

for k, iso in enumerate(CondensedIso):
  for index, iso2 in enumerate(IsoList):
    if iso == iso2 :
      CondensedData[k] += BigExList[index]

DF = pandas.DataFrame()

for k, haha in enumerate(CondensedData):

  if len(haha) == 0: continue

  data = sorted(haha, key=lambda x: x[0])

  # print(header)
  # for row in data: print(row)
  df =  pandas.DataFrame(data, columns=header)

  #===== calculate excitation energy
  ex = []
  isoName = []

  for index, row in enumerate(data):
    isoName.append(CondensedIso[k])
    ex.append(row[0] - data[0][0])

  df.insert(1, 'Ex', ex)
  df.insert(0, 'Iso', isoName)

  DF = pandas.concat([DF, df])


for iso in CondensedIso:
  print("================ " +  iso)
  if 'Iso' in DF.columns:
    print(DF[DF['Iso']==iso])

#==== single-nucleus (or no SF data): nothing to transfer, stop here
if len(CondensedIso) < 2:
  print("Only one nucleus found; skipping transfer/ESPE sections.")
  sys.exit(0)
if len(SFData) == 0:
  print("No spectroscopic factor data found; skipping transfer/ESPE sections.")
  sys.exit(0)

#======================== re-arrange SF data 

#Fine out which iso is the middle, assume only 3 iso

# if len(CondensedIso) != 3 : exit()

#============== Removal Reaction
orbtialList = header[4:]
# print(orbtialList)

# print an SF table with a full '-' separator and the sum of each SF (orbital) column at the bottom,
# with the sum row aligned to the same columns as the table (the table itself is printed unchanged)
def printSFtable(df):
  s = df.to_string()
  lines = s.splitlines()
  print(s)
  width = max(len(line) for line in lines)
  print('-' * width)
  if df.empty: return
  def dec(t):
    return len(t.split('.')[1]) if '.' in t else 0
  def right_edge(k):
    # right edge of the k-th token (0-based) over all data rows
    r = 0
    for line in lines[1:]:
      pos = 0
      for i, tok in enumerate(line.split()):
        idx = line.index(tok, pos)
        if i == k:
          r = max(r, idx + len(tok))
          break
        pos = idx + len(tok)
    return r
  row = ' ' * width
  row = 'sum of SF' + row[len('sum of SF'):]
  for k, orb in enumerate(orbtialList):
    colk = 6 + k  # data rows: [index, BE, Ex, J, pi, T] + orbitals
    coldec = max(dec(line.split()[colk]) for line in lines[1:])
    sv = ('%.6f' % df[orb].sum()).rstrip('0').rstrip('.')
    d = max(coldec, dec(sv))
    txt = '%.*f' % (d, df[orb].sum())
    start = right_edge(colk) - len(txt)
    if start < 0: start = 0
    row = row[:start] + txt + row[start + len(txt):]
  print(row.rstrip())


df = DF[DF['Iso'] == CondensedIso[1]].sort_values(by=['BE'])
dft = DF[DF['Iso'] == CondensedIso[0]].sort_values(by=['BE'])

if df.empty or dft.empty:
  print("Missing energy data for the transfer nuclei; skipping transfer/ESPE sections.")
  sys.exit(0)

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

printSFtable(removeDF)

# if len(CondensedIso) < 3 : exit()
#===================== adding reaction

df = DF[DF['Iso'] == CondensedIso[0]].sort_values(by=['BE'])
dft = DF[DF['Iso'] == CondensedIso[1]].sort_values(by=['BE'])


adding = dfSF[dfSF['E2'] == round(df['BE'].iloc[0],3)]
# print(adding)
print("================================ %s adding to %s" % (CondensedIso[0], CondensedIso[1]))

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

printSFtable(addingDF)

#============================ Calcualte the ESPEs base on the transfer reaction

#TODO use real data to see the code is correct

print("====== Ground state eergy:")
groundStateBE = []
for iso in CondensedIso:
  d_iso = DF[DF['Iso']==iso]
  if d_iso.empty:
    print("Missing energy data for %s; skipping ESPE." % iso)
    sys.exit(0)
  BE = d_iso.sort_values(by='Ex')['BE'].iloc[0]
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
    Gp += row[orb]
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

  Epp = EpAvg + groundStateBE[1] - groundStateBE[0]
  Emm = -EmAvg - groundStateBE[0] + groundStateBE[1]

  print( " Ep' :  %.3f , Ep'' : %.3f, Gp : %.3f" % (EpAvg, Epp, Gp))
  print( " Em' :  %.3f , Em'' : %.3f, Gm : %.3f" % (EmAvg, Emm, Gm))

  maxOccpancy = int(orb.split('_')[1][0]) + 1
  print( " sum G :  %.3f, 2J+1 : %.0f, quenching : %.3f" % (Gp+Gm, maxOccpancy, (Gp+Gm)/maxOccpancy))

  print(" ESPS %s : %.3f" % (orb, (Epp * Gp + Emm * Gm)/(Gp + Gm)))





