#!/usr/bin/env /usr/bin/python3

import argparse, sys, os
import pandas as pd

# the API webpage
# https://www-nds.iaea.org/relnsd/vcharthtml/api_v0_guide.html#examples

# the service URL
livechart = "https://nds.iaea.org/relnsd/v0/data?"

import urllib.request

def lc_read_csv(url):
  req = urllib.request.Request(url)
  req.add_header('User-Agent', 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0')
  return pd.read_csv(urllib.request.urlopen(req))

haha = None
def Nuclides():
  global haha
  if haha is None:
    haha = lc_read_csv(livechart + 'fields=ground_states&nuclides=all')
  return haha

mp = 938.27208816; #MeV/c^2
mn = 939.56542052;

def FindZ(AZ):
  query = livechart + "fields=ground_states&nuclides=" + AZ
  temp = lc_read_csv(query)
  try :
    return temp['z']
  except :
    return 'na'
def FindSym(Z):
  try:
    haha = Nuclides()
    return (haha['symbol'][haha['z']==Z]).iloc[0]
  except:
    return 'na'
def Mass(A, Z):
  try :
    haha = Nuclides()
    BEA = float(haha['binding'][haha['z']==Z][haha['n']==(A-Z)].iloc[0])/1000
    return (A-Z)*mn + Z*mp - A * BEA
  except :
    return -404
def MassSym(AZ):
  query = livechart + "fields=ground_states&nuclides=" + AZ
  temp = lc_read_csv(query);
  Z = temp['z']
  N = temp['n']
  try :
    return Z*mp + N*mn - (Z+N)*temp['binding']/1000
  except:
    return -404
def Sp(A,Z,a,z):
  mA = Mass(A,Z)
  mB = Mass(A-a, Z-z)
  if z == 0 :
    mb = a * mn
  elif a == z :
    mb = a * mp
  else :
    mb = Mass(a,z)
  if (mB == -404 or mb == -404 or mA == -404) :
    return -404
  else:
    return mB + mb - mA
def Ex(AZ, maxMeV):
  query = livechart + "fields=levels&nuclides=" + AZ
  tempEx = lc_read_csv(query);
  try :
    return tempEx[['energy', 'jp']][tempEx['energy']<= maxMeV * 1000]
  except:
    return -404

def Info(AZ):
  query = livechart + "fields=ground_states&nuclides=" + AZ
  temp = lc_read_csv(query);
  print("============================================== ", AZ)
  try :
    Z = temp['z'][0]
    N = temp['n'][0]
    mass = Z*mp + N*mn - (Z+N)*temp['binding']/1000
    halfLife = temp['half_life_sec'][0]
    print("  A : %3d, Z : %3d, N : %3d, Mass : %.4f MeV" % (Z+N, Z, N, mass))
    print("Jpi : %3s,    half-live : %s sec" % (temp['jp'][0], halfLife))
    print("Sn  : %8.3f MeV, Sp  : %8.3f MeV" % (Sp(Z+N,Z, 1, 0), Sp(Z+N,Z, 1, 1)))
    print("S2n : %8.3f MeV, S2p : %8.3f MeV, Sd : %8.3f MeV" % (Sp(Z+N,Z, 2, 0), Sp(Z+N,Z, 2, 2), Sp(Z+N, Z, 2, 1)))
    print("S3n : %8.3f MeV, S3p : %8.3f MeV, St : %8.3f MeV, S(3He) : %8.3f MeV" % (Sp(Z+N,Z, 3, 0), Sp(Z+N,Z, 3, 3), Sp(Z+N, Z, 3, 1), Sp(Z+N, Z, 3, 2)))
    print("S4n : %8.3f MeV, S4p : %8.3f MeV, Sa : %8.3f MeV" % (Sp(Z+N,Z, 4, 0), Sp(Z+N,Z, 4, 4), Sp(Z+N, Z, 4, 2)))
    print("   magnetic dipole : ", temp['magnetic_dipole'][0], " mu N")
    print("electric quadruple : ", temp['electric_quadrupole'][0], " barn")
    if halfLife > 0 :
      print('------------ decay mode:')
      for i in range(1, 4) :
        print("%5s  %s %%" % (temp["decay_%d" % i][0], temp["decay_%d_%%" % i][0]))
      print('-------------------------')
  except :
    print("")
  print("====================================================")


fontSize = 20
plotHeight = 500
plotWidth = 400
yMin = -1
def DrawLevelsFromData(data, name, maxEx, xShift):
  import plotly.graph_objects as go
  ex=data['energy']
  jp=data['jp']
  fig = go.Figure()
  fig.update_layout(plot_bgcolor='white', width=plotWidth, height = plotHeight, margin=dict(l=0, r=0, t=0, b=0))
  fig.update_layout(showlegend=False)
  fig.update_xaxes(showline=False,  visible= False, range=[-1, 3])
  fig.update_yaxes(showline=True,  visible= True, range=[yMin, maxEx+2])

  l=ex.last_valid_index()

  fontSizeMeV=fontSize/plotHeight*(maxEx+1-yMin)
  #print(fontSizeMeV)
  #adjust text label y-pos
  ypos = ex.copy()

  noOverlap = False
  loop = 0

  while noOverlap == False and loop < 2*l :
    #print("================= %d" % loop)
    for i in range(1, l+1) :
      diff = ypos[i] - ypos[i-1]
      #print("%2d | %.3f, %.3f | %.4f" % (i, ypos[i], ypos[i-1], diff))
      if diff < fontSizeMeV :
        ypos[i-1] += (diff - fontSizeMeV)/2
        ypos[i] += (fontSizeMeV - diff)/2
        if( ypos[i-1] < yMin + fontSizeMeV/2) :
          ypos[i-1] = yMin + fontSizeMeV/2
          ypos[i] = ypos[i-1] + fontSizeMeV
      #print("   | %.3f, %.3f" % (ypos[i], ypos[i-1]))

    #print(ypos)
    ###=======inspection
    count = 0
    for i in range(1, l+1) :
      diff = ypos[i] - ypos[i-1]
      if diff > fontSizeMeV :
        count = count +1

    if count == l :
      noOverlap = True

    loop += 1

  for i in range(0,l+1):
    fig.add_trace(go.Scatter(x=[0 + xShift,1 + xShift], y=[ex[i],ex[i]],mode='lines',line=dict(color='black', width=1)))
    fig.add_trace(go.Scatter(x=[1.03 + xShift,1.1 + xShift, 1.19 + xShift], y=[ex[i],ypos[i],ypos[i]],mode='lines',line=dict(color='gray', width=1)))
    fig.add_annotation(x=1.2 + xShift, y=ypos[i], text=("%.3f, %s" % (ex[i], jp[i])), xanchor='left', font=dict(size=fontSize), showarrow=False)

  fig.add_annotation(x=0.5 + xShift, y=-0.5, text=name, font=dict(size=1.5*fontSize), showarrow=False)

  return fig

def DrawLevels(AZ, maxEx, xShift):
  import plotly.graph_objects as go
  query = livechart + "fields=ground_states&nuclides=" + AZ
  temp = lc_read_csv(query)
  Z = temp['z'][0]
  N = temp['n'][0]
  Sym = temp['symbol'][0]
  A = Z + N
  sn = Sp(A, Z, 1, 0)
  sp = Sp(A, Z, 1, 1)

  jaja=Ex(AZ, maxEx)
  jaja['energy'] = jaja['energy']/1000
  fig = DrawLevelsFromData(jaja, "Exp", maxEx, xShift)

  fig.add_annotation(x=0.5 + xShift, y=maxEx+1, text=("<sup>%s</sup>%s" % (A, Sym)), font=dict(size=1.5*fontSize, family="Times New Roman"), showarrow=False)

  if( sn < maxEx ):
    fig.add_trace(go.Scatter(x=[-0.6 + xShift,-0.1 + xShift], y=[sn,sn],mode='lines',line=dict(color='red', width=1)))
    fig.add_annotation(x=-0.6 + xShift, y=sn, text=("Sn %.3f" % sn), xanchor='left', yanchor='bottom', font=dict(size=fontSize, color='red'), showarrow=False)
  if( sp < maxEx ):
    fig.add_trace(go.Scatter(x=[-0.6 + xShift,-0.1 + xShift], y=[sp,sp],mode='lines',line=dict(color='blue', width=1)))
    fig.add_annotation(x=-0.6 + xShift, y=sp, text=("Sp %.3f" % sp), xanchor='left', yanchor='bottom', font=dict(size=fontSize, color='blue'), showarrow=False)

  return fig


#=====================================================
#=====================================================
# Shell-model levels, read from a kshell summary_*.txt "Energy levels" table:
#   N    J prty N_Jp    T     E(MeV)  Ex(MeV)  log-file
def LabelFromFileName(fileName):
  import re
  base = os.path.splitext(os.path.basename(fileName))[0]
  m = re.search(r'([A-Za-z]+\d+|\d+[A-Za-z]+)', base)
  if m:
    return m.group(1)
  return base

def ReadShellLevels(fileName):
  rows = []
  with open(fileName) as f:
    for line in f:
      cols = line.split()
      if len(cols) != 8:
        continue
      try:
        N = int(cols[0])
      except ValueError:
        continue
      rows.append({'energy': float(cols[6]), 'jp': cols[1] + cols[2]})
  if len(rows) == 0:
    print("Error: no energy levels found in %s" % fileName)
    sys.exit(1)
  df = pd.DataFrame(rows)
  df = df.sort_values(by='energy').reset_index(drop=True)
  return df

def DrawShell(fileName, label, maxEx, xShift):
  df = ReadShellLevels(fileName)
  df = df[df['energy'] <= maxEx].reset_index(drop=True)
  fig = DrawLevelsFromData(df, "SM", maxEx, xShift)
  fig.add_annotation(x=0.5 + xShift, y=maxEx+1, text=label, font=dict(size=1.5*fontSize, family="Times New Roman"), showarrow=False)
  return fig

def BuildPanels(nuclides, shellFiles):
  if nuclides and shellFiles and len(nuclides) != len(shellFiles):
    print("Error: number of --shell files (%d) must match number of nuclides (%d)." % (len(shellFiles), len(nuclides)))
    sys.exit(1)
  panels = []
  n = max(len(nuclides), len(shellFiles))
  for i in range(n):
    if i < len(nuclides):
      panels.append(('exp', nuclides[i], nuclides[i]))
    if i < len(shellFiles):
      label = nuclides[i] if i < len(nuclides) else LabelFromFileName(shellFiles[i])
      panels.append(('shell', shellFiles[i], label))
  return panels

def PlotLevels(nuclides, shellFiles, maxEx, outName):
  import plotly.offline as pyo
  panels = BuildPanels(nuclides, shellFiles)
  fig = None

  for i, (kind, key, label) in enumerate(panels):
    if kind == 'exp':
      f2 = DrawLevels(key, maxEx, 3*i)
    else:
      f2 = DrawShell(key, label, maxEx, 3*i)
    if fig is None:
      fig = f2
    else:
      for trace in f2.data:
        fig.add_trace(trace)
      for annotation in f2.layout.annotations:
        fig.add_annotation(annotation)

  fig.update_xaxes(showline=False,  visible= False, range=[-1, 3*len(panels)])
  fig.update_layout(width=plotWidth*len(panels))

  #fig.show() # assigne random port

  pyo.plot(fig, filename=outName, auto_open=True)

def build_parser():
  parser = argparse.ArgumentParser(
    description="Plot experimental level schemes (IAEA LiveChart) and/or shell model levels (kshell summary files).",
    epilog="Examples:\n"
           "  PlotEx.py 11C 12C 13C --max 7 --out temp.html\n"
           "  PlotEx.py 11C 12C 13C --shell summary_C11_ysox.txt summary_C12_ysox.txt summary_C13_ysox.txt\n"
           "  PlotEx.py --shell summary_C11_ysox.txt",
    formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("nuclides", nargs="*", metavar="AZ",
    help="experimental nuclide(s) to plot, e.g. 11C 12C 13C")
  parser.add_argument("--shell", dest="shellFiles", nargs="+", metavar="FILE",
    help="shell model summary file(s) with the 'Energy levels' table; plotted as SM panel(s), in the same order as the nuclides")
  parser.add_argument("--max", dest="maxEx", type=float, default=7,
    help="max excitation energy in MeV (default: 7)")
  parser.add_argument("--out", dest="outName", default="temp.html",
    help="output html file (default: temp.html)")
  return parser

if __name__ == "__main__":
  parser = build_parser()
  args = parser.parse_args()
  shellFiles = args.shellFiles if args.shellFiles else []
  if len(args.nuclides) == 0 and len(shellFiles) == 0:
    parser.print_help()
    sys.exit(0)
  PlotLevels(args.nuclides, shellFiles, args.maxEx, args.outName)


