#!/usr/bin/env /usr/bin/python3

#================== Plot the spectroscopic factors (C^2*S) from kshell transit
#================== log files against the shell-model level scheme of the target
#================== nucleus, in the same style as PlotEx.py.
#
# The folder containing the log_* and summary_*.txt files is given as the first
# argument (run with -h for usage).

import argparse, sys, os, re
import pandas as pd

fontSize = 20
plotHeight = 500
plotWidth = 400
yMin = -1
yAxisBottom = yMin - 0.7  # extra room below the SF ruler at yMin

orbColors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
             '#8c564b', '#e377c2', '#17becf', '#bcbd22', '#7f7f7f']

#=====================================================
#=====================================================
# shell-model levels, read from a kshell summary_*.txt "Energy levels" table:
#   N    J prty N_Jp    T     E(MeV)  Ex(MeV)  log-file
def FindSummary(isotope, folder):
  for f in sorted(os.listdir(folder)):
    if not f.startswith('summary_'):
      continue
    parts = os.path.splitext(f)[0].split('_')
    if len(parts) > 1 and parts[1] == isotope:
      return os.path.join(folder, f)
  return None

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
      rows.append({'E': float(cols[5]), 'Ex': float(cols[6]), 'jp': cols[1] + cols[2]})
  if len(rows) == 0:
    print("Error: no energy levels found in %s" % fileName)
    sys.exit(1)
  df = pd.DataFrame(rows)
  df = df.sort_values(by='Ex').reset_index(drop=True)
  return df

#=====================================================
#=====================================================
# transit (transfer) log:
#   fn_load_wave_l = <target>_*.wav        (heavier nucleus, "final" side)
#   fn_load_wave_r = <donor>_*.wav         (lighter nucleus, "initial" side)
#   orbit :   0  1  3  1    n 0p_3/2
#   2xJf      Ef      2xJi     Ei       Ex       C^2*S
#   1(   1)  -98.851 2(   1)  -86.884  11.967    1.2408
def ReadTransitSF(fileName):
  nucleusLeft = ""
  nucleusRight = ""
  particle = "n"
  orbital = ""
  rows = []
  with open(fileName) as f:
    for line in f:
      columns = line.strip().split()
      if 'fn_load_wave_l' in line and len(columns) >= 3:
        nucleusLeft = (columns[2].split('_'))[0]
      if 'fn_load_wave_r' in line and len(columns) >= 3:
        nucleusRight = (columns[2].split('_'))[0]
      if '.wav' in line and '=' in line:
        m = re.search(r'\+\s*([A-Za-z0-9]+)\s*$', line.strip())
        if m:
          particle = m.group(1)
      if len(columns) >= 7 and columns[0] == 'orbit':
        orbital = columns[6] + "-" + columns[7]
      if len(columns) == 8 and len(columns[0]) > 1 and columns[0][1] == '(':
        rows.append({'orb': orbital,
                     'Ef': float(columns[2]),
                     'Ei': float(columns[5]),
                     'SF': float(columns[7])})
  return nucleusLeft, nucleusRight, particle, rows

def MassOf(nucleus):
  m = re.match(r'([A-Za-z]+)(\d+)', nucleus)
  if m:
    return (m.group(1), int(m.group(2)))
  return (nucleus, -1)

def NuclideLabel(nucleus):
  m = re.match(r'([A-Za-z]+)(\d+)', nucleus)
  if m:
    return ("<sup>%s</sup>%s" % (m.group(2), m.group(1)))
  return nucleus

def OrbKey(orb):
  m = re.match(r'([np])-(\d+)([a-z])_(\d+)/(\d+)', orb)
  if m:
    lorder = 'spdfgh'
    l = m.group(3)
    return (m.group(1), int(m.group(2)), lorder.index(l) if l in lorder else 9, int(m.group(4)))
  return ('z', 99, 99, 99)

#=====================================================
#=====================================================
# one reaction: target (left) built from donor (right) ground state + particle
def BuildReaction(left, right, particle, rows, folder):
  targetFile = FindSummary(left, folder)
  if targetFile is None:
    print("Error: no summary_%s_*.txt found in %s" % (left, folder))
    sys.exit(1)
  tdf = ReadShellLevels(targetFile)

  donorFile = FindSummary(right, folder)
  if donorFile is not None:
    EiGs = ReadShellLevels(donorFile)['E'].min()
    print("%s -> %s : donor %s ground state E = %.3f MeV" % (right, left, right, EiGs))
  else:
    EiGs = min(r['Ei'] for r in rows)
    print("Warning: no summary for %s, donor ground state E = %.3f MeV from log" % (right, EiGs))

  levels = []
  for idx, row in tdf.iterrows():
    levels.append({'E': row['E'], 'Ex': row['Ex'], 'jp': row['jp'], 'orbs': {}})

  nMatched = 0
  for r in rows:
    if abs(r['Ei'] - EiGs) > 0.01 or r['SF'] <= 0:
      continue
    best = None
    bestD = 1e9
    for lv in levels:
      d = abs(lv['E'] - r['Ef'])
      if d < bestD:
        bestD = d
        best = lv
    if best is not None and bestD < 0.05:
      best['orbs'][r['orb']] = best['orbs'].get(r['orb'], 0.0) + r['SF']
      nMatched += 1

  print("%s -> %s : %d SF rows matched to %s levels" % (right, left, nMatched, left))
  return {'left': left, 'right': right, 'particle': particle, 'levels': levels}

#=====================================================
#=====================================================
def niceStep(x):
  import math
  p = 10 ** math.floor(math.log10(x))
  for m in (1, 2, 2.5, 5, 10):
    if m * p >= x:
      return m * p
  return 10 * p

def DrawSFPanel(react, maxEx, xShift, orbColorMap, orbOffset):
  import plotly.graph_objects as go
  levels = [lv for lv in react['levels'] if lv['Ex'] <= maxEx and lv['orbs']]
  orbs = [o for o in orbColorMap if any(o in lv['orbs'] for lv in levels)]
  maxSF = 0.0
  for lv in levels:
    for o, s in lv['orbs'].items():
      maxSF = max(maxSF, s)
  if maxSF <= 0:
    maxSF = 1.0

  fig = go.Figure()
  fig.update_layout(plot_bgcolor='white', width=plotWidth, height=plotHeight, margin=dict(l=0, r=0, t=0, b=0))
  fig.update_layout(showlegend=False)
  fig.update_xaxes(showline=False,  visible= False, range=[-1, 3])
  fig.update_yaxes(showline=True,  visible= True, range=[yAxisBottom, maxEx+2])

  ex = pd.Series([lv['Ex'] for lv in levels])
  l = len(ex) - 1

  fontSizeMeV = fontSize/plotHeight*(maxEx+2-yAxisBottom)
  ypos = ex.copy()
  noOverlap = False
  loop = 0
  while noOverlap == False and loop < 2*l :
    for i in range(1, l+1) :
      diff = ypos[i] - ypos[i-1]
      if diff < fontSizeMeV :
        ypos[i-1] += (diff - fontSizeMeV)/2
        ypos[i] += (fontSizeMeV - diff)/2
        if( ypos[i-1] < yAxisBottom + fontSizeMeV/2) :
          ypos[i-1] = yAxisBottom + fontSizeMeV/2
          ypos[i] = ypos[i-1] + fontSizeMeV
    count = 0
    for i in range(1, l+1) :
      if (ypos[i] - ypos[i-1]) > fontSizeMeV :
        count = count + 1
    if count == l :
      noOverlap = True
    loop += 1

  legendDone = {}
  for i, lv in enumerate(levels):
    # level label connector + label
    fig.add_trace(go.Scatter(x=[1.03 + xShift,1.1 + xShift, 1.19 + xShift], y=[lv['Ex'],ypos[i],ypos[i]],mode='lines',line=dict(color='gray', width=1),showlegend=False))
    fig.add_annotation(x=1.2 + xShift, y=ypos[i], text=("%.3f, %s" % (lv['Ex'], lv['jp'])), xanchor='left', font=dict(size=fontSize), showarrow=False)

    nOrbs = len(orbs)
    for j, o in enumerate(orbs):
      if o not in lv['orbs']:
        continue
      sf = lv['orbs'][o]
      y = lv['Ex'] + (j - (nOrbs-1)/2)*orbOffset
      x1 = 0.9 * sf / maxSF
      color = orbColorMap[o]
      show = (o not in legendDone)
      legendDone[o] = True
      fig.add_trace(go.Scatter(x=[0 + xShift, x1 + xShift], y=[y, y], mode='lines',
                               line=dict(color=color, width=4),
                               legendgroup=o, name=o, showlegend=show,
                               hovertext=["%s  SF=%.4f" % (o, sf)]*2, hoverinfo='text'))

  # SF ruler at yMin: bar length = 0.9 * SF / maxSF, so SF = x/0.9 * maxSF
  step = niceStep(maxSF/4)
  ticks = []
  v = 0.0
  while v <= maxSF + 1e-9:
    ticks.append(v)
    v += step
  fig.add_trace(go.Scatter(x=[0 + xShift, 0.9 + xShift], y=[yMin, yMin], mode='lines', line=dict(color='black', width=1), showlegend=False))
  for v in ticks:
    xt = 0.9 * v / maxSF
    fig.add_trace(go.Scatter(x=[xt + xShift, xt + xShift], y=[yMin, yMin+0.12], mode='lines', line=dict(color='black', width=1), showlegend=False))
    fig.add_annotation(x=xt + xShift, y=yMin-0.12, text=('%g' % v), xanchor='center', yanchor='top', font=dict(size=14), showarrow=False)
  fig.add_annotation(x=-0.08 + xShift, y=yMin, text="C<sup>2</sup>*S", xanchor='right', yanchor='bottom', font=dict(size=14), showarrow=False)

  fig.add_annotation(x=0.5 + xShift, y=maxEx+1, text=NuclideLabel(react['left']), font=dict(size=1.5*fontSize, family="Times New Roman"), showarrow=False)
  fig.add_annotation(x=0.5 + xShift, y=-0.5, text=("%s + %s" % (react['right'], react['particle'])), font=dict(size=fontSize), showarrow=False)

  return fig

#=====================================================
#=====================================================
# legend column: one colored line sample + orbital name per orbital
def DrawOrbitalPanel(orbList, orbColorMap, maxEx, xShift):
  import plotly.graph_objects as go
  fig = go.Figure()
  fig.update_layout(plot_bgcolor='white', width=plotWidth, height=plotHeight, margin=dict(l=0, r=0, t=0, b=0))
  fig.update_layout(showlegend=False)
  fig.update_xaxes(showline=False,  visible= False, range=[-1, 3])
  fig.update_yaxes(showline=True,  visible= True, range=[yAxisBottom, maxEx+2])

  n = len(orbList)
  spacing = (maxEx - yMin) / (n + 1)
  for i, o in enumerate(orbList):
    y = yMin + (i+1)*spacing
    fig.add_trace(go.Scatter(x=[0 + xShift, 0.25 + xShift], y=[y, y], mode='lines',
                             line=dict(color=orbColorMap[o], width=4), showlegend=False))
    fig.add_annotation(x=0.32 + xShift, y=y, text=o, xanchor='left', font=dict(size=fontSize), showarrow=False)

  fig.add_annotation(x=0.5 + xShift, y=maxEx+1, text="orbital", font=dict(size=1.5*fontSize, family="Times New Roman"), showarrow=False)

  return fig

#=====================================================
#=====================================================
def PlotSF(folder, maxEx, outName, orbOffset):
  try:
    import plotly.offline as pyo
  except ModuleNotFoundError:
    print("Error: 'plotly' is not installed for this python interpreter.")
    print("Please use the venv, e.g.:")
    print("    /home/ryan/venv/bin/python PlotSF.py %s --out %s" % (folder, outName))
    print("or activate it first:  source /home/ryan/venv/bin/activate")
    sys.exit(1)

  all_files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]
  fileList = [f for f in all_files if 'log' in f]
  fileList = [item for item in fileList if "summary" not in item]
  fileList = [item for item in fileList if ".py" not in item]

  reactions = {}
  for f in fileList:
    path = os.path.join(folder, f)
    try:
      with open(path, 'r') as file:
        content = file.read()
    except (FileNotFoundError, IOError):
      continue
    if 'Spectroscopic factor' not in content:
      continue
    print("!" + f)
    left, right, particle, rows = ReadTransitSF(path)
    if not left or not right or len(rows) == 0:
      continue
    if (left, right) not in reactions:
      reactions[(left, right)] = (particle, [])
    reactions[(left, right)][1].extend(rows)

  if len(reactions) == 0:
    print("No transit (Spectroscopic factor) log files found in %s" % folder)
    sys.exit(1)

  pairs = sorted(reactions.keys(), key=lambda pr: MassOf(pr[0]))

  built = []
  for (left, right) in pairs:
    particle, rows = reactions[(left, right)]
    built.append(BuildReaction(left, right, particle, rows, folder))

  orbList = sorted(set(o for r in built for lv in r['levels'] for o in lv['orbs']), key=OrbKey)
  orbColorMap = {o: orbColors[i % len(orbColors)] for i, o in enumerate(orbList)}
  print("orbitals: %s" % ", ".join(orbList))

  fig = None
  for i, react in enumerate(built):
    f2 = DrawSFPanel(react, maxEx, 3*i, orbColorMap, orbOffset)
    if fig is None:
      fig = f2
    else:
      for trace in f2.data:
        fig.add_trace(trace)
      for annotation in f2.layout.annotations:
        fig.add_annotation(annotation)

  if orbList:
    f2 = DrawOrbitalPanel(orbList, orbColorMap, maxEx, 3*len(built))
    if fig is None:
      fig = f2
    else:
      for trace in f2.data:
        fig.add_trace(trace)
      for annotation in f2.layout.annotations:
        fig.add_annotation(annotation)
    nPanels = len(built) + 1
  else:
    nPanels = len(built)

  fig.update_xaxes(showline=False,  visible= False, range=[-1, 3*nPanels])
  fig.update_layout(width=plotWidth*nPanels)
  fig.update_layout(showlegend=False)

  pyo.plot(fig, filename=outName, auto_open=True)

def build_parser():
  parser = argparse.ArgumentParser(
    description="Plot spectroscopic factors (C^2*S) from kshell transit logs "
                "against the shell-model level scheme of the target nucleus.",
    epilog="Usage notes:\n"
           "  * Any file in the folder whose name contains 'log' (and does not\n"
           "    contain 'summary', and is not a .py file) is inspected; files\n"
           "    containing 'Spectroscopic factor' are transit logs.\n"
           "  * For each transit log, the left (heavier) nucleus is the target\n"
           "    and the right (lighter) nucleus is the donor. Only SF rows from\n"
           "    the donor ground state are plotted, and each final-state energy\n"
           "    Ef is matched to the closest target level in summary_<target>_*.txt.\n"
           "  * One panel per reaction, one colored line per orbital; line length\n"
           "    is proportional to the SF.",
    formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("folder", help="folder containing the log files")
  parser.add_argument("--max", dest="maxEx", type=float, default=12,
    help="max excitation energy in MeV (default: 12)")
  parser.add_argument("--offset", dest="orbOffset", type=float, default=0.05,
    help="vertical offset between orbital lines in MeV (default: 0.05)")
  parser.add_argument("--out", dest="outName", default="temp_sf.html",
    help="output html file (default: temp_sf.html)")
  return parser

if __name__ == "__main__":
  parser = build_parser()
  if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(0)
  args = parser.parse_args()
  PlotSF(args.folder, args.maxEx, args.outName, args.orbOffset)
