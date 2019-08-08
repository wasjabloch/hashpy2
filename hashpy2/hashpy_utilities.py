#! /usr/bin/python3
from subprocess import Popen, PIPE
from struct import unpack
from glob import glob
from warnings import warn
from copy import deepcopy
import numpy as np
from obspy import read_events, read
from obspy.signal.util import util_geo_km
from obspy.signal.trigger import recursive_sta_lta, pk_baer
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Cursor
from pyrocko import moment_tensor as mtt


def MakeEvent(sID, hypdir):
    year = '20' + str(sID[0:2])
    sIDnew = sID[0:6]+'.'+sID[6:12]
    print(sIDnew)
    try:
        hypfile = glob(hypdir + '*'+sIDnew+'*.qml.hyp')[0]
    except:
        print(hypdir + '*'+sIDnew+'*.hyp')
        print(glob(hypdir + '*'+sIDnew+'*.hyp'))
        oldhypfile = glob(hypdir + '*'+sIDnew+'*.hyp')[0]
        hypfile = oldhypfile.replace('.hyp', '.qml.hyp')
        FixBrokenHyp(oldhypfile, hypfile)
    event = read_events(hypfile)[0]
    with open(hypfile, 'r') as hf:
        for line in hf:
            if line.startswith('GEOGRAPHIC'):
                event.origins[0].latitude = float(line.split()[9])
                event.origins[0].longitude = float(line.split()[11])
    return event, hypfile


def FixBrokenHyp(oldhypfile, newhypfile):
    """ Open NLL hypfile and inster QML keywords, if absent"""
    qmlline = []
    with open(oldhypfile, 'r') as hf:
        data = hf.read()
    lines = data.splitlines()
    for j in range(len(lines)):
      lines[j]+='\n'
    lined = dict([line.split(None, 1) for line in lines[0:17]])
    assocPhCt = usedPhCt = str(0)
    assocStaCt = usedStaCt =  depthPhCt = str(0)
    stdErr = str(0)
    azGap = secAzGap = str(0)
    gtLevel = str(0)
    minDist = maxDist = medDist = str(0)
    horUnc = minHorUnc = maxHorUnc = azMaxHorUnc = str(0)
    try:
        lined['QML_OriginQuality']
    except:
        qmlline = 'QML_OriginQuality ' + 'assocPhCt ' + assocPhCt + ' usedPhCt ' +  usedPhCt + ' assocStaCt ' +\
                  assocStaCt  + ' usedStaCt ' + usedStaCt  + ' depthPhCt ' + depthPhCt  + ' stdErr ' +\
                  stdErr  + ' azGap ' + azGap  + ' secAzGap ' + secAzGap  + ' gtLevel ' + gtLevel + ' minDist ' +\
                  minDist + ' maxDist ' + maxDist + ' medDist ' + medDist + '\n'
        lines.insert(13, qmlline)
    try:
        lined['QML_OriginUncertainty']
    except:
        qmlline = 'QML_OriginUncertainty' + ' horUnc ' + horUnc + ' minHorUnc ' +  minHorUnc + ' maxHorUnc ' +  maxHorUnc + ' azMaxHorUnc ' + azMaxHorUnc + '\n'
        lines.insert(13, qmlline)
    if qmlline:
        with open(newhypfile, 'w') as hf:
            for line in lines:
                hf.write(line)


def MakeTraces(sID, wvdir):
    year = '20' + str(sID[0:2])
    wavedir = wvdir + year + '/' + sID 
    traces = read(wavedir + '/*HH*')
    return traces


def dist_calc(lon1, lat1, lon2, lat2):
    lon1 = np.radians(lon1)
    lon2 = np.radians(lon2)
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2.)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.)**2
    r = 6371.
    distance = r*2.*np.arcsin(np.sqrt(a))

    return distance


def dist2stations(stations, event, filename):
    stationdict = {}
    with open(filename,'r') as stationfile:
        for line in stationfile:
            line = line.rsplit()
            station = line[1]

            if station in stations:
                evLat = np.float(event.origins[0].latitude)
                evLon = np.float(event.origins[0].longitude)
                lat = np.float(line[3])
                lon = np.float(line[4])

                dist = dist_calc(evLon, evLat, lon, lat)
                stationdict[station] = dist 

    return stationdict


def AskForGrids(directory):
    griddirs = directory+'data/nlloc/ttables/ttt'
    griddirs = griddirs.split(" ")
    return griddirs


def ReadNLLAngle(station, lon, lat, z, grid):
    """
    Reads a NonLinLoc angle files and returns take-off angles at specific point
    Input:
    basename: NLL angle file
    station: station name
    lon: longitude...
    lat: latitude of the event 
    z: depth positive down
    Output:
    azi: Azimut clockwise from North
    dip: Dip, 0.0 (down) to (180.0) up
    qual: quality, 0 (low) to 10 (high)
    """
    
    hfile = grid
    bfile = grid[:-3] + 'buf'

    with open(hfile, 'r') as hf:
        for l, line in enumerate(hf):
            if l==0:
                nx, ny, nz, x0, y0, z0, d  = list(map(float, line.split(None)[0:7]))
            if line.startswith('TRANSFORM'):
                lat0 = float(line.split()[5])
                lon0 = float(line.split()[7])

    x, y = util_geo_km(lon0, lat0, lon, lat)

    xs = np.arange(x0, nx*d, d)
    ys = np.arange(y0, ny*d, d)
    zs = np.arange(z0, nz*d, d)

    for dim in ['x', 'y', 'z']:
        dims = dim + 's'
        ma = max(eval(dims))
        mi = min(eval(dims))
        if eval(dim) < mi or eval(dim) > ma:
            msg = dim + ' out of range! (is ' + str(eval(dim)) + ', but range is ' + str(mi) + ' to ' + str(ma) + ')'
            warn(msg)

    xx = (abs(xs-x)).argmin()
    yy = (abs(ys-y)).argmin()
    zz = (abs(zs-z)).argmin()

    i = int(((xx*ny + yy)*nz + zz)*4)

    with open(bfile, 'rb') as buf:
        buf.seek(i)
        fval = buf.read(4) # contains quality, azimuth and dip. 4 = short int

    qualdip, azi = unpack('hh', fval)
    azi = azi / 10.
    dip = (qualdip // 16) / 10.
    qual = qualdip % 16 # see http://alomax.free.fr/nlloc/ -> Formats

    return azi, dip, qual


# Picker
################################

def PickRequest(traces, station, event, dist):
    for pick in event.picks:
        if (pick.waveform_id['station_code'] == station) and (pick.phase_hint == 'P'):
            pickvalid = True
            picktime = pick.time
            hint = 'b'
            pol = ''
            break
        else: 
            pickvalid = False
    if not pickvalid:
        picktime, pol = AutoPicker(traces,station, event, dist)
        hint = 'a'
    pick = [picktime, hint, pol]

    return pick


def AutoPicker(traces, station, event, dist, component = 'HHZ'):
    starttime, endtime = GetTimeWindow(event, station) #.filter('highpass', freq=1.0)
    tr = traces.select(station=station, channel=component).slice(starttime=starttime, endtime=endtime).detrend('demean')
    if len(tr) == 0:
        return 0, ''
        
    df = tr[0].stats.sampling_rate
    p_pick, phase_info = pk_baer(tr[0].data,df,20,60,7.0,12.0,100,100) # default parameters
    picktime = starttime + (p_pick / df) 
    print('Phaseinfo', picktime, phase_info)
    if phase_info == [] or phase_info == '':
        pol = ''
    else:
        pol = phase_info[2] 

    return picktime, pol


def PolarityPicker(picktime, trace, gs):
    poltrace = trace.copy()
    fbmode = 'gradient'

    pol, poltime = PolarityStd(picktime, poltrace.filter('highpass', freq=1.0), gs[3,1])
    polhint = 'Std'
    poltrace = trace.copy()
    if pol in ['+', '-']:
        PolarityGradient(picktime, poltrace, gs[3,0])
        # PolarityFB(picktime, poltrace, gs, mode=fbmode)
    else:
        pol, poltime = PolarityGradient(picktime, poltrace, gs[3,0])
        polhint = 'Gradient'
        if pol in ['+', '-']:
            # PolarityFB(picktime, poltrace, gs, mode=fbmode)
            pass
        else:
            # pol, poltime = PolarityFB(picktime, poltrace, gs, mode=fbmode)
            # polhint = 'FB'
            pass

    return pol, poltime, polhint


def PolarityGradient(picktime, trace, gs, color='black'):
    print('GradientPicker')
    timewindow = 0.5
    poltrace = trace.copy()


    tr = poltrace.slice(starttime=picktime, endtime=picktime+timewindow)#.detrend('linear')
    try:
        starttime = tr.stats.starttime
    except:
        return '', 0

    tmpstd = np.std(tr.data)
    if np.isnan(tmpstd):
        pol = ''
        poltime = 0
        return '', poltime
 
    df = tr.stats.sampling_rate

    print(len(tr.data), tmpstd)
    ### simple gradient picker
    if False:
        for ii in range(len(tr.data)-2):
            d1 = tr.data[ii]
            d2 = tr.data[ii+1]
            d3 = tr.data[ii+2]

            if np.sign(d2-d1) != np.sign(d3-d2):

                poltime = tr.times()[ii+1]
                if np.sign(d2-d1) == 1.0:
                    pol = '+'
                elif np.sign(d2-d1) == -1.0:
                    pol = '-'
                break

            else:
                pol = ''
                poltime = 0

    ### advanced gradient picker, measuring the gradient change 
    ### and check if the gradient before the polarity pick has the same sign for some samples
    elif True:
        for ii in range(len(tr.data)-2):
            d1 = tr.data[ii]
            d2 = tr.data[ii+1]
            d3 = tr.data[ii+2]

            if abs(d2-d1) < abs(d3-d2):
                for nn in range(len(tr.data)-3-ii):
                    nn +=1
                    d1 = tr.data[ii+nn]
                    d2 = tr.data[ii+1+nn]
                    d3 = tr.data[ii+2+nn]
                    pol = ''
                    poltime = 0
                    if np.sign(d2-d1) != np.sign(d3-d2):
                        poltime = tr.times()[ii+1+nn]
                        if np.sign(d2-d1) == 1.0:
                            pol = '+'
                        elif np.sign(d2-d1) == -1.0:
                            pol = '-'
                        break

                if pol in ['+', '-']:
                    dpp2 = tr.data[ii+3+nn]
                    dpp1 = tr.data[ii+2+nn]
                    dp = tr.data[ii+1+nn]
                    dpm1 = tr.data[ii+nn]
                    dpm2 = tr.data[ii-1+nn]
                    if np.sign(dp-dpm1) != np.sign(dpm1-dpm2) and np.sign(dpp2-dpp1) != np.sign(dpp1-dp):
                        continue
                    else:
                        break
            else:
                pol = ''
                poltime = 0
        else:
            pol = ''
            poltime = 0


    ### gradient picker, the gradient has to be stable for few more samples
    else:
        # now it will look every sample after the intial polarity pick and check if the sign is reversed
        # another idea would be to check if the amplitude is not changing that much for the next samples...
        checkLength = 2
        for ii in range(len(tr.data)-2-checkLength):
            d1 = tr.data[ii]
            d2 = tr.data[ii+1]
            d3 = tr.data[ii+2]

            if np.sign(d2-d1) != np.sign(d3-d2):
                nn = 0
                while np.sign(d2-d1) != np.sign(tr.data[ii+2+nn]-tr.data[ii+1+nn]):
                    nn += 1
                    if nn > checkLength:
                        poltime = tr.times()[ii+1]
                        if np.sign(d2-d1) == 1.0:
                            pol = '+'
                        elif np.sign(d2-d1) == -1.0:
                            pol = '-'
 
                        break
                    else:
                        pol = ''
                        poltime = 0
                else:
                    pol = ''
                    poltime = 0

                if pol != '':
                    break
            else:
                pol = ''
                poltime = 0

    ax0 = plt.subplot(gs)
    ax0.set_title('Polarity Picker Gradient')
    ax0.plot(tr.times(), tr.data, color=color)
    ax0.axvline(x=picktime - starttime + poltime, color=color, zorder=100, label='Gr')
    ax0.plot(picktime - starttime + poltime, np.median(tr.data), GetMarker(pol), color=color, zorder=100)

    return pol, poltime


def PolarityStd(picktime, trace, gs, color='black'):
    print('StdPicker')
    ax0 = plt.subplot(gs)
    ax0.set_title('STD Picker')

    poltrace = trace.copy()
    noisestart = 1.5 
    noiseend = 0.5
    stdtr = poltrace.slice(starttime=picktime-noisestart, endtime=picktime-noiseend)
    std = np.std(stdtr.data)
    
    poltrace = trace.copy()
    signalstart = 0.5
    signalend = 1.5
    stdSigtr = poltrace.slice(starttime=picktime+signalstart, endtime=picktime+signalend)
    stdSig = np.std(stdSigtr.data)

    if np.isnan(std) or np.isnan(stdSig):
        pol = ''
        poltime = 0
        return '', poltime

    ratioStd = stdSig/std

    timewindow = 1.
    timebeforepick = 0.1
    poltrace = trace.copy()
    tr1 = poltrace.slice(starttime=picktime-timebeforepick, endtime=picktime+timewindow)

    starttime = tr1.stats.starttime
    df = tr1.stats.sampling_rate

    if len(tr1.data) == 0:
        return '', 0

    ax0.plot(tr1.times(), tr1.data, color='black')
    ax0.axvline(x=picktime - starttime, color='red', zorder=100, label='Pick')
    print(picktime, starttime, picktime- starttime)

    stdmin = 0.
    stdmax = 100

    if ratioStd > 150:
        minCnt = 15
        stdsinglef = 20
    elif ratioStd > 30:
        minCnt = 10
        stdsinglef = 8
    elif ratioStd > 20:
        minCnt = 7
        stdsinglef = 5
    else:
        minCnt = 3
        stdsinglef = 2
    print(minCnt, ratioStd, stdsinglef)


    ### minCounts
    stdfacs = np.arange(stdmin,stdmax,1)

    norm = mpl.colors.Normalize(vmin=stdmin, vmax=stdmax)
    cmap = mpl.cm.get_cmap('Spectral')

    colors = ['black', 'red', 'blue', 'green', 'gray']
    cnt = 0

    fpol = ''
    for stdfac in stdfacs:
        color = cmap(norm(stdfac)) 
        stdpos = tr1.data[0] + stdfac*std
        stdneg = tr1.data[0] - stdfac*std
        for ii in range(len(tr1.data)-1):   
            ii += 1
            if tr1.data[ii] > stdpos or tr1.data[ii] < stdneg :
                poltime = tr1.times()[ii]
                if np.sign(tr1.data[ii]-tr1.data[ii-1]) == 1.0:
                    pol = '+'
                elif np.sign(tr1.data[ii]-tr1.data[ii-1]) == -1.0:
                    pol = '-'
                break

            else:
                pol = ''
                poltime = 0
        
        if pol == fpol:
            cnt += 1       
        else:
            cnt = 0
            fpol = pol

        if fpol == '':
            fpol = pol


        if minCnt <= cnt:
            break

        ax0.axvline(x=poltime, color=color, zorder=10, label='STD')
        ax0.plot(poltime, 0, GetMarker(pol), color=color, zorder=10)
        ax0.axhline(y=stdpos, color=color, zorder=-100, alpha=0.5)
        ax0.axhline(y=stdneg, color=color, zorder=-100, alpha=0.5)

    if minCnt > cnt:
        pol = ''
        poltime = 0.

    ### stdsinglef
    stdpos = tr1.data[0] + stdsinglef*std
    stdneg = tr1.data[0] - stdsinglef*std
    for ii in range(len(tr1.data)-1):   
        ii += 1
        if tr1.data[ii] > stdpos or tr1.data[ii] < stdneg :
            poltime = tr1.times()[ii]
            if np.sign(tr1.data[ii]-tr1.data[ii-1]) == 1.0:
                pol = '+'
            elif np.sign(tr1.data[ii]-tr1.data[ii-1]) == -1.0:
                pol = '-'
            break

        else:
            pol = ''
            poltime = 0

    ax0.axvline(x=poltime, color='green', zorder=100, label='STD')
    ax0.plot(poltime, 0, GetMarker(pol), color='green', zorder=100)
    ax0.axhline(y=stdpos, color='green', zorder=-10, alpha=0.5)
    ax0.axhline(y=stdneg, color='green', zorder=-10, alpha=0.5)

    poltime = poltime-(picktime-starttime)

    return pol, poltime


def PolarityFB(picktime, trace, gs, mode='gradient'):
    ## really not tested enough

    freqs = [0.5, 1, 2, 3]
    colors = ['red', 'blue', 'green', 'gray']
    weights = [1.5,1,1,0.5]
    timewindow = 0.5
    pols = 0

    if mode == 'std':
        PolarityFunc = deepcopy(PolarityStd)
        gs = gs[3,1]

    elif mode == 'gradient':
        PolarityFunc = deepcopy(PolarityGradient)
        gs = gs[3,0]

        poltrace = trace.copy()
        tr = poltrace #.slice(starttime=picktime, endtime=picktime+timewindow)
        pol, poltime = PolarityFunc(picktime, tr, gs, 'black')
        
        if pol == '+':
            pols += 1.5*1.0
        elif pol == '-':
            pols += 1.5*-1.0
        else:
            pass

    else:
        print('False FB PolarityPicker Mode')
        return '', 0,
  
    for freq, color, weight in zip(freqs, colors, weights):
        poltrace = trace.copy()
        tr = poltrace.filter('highpass', freq=freq) #.slice(starttime=picktime, endtime=picktime+timewindow)
        pol, poltime = PolarityFunc(picktime, tr, gs, color)
        if pol == '+':
            pols += weight*1.0
        elif pol == '-':
            pols += weight*-1.0
        else:
            pass

    if np.sign(pols) == 1.0:
        pol = '+'
    elif np.sign(pols) == -1.0:
        pol = '-'
    
    ax0 = plt.subplot(gs)
    ax0.plot(0, 0, GetMarker(pol), markersize=10, color='purple', zorder=100)
    poltime = 0

    return pol, poltime


def GetMarker(pol):
    if pol == '+':
        return '^'
    elif pol == '-':
        return 'v'
    else:
        return ''


# Polarities
################################

def MakePolarityRequest(traces, station, event, pick, distance):
    isvalid = True

    pol = AskForPolarity(traces, station, event, pick, distance)
    print(pol)
    if pol is 'x':
        isvalid = False

    return pol, isvalid


def AskForPolarity(traces, station, event, pick, distance):
    fig, pol = ShowPolarity(traces, station, event, pick, distance)
    fig.show()

    ans = ''
    if pol == '':
        while ans not in ['+', '-', 'x']:
            ans = input('What P-wave polarity do you see? (+, -, x)')
    else:
        ans = input('Press \'a\' or ENTER for accepting the suggested Polarity: %s or correct it.' %pol)
        while ans not in ['+', '-', 'x', 'a', '']:
            if ans in ['a','']:
                break
            else:
                ans = input('Polarity must be one of +, -, x or a.\n')
    if ans in ['+', '-', 'x']:
        pol = ans
    plt.close(fig)

    return pol


def ShowPolarity(traces, station, event, pick, distance, component = 'HHZ'):
        
    picktime, hint, polPick = pick
    starttime, endtime = GetTimeWindow(event, station, distance)

    rawtr = deepcopy(traces.select(station=station, channel=component).slice(starttime=starttime, endtime=endtime).detrend('demean').integrate())
    tr = rawtr.copy()
    try:
        df = rawtr[0].stats.sampling_rate
        dt = rawtr[0].stats.starttime - rawtr[0].stats.endtime
    except:
        print('Something is broken')
        fig = plt.figure()
        return fig,'x'

    fig = plt.figure(figsize=(16, 9))
    wm = plt.get_current_fig_manager()
    wm.window.attributes('-topmost', 0)
    gs = gridspec.GridSpec(4,2)

    ## Polarity window
    pol, poltime, polhint = PolarityPicker(picktime,tr[0], gs)
    tr = rawtr.copy()

    ## Raw data
    ax1 = plt.subplot(gs[0,:])
    ax1.plot(tr[0].times(), tr[0].data, 'k-')
    ax1.axvline(x=picktime - starttime, color='red', zorder=100)
    plt.text(picktime - starttime + 0.25, 0.9*(max(tr[0].data)), '%s,%s' %(hint,polPick), color='red')
    ax1.set_title('Raw trace %s.%s \n at %.5s km distance ' %(station,component,distance))

    ## Zoomed data raw 
    startzoom = picktime - 1.
    endzoom = picktime + 1.5
    tr = rawtr.copy()
    trzoom = tr.slice(starttime=startzoom, endtime=endzoom).detrend('linear')
    try:
        tmp = trzoom[0]
    except:
        print('Something is broken')
        fig = plt.figure()
        return fig, 'x'

    ax3 = plt.subplot(gs[2,0])
    ax3.plot(trzoom[0].times(), trzoom[0].data, 'k-')
    ax3.axvline(x=picktime - startzoom, color='red', zorder=100)

    ax3.axvline(x=picktime - startzoom + poltime, color='blue', zorder=100)
    ax3.plot(picktime - startzoom + poltime, 0, GetMarker(pol), color='blue')
    plt.text(picktime - startzoom + poltime + 0.05, 0.9*(max(trzoom[0].data)), '%s' %(polhint), color='blue')

    ax3.set_title('Zoomed raw trace')

    ## Filtered data
    freq = 1.
    tr = rawtr.copy()
    trfil = tr.filter('highpass', freq=freq)
    ax2 = plt.subplot(gs[1,:])
    ax2.plot(trfil[0].times(), trfil[0].data, 'k-')
    ax2.axvline(x=picktime - starttime, color='red', zorder=100)
    ax2.set_title('Filtered trace, HP %s Hz' %freq)

    ## Zoomed data filtered
    trzoomfil = trfil.slice(starttime=startzoom, endtime=endzoom)
    ax4 = plt.subplot(gs[2,1])
    ax4.plot(trzoomfil[0].times(), trzoomfil[0].data, 'k-')
    ax4.axvline(x=picktime - startzoom, color='red', zorder=100)
    ax4.set_title('Zoomed filtered trace')

    ax4.axvline(x=picktime - startzoom + poltime, color='blue', zorder=100)
    ax4.plot(picktime - startzoom + poltime, 0, GetMarker(pol), color='blue')
    plt.text(picktime - startzoom + poltime + 0.05, 0.9*(max(trzoomfil[0].data)), '%s' %(polhint), color='blue')

    plt.tight_layout()

    return fig, pol


def GetTimeWindow(event, station, distance=0):
    depth = event.origins[0].depth / 1000.
    t0 = event.origins[0].time
    v = 15.
    starttime = t0 + depth/v + distance/v 
    endtime = starttime + 60 
    return starttime, endtime


def AskIfPolarityIsGood(traces, station, event, pick, distance, pol):
    msg = '\nThis is station ' + station + '.\n' + \
          'Is this polarity: ' + pol + ' (ENTER (for skipping)| a | x / + / -)?'
    try:
        fig, newpol = ShowPolarity(traces, station, event, pick, distance)
    except IndexError:
        warn('Data unavailable')
        return True
    ax = fig.get_axes()[2]
    ax.text(0.5, 0.9, 'Polarity:%s' %pol, verticalalignment='center', horizontalalignment='center', transform=ax.transAxes)
    fig.show()

    ans = False
    # should be overworked
    while ans == False:
        ans = input(msg)
        if ans in ['a']:
            newpol = pol
        elif ans is '':
            newpol = ''
        elif ans in ['x', '+', '-']:
            newpol = ans
        else:
            ans = False

    plt.close(fig)
    
    return newpol



# Amplitudes
################################

def MakeAmplitudeRequest(traces, station, pick, dist):
    isvalid = True
    p, s = AskForPSAmplitudes(traces, station, pick, dist)
    if p is 'x':
        isvalid = False

    return p, s, isvalid


def AskForPSAmplitudes(traces, station, pick, dist):
    picktime, hint, polPick = pick
    cuttime = 3.
    starttime = picktime - 5
    endtime = starttime + 50+dist/8.

    tr = traces.select(station=station, channel='HH?').detrend('demean').integrate().detrend('linear').slice(starttime=starttime, endtime=endtime)
    df = tr[0].stats.sampling_rate
    pickindex = picktime - starttime-cuttime
    try: 
        sumtrace = np.sqrt(tr[0].data**2 + tr[1].data**2 + tr[2].data**2)
    except:
        print('Could not Sum trace. Maybe one trace is missing.')
        return 'x', 'x'
    sumtrace = sumtrace / max(sumtrace)*100

    fig = plt.figure(figsize=(16, 9))
    ax1 = plt.subplot(211)
    ax1.plot(tr[0].times()[int(cuttime*df):]-cuttime, sumtrace[int(cuttime*df):], color='black')
    ax1.axvline(x=pickindex, color='red', zorder=100)
    ax1.set_ylabel('Normalized amplitude')
    ax1.set_title('SquaredSum trace')

    trfil = tr.filter('highpass', freq=1).slice(starttime=starttime+cuttime, endtime=endtime)
    try: 
        sumtracefil = np.sqrt(trfil[0].data**2 + trfil[1].data**2 + trfil[2].data**2)
    except:
        print('Could not Sum trace. Maybe one trace is missing.')
        return 'x', 'x'
    sumtracefil = sumtracefil / max(sumtracefil)*100

    ax2 = plt.subplot(212)
    ax2.plot(trfil[0].times(), sumtracefil, color='black')    
    ax2.axvline(x=pickindex, color='red', zorder=100)
    ax2.set_ylabel('Normalized amplitude')
    ax2.set_xlabel('Normalized time')
    ax2.set_title('Filtered SquareSum trace')
    cursor = Cursor(ax2, useblit=True, color='gray', linewidth=1)

    plt.tight_layout()
    fig.show()
    
    p,s = AmplitudePicker(sumtracefil, pickindex, df, dist)

    valid = False
    psamp = input('What is the P- and the S-wave amplitude? (Two values or x or a/ENTER for accept)\n').split(' ')
    while not valid:
        if 'a' in psamp or psamp == ['']:
            valid = True
        else:
            try: 
                p = float(psamp[0])
                s = float(psamp[1])
                valid = True
            except (ValueError, IndexError):
                if 'x' in psamp:
                    p = s = 'x'
                    valid = True
                else:
                    psamp = input('No valid input. Give "x" to abort or exactly two values separated with a SPACE\n').split(' ')
    plt.close()
    print(p, s)
    return p, s


def AmplitudePicker(trace, pickindex, df, dist):
    # vp = 6.
    # vs = 6./np.sqrt(3)
    # dt = (dist/vs) - (dist/vp) -0.5 # in seconds
    dt = 1+dist/50.
    p = max(trace[int((pickindex-1)*df):int((pickindex+dt)*df)])
    s = max(trace[int((pickindex+dt)*df):])
    plt.axvline(x=pickindex-1, color='gray', alpha=0.1)
    plt.axvline(x=pickindex+dt, color='gray', alpha=0.1)
    print('The automatic picked P and S maximum:\n', p,s)
    return p, s


def AskIfSPRatioIsGood(traces, station, pick, dist, soverp):
    print('Old ratio is %s' %soverp)
    p,s = AskForPSAmplitudes(traces, station, pick, dist)

    return p,s


# Run HASH 
################################

def GetQuality(ampmisfit, pmisfit, stdr, pcnt):
    if pcnt < 8:
        return 'D'
    stdr = stdr*pcnt
    if ampmisfit <= 50 and pmisfit <= 10 and stdr >= 1000:
        quality = 'A'
    elif ampmisfit <= 70 and pmisfit <= 15 and stdr >= 700:
        quality = 'B'
    elif ampmisfit <= 115 and pmisfit <= 20 and stdr >= 400:
        quality = 'C'
    else:
        quality = 'D'
    return quality


def GetKaganMisfit(angles, pcnt):
    ## has to be reworked
    if pcnt < 8:
        return 1.

    angs = np.array(angles)
    std = np.std(angs)
    avg = np.mean(angs)
    maxang = np.max(angs)

    # # print avg, std, maxang
    # if maxang < 30:
    #     quality = 'A'
    # elif avg < 5:
    #     quality = 'A'
    # elif avg < 15:
    #     quality = 'B'
    # elif avg < 25:
    #     quality = 'C'
    # else:
    #     quality = 'D'
    # return quality
    maxkagan = 120. # what is the max kagan angle?
    misfit = avg/maxkagan
    return misfit


def RunHASH(eventID, workdir):
    """Call an instance of Hash, """
    controlfile = workdir + eventID + '.inp'
    polfile = workdir + eventID + '.pol.hash'
    allresultsfile = workdir + eventID + '.all.fps'
    resultfile = workdir + eventID + '.fps'
    with open('hashpy.inp', 'r') as hashin:
        with open(controlfile, 'w') as ctrlfile:
            for l, line in enumerate(hashin):
                if l==0: ctrlfile.write(polfile + '\n')
                elif l==1: ctrlfile.write(allresultsfile + '\n')
                elif l==2: ctrlfile.write(resultfile + '\n')
                else: ctrlfile.write(line)
    out, err, ret = runBash("./hashchilesp < " + controlfile) #run HASH
    ret = 0
    if ret != 0:
        msg = 'HASH endend with an error:\n' + str(err)
        warn(msg)
    else:
        print(('Done! Results written to ' + resultfile))
        with open(resultfile, 'r') as rf:
            result = rf.read()
        strike, dip, rake, polaritymisfit, stationdistributionratio, amplitudemisfit = result.split()
        nobs = sum(1 for line in open(polfile)) - 1
        # misfit = (float(amplitudemisfit) + float(polaritymisfit) + (1-float(stationdistributionratio)))/100

        ## quality with the kagans angles between preferred and best solutions
        m0 = mtt.magnitude_to_moment(3.)
        mt = mtt.MomentTensor(strike=float(strike), dip=float(dip), rake=float(rake), scalar_moment=m0)

        kagans = []
        with open(allresultsfile) as file:
            for line in file:
                line = line.rsplit()
                str2 = line[0]
                dip2 = line[1]
                rake2 = line[2]
                mt2 = mtt.MomentTensor(strike=float(str2), dip=float(dip2), rake=float(rake2), scalar_moment=m0)
                
                kagan = mtt.kagan_angle(mt, mt2)
                kagans.append(kagan)

        misfit = GetKaganMisfit(kagans, nobs)
        print('Misfit: %s' %(misfit))

    return strike, dip, rake, misfit


def runBash(cmd):
    """Run cmd as a shell comand, return stdout, stderr and return code"""
    p   = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    out = p.stdout.read().strip()
    err = p.stderr.read().strip()
    p.poll()
    ret = p.returncode
    # print out, err, ret
    return out, err, ret


def WriteMechanismToHypfile(hypfile, strike, dip, rake, quality, nobs, workdir):
    basename = hypfile.split('/')[-1].replace('.qml.', '.fps.')
    hypfileout = workdir + '/' + basename
    print(hypfileout)
    with open(hypfile, 'r') as infile:
        with open(hypfileout, 'w') as outfile:
            for inline in infile:
                if inline.startswith('GEOGRAPHIC'):
                    lat = inline.split()[9]
                    lon = inline.split()[11]
                    dep = inline.split()[13]
                if inline.startswith('FOCALMECH'):
                    outline = 'FOCALMECH {} {} {} Mech {} {} {} mf {} nObs {}\n'.format(lat, lon, dep, strike, dip, rake, quality, nobs)
                else:
                    outline = inline
                outfile.write(outline)


def PlotMechanism(eventID, workdir):
    polfile = workdir + eventID + '.pol.hash'
    allresultsfile = workdir + eventID + '.all.fps'
    resultfile = workdir + eventID + '.fps'
    plotfile = workdir + eventID + '.ps'
    _, _, _ = runBash("./plot_mechanism.sh " + allresultsfile + " " + resultfile + " " + polfile + " " + plotfile) #run plotting script
    # _, _, _ = runBash("gv " + plotfile) #run gv
    # print("Made plot. Execute:\nevince " + plotfile + "\nto view.")
