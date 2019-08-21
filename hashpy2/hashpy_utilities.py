#! /usr/bin/python3
from subprocess import Popen, PIPE
from struct import unpack
from glob import glob
import logging
from warnings import warn, filterwarnings
from copy import deepcopy
import numpy as np
from obspy import read_events, read, UTCDateTime
from obspy.signal.util import util_geo_km
from obspy.signal.trigger import recursive_sta_lta, pk_baer
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Cursor


def print_config_tipps():
    """ Print meaning of configuration parameters to screen """
    print('')
    print('----------------------------------')
    print('Control parameters of hashpy2 are:')
    print('')
    print('CATALOG:')
    print('    Regular expression pattern for NonLinLoc .hyp files holding')
    print('    hypocenter, location error and phase information.')
    print('')
    print('RESULTS:')
    print('    Directory in which the results will be stored. Must be present.')
    print('')
    print('WAVEFORMS:')
    print('    Directory in which the waveform files are stored. Directory')
    print('    structure is:')
    print('    WAVEFORMS/YYYY/YYMMDDhhmmss, where')
    print('    - YYYY is four digit year,')
    print('    - YY, MM, DD, hh, mm, ss are two digit year, month, day,')
    print('      hour, minute, and second.')
    print('    Make sure that all three components are present and that they')
    print('    are present only once.')
    print('')
    print('VELOCITIES:')
    print('    List of velocity model files holding two columns:')
    print('    - depth (km), velocity (km/s)')
    print('    Velocities are interpolated linearly. Use zero thickness layers')
    print('    to model discontinuities.')
    print('')
    print('LOGGING:')
    print('    Choose verbosity of program output.')
    print('')
    print('-------------------------------')
    print('Control parameters of HASH are:')
    print('')
    print('dang:')
    print('    Angle increment for grid search.')
    print('')
    print('nmc:')
    print('    Number of perutbations of take-off angles for different source')
    print('    depths and velocity models.')
    print('')
    print('maxout:')
    print('    Maximum number focal mechanisms that match misfit critria to')
    print('    return.')
    print('')
    print('cangle:')
    print('    Angular distance between different families of focal mechanisms.')
    print('')
    print('prob_max:')
    print('    Fraction of focal mechanisms that need to be within cangle')
    print('    to make up a new famlily of focal mechanisms.')
    print('')
    print('qbadfac:')
    print('    log10 of uncertainty factor for s/p ratios.')
    print('')


def write_default_config(filename):
    """
    Write default configuration parameters to file.
    filename: (str) Name of the file
    """
    with open(filename, 'w') as f:
        f.write('CATALOG: ./hyp/*.hyp\n')
        f.write('RESULTS: ./results\n')
        f.write('WAVEFORMS: ./events\n')
        f.write('STATIONS: stations.nll\n')
        f.write('VELOCITIES:\n')
        f.write('    - vmodel.zv\n')
        f.write('LOGGING: WARNING\n')
        f.write('dang: 1\n')
        f.write('nmc: 50\n')
        f.write('maxout: 100\n')
        f.write('cangle: 45\n')
        f.write('prob_max: 0.2\n')
        f.write('qbadfac: 0.2\n')


def get_event(t0, hypfiles):
    """
    Retrun obspy.Event near time t0 from list of hypfiles
    """
    hypfile = find_hypfile_from_t0(t0, hypfiles)
    event = read_events(hypfile)[0]
    return event, hypfile


def get_t0_from_hyp(hypfile):
    """ Get origin time from hypfile """
    with open(hypfile, 'r') as hf:
        for line in hf:
            if line.startswith('GEOGRAPHIC'):
                f = line.split()
                t0 = UTCDateTime(int(f[2]), int(f[3]), int(f[4]), int(f[5]),
                                 int(f[6])) + float(f[7])
                return t0


def find_hypfile_from_t0(t0, hypfiles):
    """
    Search for hypfile that holds event with origin time closest to t0
    hypfiles: (list) list of NonLinLoc hypfiles
    t0: (UTCDateTime) origin time
    """
    t0 = float(t0)
    hfs = {float(get_t0_from_hyp(hypfile)): hypfile for hypfile in hypfiles}
    dtmin = np.inf
    for t in sorted(hfs):
        dt = abs(t - t0)
        if dt < dtmin:
            dtmin = dt
            thishyp = hfs[t]
        else:
            break
    return thishyp


def find_ID_from_t0(t0, sIDs):
    """
    Search for the ID represents the event with origin time closest to t0
    sIDs: (list) list of short IDs of format YYMMDDhhmmss
    t0: (UTCDateTime) origin time
    """
    t0 = float(t0)
    ids = {float(UTCDateTime.strptime(thisID, '%y%m%d%H%M%S')): thisID
           for thisID in sIDs}
    dtmin = np.inf
    for t in sorted(ids):
        dt = abs(t - t0)
        if dt < dtmin:
            dtmin = dt
            ID = ids[t]
        else:
            break
    return ID


def time_from_ID(ID):
    """Decode timestring YYYYMMDDhhmmss to UTCDateTime"""
    return UTCDateTime.strptime(ID, '%Y%m%d%H%M%S')


def make_traces(timestamp, wvdir):
    """
    Read seismic traces stored in directory structre YYYY/YYMMDDhhmmss
    """
    filterwarnings("ignore", category=UserWarning,
                   message='Record contains a fractional second')
    year = timestamp.strftime('%Y')
    IDs = [d.split('/')[-1] for d in glob(wvdir + year + '/*')]
    ID = find_ID_from_t0(timestamp, IDs)
    wavedir = wvdir + year + '/' + ID
    traces = read(wavedir + '/*[BH]H*')
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
    r = 6371.  # km
    distance = r*2.*np.arcsin(np.sqrt(a))

    return distance


def dist2stations(stations, event, filename):
    """
    Compute dictionary that holds event - station distances
    station -> distance (km)
    """
    stationdict = {}
    with open(filename, 'r') as stationfile:
        for line in stationfile:
            if not line.startswith('GTSRCE'):
                continue
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

# Picker
################################


def PickRequest(traces, station, event, dist):
    for pick in event.picks:
        pickvalid = False
        if ((pick.waveform_id['station_code'] == station) and
           (pick.phase_hint == 'P')):
            pickvalid = True
            picktime = pick.time
            hint = 'b'
            pol = 'x'
            break
    if not pickvalid:
        picktime, pol = AutoPicker(traces, station, event, dist)
        hint = 'a'
    pick = [picktime, hint, pol]

    return pick


def AutoPicker(traces, station, event, dist, component='?HZ'):
    starttime, endtime = GetTimeWindow(event, station)
    tr = traces.select(station=station, channel=component).slice(
        starttime=starttime, endtime=endtime).detrend('demean')
    if len(tr) == 0:
        return 0, ''
    df = tr[0].stats.sampling_rate
    p_pick, phase_info = pk_baer(tr[0].data, df, 20, 60, 7.0, 12.0, 100, 100)
    picktime = starttime + (p_pick / df)
#    print('Phaseinfo', picktime, phase_info)
    if phase_info == [] or phase_info == '':
        pol = 'x'
    else:
        pol = phase_info[2]

    return picktime, pol


def PolarityPicker(picktime, trace, gs):
    poltrace = trace.copy()

    pol, poltime = PolarityStd(picktime, poltrace.filter('highpass', freq=1.0),
                               gs[3, 1])
    polhint = 'Std'
    poltrace = trace.copy()
    if pol in ['+', '-']:
        PolarityGradient(picktime, poltrace, gs[3, 0])
        # PolarityFB(picktime, poltrace, gs, mode='gradient')
    else:
        pol, poltime = PolarityGradient(picktime, poltrace, gs[3, 0])
        polhint = 'Gradient'
        if pol in ['+', '-']:
            # PolarityFB(picktime, poltrace, gs, mode='gradient')
            pass
        else:
            # pol, poltime = PolarityFB(picktime, poltrace, gs, mode='gradient')
            # polhint = 'FB'
            pass

    return pol, poltime, polhint


def PolarityGradient(picktime, trace, gs, color='black'):
    timewindow = 0.5
    poltrace = trace.copy()

    tr = poltrace.slice(starttime=picktime, endtime=picktime+timewindow)
    try:
        starttime = tr.stats.starttime
    except AttributeError:
        return '', 0

    tmpstd = np.std(tr.data)
    if np.isnan(tmpstd):
        pol = 'x'
        poltime = 0
        return pol, poltime

    # simple gradient picker
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
                pol = 'x'
                poltime = 0

    # advanced gradient picker, measuring the gradient change
    # and check if the gradient before the polarity pick has the same sign for
    # some samples
    elif True:
        for ii in range(len(tr.data)-2):
            d1 = tr.data[ii]
            d2 = tr.data[ii+1]
            d3 = tr.data[ii+2]

            if abs(d2-d1) < abs(d3-d2):
                for nn in range(len(tr.data)-3-ii):
                    nn += 1
                    d1 = tr.data[ii+nn]
                    d2 = tr.data[ii+1+nn]
                    d3 = tr.data[ii+2+nn]
                    pol = 'x'
                    poltime = 0
                    if np.sign(d2-d1) != np.sign(d3-d2):
                        poltime = tr.times()[ii+1+nn]
                        if np.sign(d2-d1) == 1.0:
                            pol = '+'
                        elif np.sign(d2-d1) == -1.0:
                            pol = '-'
                        break

                if pol in ['+', '-']:
                    try:
                        dpp2 = tr.data[ii+3+nn]
                        dpp1 = tr.data[ii+2+nn]
                        dp = tr.data[ii+1+nn]
                        dpm1 = tr.data[ii+nn]
                        dpm2 = tr.data[ii-1+nn]
                    except IndexError:
                        continue
                    if (np.sign(dp-dpm1) != np.sign(dpm1-dpm2) and
                       np.sign(dpp2-dpp1) != np.sign(dpp1-dp)):
                        continue
                    else:
                        break
            else:
                pol = 'x'
                poltime = 0
        else:
            pol = 'x'
            poltime = 0

    # gradient picker, the gradient has to be stable for few more samples
    else:
        # now it will look every sample after the intial polarity pick and
        # check if the sign is reversed another idea would be to check if
        # the amplitude is not changing that much for the next samples...

        checkLength = 2
        for ii in range(len(tr.data)-2-checkLength):
            d1 = tr.data[ii]
            d2 = tr.data[ii+1]
            d3 = tr.data[ii+2]

            if np.sign(d2-d1) != np.sign(d3-d2):
                nn = 0
                while np.sign(d2-d1) != (np.sign(tr.data[ii+2+nn] -
                                         tr.data[ii+1+nn])):
                    nn += 1
                    if nn > checkLength:
                        poltime = tr.times()[ii+1]
                        if np.sign(d2-d1) == 1.0:
                            pol = '+'
                        elif np.sign(d2-d1) == -1.0:
                            pol = '-'

                        break
                    else:
                        pol = 'x'
                        poltime = 0
                else:
                    pol = 'x'
                    poltime = 0

                if pol != 'x':
                    break
            else:
                pol = 'x'
                poltime = 0

    ax0 = plt.subplot(gs)
    ax0.set_title('Polarity Picker Gradient')
    ax0.plot(tr.times(), tr.data, color=color)
    ax0.axvline(x=picktime - starttime + poltime, color=color, zorder=100,
                label='Gr')
    ax0.plot(picktime - starttime + poltime, np.median(tr.data),
             GetMarker(pol), color=color, zorder=100)

    return pol, poltime


def PolarityStd(picktime, trace, gs, color='black'):
    ax0 = plt.subplot(gs)
    ax0.set_title('STD Picker')

    poltrace = trace.copy()
    noisestart = 1.5
    noiseend = 0.5
    stdtr = poltrace.slice(starttime=picktime-noisestart,
                           endtime=picktime-noiseend)
    std = np.std(stdtr.data)

    poltrace = trace.copy()
    signalstart = 0.5
    signalend = 1.5
    stdSigtr = poltrace.slice(starttime=picktime+signalstart,
                              endtime=picktime+signalend)
    stdSig = np.std(stdSigtr.data)

    if np.isnan(std) or np.isnan(stdSig):
        pol = 'x'
        poltime = 0
        return pol, poltime

    ratioStd = stdSig/std

    timewindow = 1.
    timebeforepick = 0.1
    poltrace = trace.copy()
    tr1 = poltrace.slice(starttime=picktime-timebeforepick,
                         endtime=picktime+timewindow)

    starttime = tr1.stats.starttime

    if len(tr1.data) == 0:
        return '', 0

    ax0.plot(tr1.times(), tr1.data, color='black')
    ax0.axvline(x=picktime - starttime, color='red', zorder=100, label='Pick')
#    print(picktime, starttime, picktime- starttime)

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
#    print(minCnt, ratioStd, stdsinglef)

    # minCounts
    stdfacs = np.arange(stdmin, stdmax, 1)

    norm = mpl.colors.Normalize(vmin=stdmin, vmax=stdmax)
    cmap = mpl.cm.get_cmap('Spectral')

    cnt = 0

    fpol = ''
    for stdfac in stdfacs:
        color = cmap(norm(stdfac))
        stdpos = tr1.data[0] + stdfac*std
        stdneg = tr1.data[0] - stdfac*std
        for ii in range(len(tr1.data)-1):
            ii += 1
            if tr1.data[ii] > stdpos or tr1.data[ii] < stdneg:
                poltime = tr1.times()[ii]
                if np.sign(tr1.data[ii]-tr1.data[ii-1]) == 1.0:
                    pol = '+'
                elif np.sign(tr1.data[ii]-tr1.data[ii-1]) == -1.0:
                    pol = '-'
                break

            else:
                pol = 'x'
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
        pol = 'x'
        poltime = 0.

    # stdsinglef
    stdpos = tr1.data[0] + stdsinglef*std
    stdneg = tr1.data[0] - stdsinglef*std
    for ii in range(len(tr1.data)-1):
        ii += 1
        if tr1.data[ii] > stdpos or tr1.data[ii] < stdneg:
            poltime = tr1.times()[ii]
            if np.sign(tr1.data[ii]-tr1.data[ii-1]) == 1.0:
                pol = '+'
            elif np.sign(tr1.data[ii]-tr1.data[ii-1]) == -1.0:
                pol = '-'
            break

        else:
            pol = 'x'
            poltime = 0

    ax0.axvline(x=poltime, color='green', zorder=100, label='STD')
    ax0.plot(poltime, 0, GetMarker(pol), color='green', zorder=100)
    ax0.axhline(y=stdpos, color='green', zorder=-10, alpha=0.5)
    ax0.axhline(y=stdneg, color='green', zorder=-10, alpha=0.5)

    poltime = poltime-(picktime-starttime)

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

def AskForPolarity(traces, station, event, pick, distance):
    """ Ask the user for polarity observation """
    fig, pol = ShowPolarity(traces, station, event, pick, distance)
    fig.show()

    quality = {'++': 1, '--': 1, '': 1, '+': 0, '-': 0, 'x': 0, 'L': 0}
    polarity = {'++': '+', '+': '+', '--': '-', '-': '-', 'x': 'x', 'L': 'x'}
    ans = None
    leave = False
    while ans not in ['++', '+', 'x', '-', '--', '', 'L']:
        print('(ENTER) to accept polarity: {:}{:}'.format(pol, pol))
        print('(++, +, x, -, --) to correct')
        print('(L)eave picking routine')
        ans = input('')
    if ans != '':
        pol = polarity[ans]
    if ans == 'L':
        leave = True
    qpol = quality[ans]
    plt.close(fig)

    return pol, qpol, leave


def ShowPolarity(traces, station, event, pick, distance, component='?HZ'):
    picktime, hint, polPick = pick
    starttime, endtime = GetTimeWindow(event, station, distance)

    rawtr = deepcopy(traces.select(station=station, channel=component).slice(
        starttime=starttime, endtime=endtime).detrend('demean').integrate())
    tr = rawtr.copy()
    dt = rawtr[0].stats.starttime - rawtr[0].stats.endtime

    fig = plt.figure(figsize=(8, 10))
    wm = plt.get_current_fig_manager()
    wm.window.attributes('-topmost', 0)
    gs = gridspec.GridSpec(4, 2)

    # Polarity window
    pol, poltime, polhint = PolarityPicker(picktime, tr[0], gs)
    tr = rawtr.copy()

    # Raw data
    ax1 = plt.subplot(gs[0, :])
    ax1.plot(tr[0].times(), tr[0].data, 'k-')
    ax1.axvline(x=picktime - starttime, color='red', zorder=100)
    ax1.set_title('Raw trace {:}.{:}\nDistance: {:.0f} km'.format(
                  station, component, distance))

    # Zoomed data raw
    dt = 0
    while True:
        if dt > 1:
            logging.warning('Could not make trace closeup.')
            fig = plt.figure()
            return fig, 'x'
        startzoom = picktime - 1. + dt
        endzoom = picktime + 1.5 - dt
        tr = rawtr.copy()
        trzoom = tr.slice(starttime=startzoom, endtime=endzoom).detrend(
                          'linear')
        try:
            trzoom[0]
        except IndexError:
            dt += 0.1  # Zoom closer to avoid empty trace
        else:
            break

    ax3 = plt.subplot(gs[2, 0])
    ax3.set_title('Zoomed raw trace')
    ax3.plot(trzoom[0].times(), trzoom[0].data, 'k-')
    ax3.axvline(x=picktime - startzoom, color='red', zorder=100)

    ax3.axvline(x=picktime - startzoom + poltime, color='blue', zorder=100)
    ax3.plot(picktime - startzoom + poltime, 0, GetMarker(pol), color='blue')
    plt.text(picktime - startzoom + poltime + 0.05, 0.9*(max(trzoom[0].data)),
             '%s' % (polhint), color='blue')

    # Filtered data
    freq = 1.
    tr = rawtr.copy()
    trfil = tr.filter('highpass', freq=freq)
    ax2 = plt.subplot(gs[1, :], sharex=ax1)
    ax2.plot(trfil[0].times(), trfil[0].data, 'k-')
    ax2.axvline(x=picktime - starttime, color='red', zorder=100)
    ax2.set_title('Filtered trace, HP %s Hz' % freq)

    # Zoomed data filtered
    trzoomfil = trfil.slice(starttime=startzoom, endtime=endzoom)
    ax4 = plt.subplot(gs[2, 1])
    ax4.plot(trzoomfil[0].times(), trzoomfil[0].data, 'k-')
    ax4.axvline(x=picktime - startzoom, color='red', zorder=100)
    ax4.set_title('Zoomed filtered trace')

    ax4.axvline(x=picktime - startzoom + poltime, color='blue', zorder=100)
    ax4.plot(picktime - startzoom + poltime, 0, GetMarker(pol), color='blue')
    plt.text(picktime - startzoom + poltime + 0.05,
             0.9*(max(trzoomfil[0].data)),
             '%s' % (polhint), color='blue')

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
          'Is this polarity: ' + pol + \
          ' (ENTER (for skipping)| a | x / + / -)?'
    try:
        fig, newpol = ShowPolarity(traces, station, event, pick, distance)
    except IndexError:
        warn('Data unavailable')
        return True
    ax = fig.get_axes()[2]
    ax.text(0.5, 0.9, 'Polarity:%s' % pol, verticalalignment='center',
            horizontalalignment='center', transform=ax.transAxes)
    fig.show()

    ans = False
    # should be overworked
    while not ans:
        ans = input(msg)
        if ans in ['a']:
            newpol = pol
        elif ans is 'x':
            newpol = 'x'
        elif ans in ['x', '+', '-']:
            newpol = ans
        else:
            ans = False

    plt.close(fig)

    return newpol


# Amplitudes
################################

def get_spratio(traces, station, pick, dist):
    """
    Consumes seismic traces, station, pick and distance information
    Automatically dertmines the amplitude of P and S amplitudes
    Asks the user for confirmation or input
    """
    picktime, _, _ = pick
    cuttime = 3.
    starttime = picktime - 5
    endtime = starttime + 50+dist/8.

    tr = traces.select(station=station, channel='[BH]H?').detrend(
        'demean').integrate().detrend('linear').slice(
            starttime=starttime, endtime=endtime)
    df = tr[0].stats.sampling_rate
    pickindex = picktime - starttime-cuttime
    try:
        sumtrace = np.sqrt(tr[0].data**2 + tr[1].data**2 + tr[2].data**2)
    except (IndexError, ValueError):
        logging.warning('Could not Sum trace. Maybe one trace is missing or ' +
                        'corrupt.')
        return 0, 0
    sumtrace = sumtrace / max(sumtrace)*100

    fig = plt.figure(figsize=(8, 10))
    ax1 = plt.subplot(311)
    ax1.plot(tr[0].times()[int(cuttime*df):]-cuttime,
             sumtrace[int(cuttime*df):], color='black')
    ax1.set_ylabel('Normalized amplitude')
    ax1.set_title('Cartesian sum of traces\nDistance: {:.0f}km'.format(dist))

    trfil = tr.filter('highpass', freq=1).slice(starttime=starttime+cuttime,
                                                endtime=endtime)
    try:
        sumtracefil = np.sqrt(trfil[0].data**2 + trfil[1].data**2 +
                              trfil[2].data**2)
    except (IndexError, ValueError):
        logging.warning('Could not Sum trace. Maybe one trace is missing or ' +
                        'corrupt.')
        return 0, 0
    sumtracefil = sumtracefil / max(sumtracefil)*100

    ax2 = plt.subplot(312, sharex=ax1)
    ax2.plot(trfil[0].times(), sumtracefil, color='black')
    ax2.set_ylabel('Normalized amplitude')
    ax2.set_xlabel('Normalized time')
    ax2.set_title('Filtered Cartesian sum of traces')
    Cursor(ax2, useblit=True, color='gray', linewidth=1)

    ax3 = plt.subplot(313, sharex=ax2)
    ax3.plot(trfil[0].times(), trfil[0].data, linewidth=.5, color='black')
    ax3.plot(trfil[0].times(), trfil[1].data, linewidth=.5, color='steelblue')
    ax3.plot(trfil[0].times(), trfil[2].data, linewidth=.5, color='purple')
    ax3.plot(pickindex, 0, marker='|')

    p, s = AmplitudePicker(sumtracefil, pickindex, df, dist)

    for a in [ax1, ax2]:
        a.axhline(y=p, zorder=-100, color='maroon')
        a.axhline(y=s, zorder=-100, color='blue')
        a.text(0, p, '{:3.0f}'.format(p), va='bottom', ha='left',
               color='maroon')
        a.text(0, s, '{:3.0f}'.format(s), va='bottom', ha='left', color='blue')
        a.plot(pickindex, 0, marker='|')

    plt.tight_layout()
    fig.show()

    p, s = ask_for_amplitudes(p, s)

    if p == 'x' or s == 'x':
        plt.close()
        return 0

    sp = s/p

    plt.close()
    return sp


def ask_for_amplitudes(p='x', s='x'):
    """
    Print the amplidue question to console
    p: (float) a priori p wave amplitude)
    s: (float) a priori s wave amplitude)
    """
    while True:
        print('What are the P and S wave amplitudes?')
        print('(ENTER) accept {:.0f} {:.0f}'.format(p, s))
        print('(x) for undefined')
        print('(two values) best estimate and error estimate')
        ans = input('').split(' ')
        if ans == ['']:
            break
        elif 'x' in ans:
            p = 'x'
            s = 'x'
            break
        else:
            try:
                p = float(ans[0])
                s = float(ans[1])
            except (ValueError, IndexError):
                print('Invalid input: ' + ans)
                continue
            if p <= 0 or s <= 0:
                print('Amplitudes must be greater 0')
            else:
                break
    return p, s


def AmplitudePicker(trace, pickindex, df, dist):
    dt = (dist/3.5) - (dist/6)  # in seconds
    p_wind = 2  # seconds, P window
    s_wind = 10  # seconds, s window
    p = max(trace[int((pickindex-1)*df):int((pickindex+p_wind)*df)])
    s = max(trace[int((pickindex+dt)*df):int((pickindex+dt+s_wind)*df)])
    return p, s


def AskIfSPRatioIsGood(traces, station, pick, dist, soverp):
    print('Old ratio is %s' % soverp)
    sp = get_spratio(traces, station, pick, dist)

    return sp


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


def create_HASH_runfile(hashfile, polfile, params, ID):
    """
    Make inputfile for HASH
    hashfile: (str) Name of HASH input file
    polfile: (str) Name of polarity file
    params: parameter dictionary
    ID: (str) Event identification number
    """
    goodpol = 0
    with open(polfile, 'r') as pf:
        for n, line in enumerate(pf):
            if n > 0:
                line = line.split()
                try:
                    goodpol += int(line[2])
                except IndexError:
                    continue

    badpol = n - goodpol

    #  Ask to overwrite bad polarity readings
    while True:
        print('(ENTER) confirm {:d} bad polarity readings, or type'.format(
              badpol))
        print('(integer) number of bad polartiy readings.')
        ans = input('')
        if ans == '':
            break
        try:
            badpol = int(ans)
            break
        except ValueError:
            continue

    # Ask to override maximum distance
    delmax = 20040  # half earth circumference
    while True:
        print(('(ENTER) confirm {:.0f} km maximum source receiver distance, ' +
               'or type').format(delmax))
        print('(float) km maximum source receiver distance.')
        ans = input('')
        if ans == '':
            break
        try:
            delmax = float(ans)
            break
        except ValueError:
            continue


    with open(hashfile, 'w') as hf:
        hf.write('{:}/{:}.pol.hash\n'.format(params['RESULTS'], ID))
        hf.write(params['STATIONS'] + '\n')
        hf.write('{:}/{:}.fps\n'.format(params['RESULTS'], ID))
        hf.write('{:}/{:}.all.fps\n'.format(params['RESULTS'], ID))
        hf.write('{:}/{:}.rays\n'.format(params['RESULTS'], ID))
        hf.write('{:}\n'.format(params['dang']))
        hf.write('{:}\n'.format(params['nmc']))
        hf.write('{:}\n'.format(params['maxout']))
        hf.write('{:}\n'.format(delmax))
        hf.write('{:}\n'.format(badpol))
        hf.write('{:}\n'.format(params['qbadfac']))
        hf.write('{:}\n'.format(params['cangle']))
        hf.write('{:}\n'.format(params['prob_max']))
        hf.write('{:}\n'.format(len(params['VELOCITIES'])))
        for vz in params['VELOCITIES']:
            hf.write('{:}\n'.format(vz))


def RunHASH(controlfile):
    """Call an instance of HASH"""
    with open(controlfile, 'r') as cf:
        for n, line in enumerate(cf):
            if n == 2:
                resultfile = line.strip()
    out, err, ret = runBash("hash_hashpy1D < " + controlfile)
    ret = 0
    if ret != 0:
        msg = 'HASH endend with an error:\n' + str(err)
        logging.warning(msg)
    else:
        logging.info('HASH endded successfully')
    logging.info('HASH output:\n' + out.decode('utf8'))
    ang_rms = np.inf
    qual = 'Z'  # worst quality ever
    with open(resultfile, 'r') as rf:
        for nn, line in enumerate(rf):
            if line.split()[4] <= qual and float(line.split()[3]) <= ang_rms:
                strike, dip, rake, ang_rms = map(int, line.split()[0:4])
                qual = line.split()[4]
    if nn > 0:
        logging.warning('Found multiple solutions. Reporting best one.')

    return strike, dip, rake, ang_rms


def runBash(cmd):
    """Run cmd as a shell comand, return stdout, stderr and return code"""
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    out = p.stdout.read().strip()
    err = p.stderr.read().strip()
    p.poll()
    ret = p.returncode
    return out, err, ret


def WriteMechanismToHypfile(hypfile, strike, dip, rake, quality, nobs,
                            workdir):
    basename = hypfile.split('/')[-1]
    hypfileout = workdir + '/' + basename
    with open(hypfile, 'r') as infile:
        with open(hypfileout, 'w') as outfile:
            for inline in infile:
                if inline.startswith('GEOGRAPHIC'):
                    lat = inline.split()[9]
                    lon = inline.split()[11]
                    dep = inline.split()[13]
                if inline.startswith('FOCALMECH'):
                    outline = ('FOCALMECH {} {} {} Mech ' +
                               '{} {} {} mf {} nObs {}\n').format(
                               lat, lon, dep, strike, dip, rake, quality, nobs)
                else:
                    outline = inline
                outfile.write(outline)
    return hypfileout


def PlotMechanism(resultdir, hypfile, ID):
    allfpsf = resultdir + ID + '.all.fps'
    angf = resultdir + ID + '.rays'
    plotfile = resultdir + ID + '.ps'
    cmd = 'plot_mechanism.sh {:} {:} {:} {:}'.format(
            hypfile, allfpsf, angf, plotfile)
    out, _, _ = runBash(cmd)  # run plotting script
    logging.info('The plotting script returned:\n' + out.decode('utf-8'))
    _, _, _ = runBash("gv " + plotfile)  # run gv
