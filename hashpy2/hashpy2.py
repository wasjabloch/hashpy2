#! /usr/bin/env python3

# Script to extract take-off angles from NLL-angle.buf files
import logging
from glob import glob
from sys import argv
from subprocess import Popen, PIPE
from warnings import warn
from os.path import isfile
from hashpy_utilities import ReadNLLAngle, FixBrokenHyp, RunHASH, PlotMechanism, WriteMechanismToHypfile, ShowPolarity, PickRequest
from hashpy_utilities import runBash, MakeTraces, MakeEvent, MakePolarityRequest, AskForPSAmplitudes
from hashpy_utilities import AskForGrids, AskIfPolarityIsGood, AskIfSPRatioIsGood, MakeAmplitudeRequest, dist2stations

logging.basicConfig(level=logging.INFO)

def main(sID):
    review = False
    dowaveforms = True
    doadd = False
    griddirs = ['/home/wasja/data/NLL/model/pamir-trialloc-2d  /home/wasja/data/NLL/time/'] 
    # Directories where NLL angle grids are stored
    resultdir = '/home/wasja/data/hashpy/results/' # saving directory
    hashfile = resultdir + sID + '.pol.hash'
    plotfile = resultdir + sID + '.ps'
    wvdir = '/mnt/data/events/'  # Directory where waveform data is stored
    hypdir = '/home/wasja/FU2017/local/wasja/maps/LOCATIONS/' #+ 'data/nlloc/loc/'
    # Station coordinates as NonLinLoc GTSRCE keyword
    stationfile = '/home/wasja/data/network/chile.nll'

    event, hypfile = MakeEvent(sID, hypdir)
    traces = MakeTraces(sID, wvdir)
    stations = set([str(t.stats.station) for t in traces])
    stationdict = dist2stations(stations, event, stationfile)

    logging.info('Found %s available stations' %(len(stations)))
    logging.info('and %s of them have available coordinates\n' %(len(stationdict)))

    if isfile(hashfile):
        if isfile(plotfile):
            msg = 'Nothing (n), show beachball (b), overwrite (o), review (r), or add missing stations (a)?'
        else:
            msg = 'Nothing (n), overwrite (o), review (r), or add missing stations (a)?'
        ans = []
        nobs = sum(1 for line in open(hashfile)) - 1
        while ans not in ['o', 'n', 'a', 'r']:
            print((hashfile + " exists! It has " + str(nobs) + " observations."))
            print("\nWhat should I do?")
            ans = input(msg)
            if ans == 'n':
                dowaveforms = False
            elif ans == 'a':
                doadd = True
            elif ans == 'r':
                review = True
                dowaveforms = False
            elif ans == 'b' and isfile(plotfile):
                runBash('gv ' + plotfile + '&')

    # Review #
    ##########
    if review:
        outbuffer = []
        with open(hashfile, 'r') as hf:
            filedata = hf.readlines()
        try:
            for l, line in enumerate(filedata):
                if l == 0:
                    pass
                else:
                    station, pol, dist, azi, dip, aerr, derr, soverp = line.split()
                    dist = float(dist)
                    pick = PickRequest(traces, station, event, dist)
                    ans = AskIfPolarityIsGood(traces, station, event, pick, dist, pol)
                    if ans in ['+', '-']:
                        pol = ans
                        isvalid = True
                    elif ans in ['x']:
                        filedata[l] = ''
                        isvalid = False
                    elif ans in ['a']:
                        isvalid = True
                        pol = line.rsplit()[1]
                    else:
                        isvalid = False

                    if isvalid:
                        p, s = AskIfSPRatioIsGood(traces, station, pick, dist, soverp)
                        if p is 'x':
                            isvalid = False
                    if isvalid:
                        filedata[l] = '{:5s} {:1s} {:6.4f} {:03.0f} {:2.0f} {:03.0f} {:03.0f} {:08.3f}\n'.format(station, pol, float(dist), float(azi), float(dip), 5, 5, s/p)
                        
        except:
            with open(hashfile, 'w') as hf:
                hf.writelines(filedata)
            print(KillingTheProcessOnPurpose)

        with open(hashfile, 'w') as hf:
            hf.writelines(filedata)


    # Do waveforms #
    ################
    if dowaveforms:
        lon = event.origins[0].longitude
        lat = event.origins[0].latitude
        z = event.origins[0].depth/1000. # in km
  
        if doadd:
            with open(hashfile, 'r') as hf:
                lines = hf.readlines()
            oldstations = [l.split()[0] for l in lines if l != 0]
        else:
            header = '{:} {:05.2f} {:010.6f} {:010.6f} {:09.6f} {:05.2f} {:05.2f} {:}' .format(sID, 0, lon+360, lat+360, z, 0, 0, hypfile + '\n')
            with open(hashfile, 'w') as hf:
                hf.writelines(header)
            oldstations = []

        nsta = len(stationdict)
        onsta = len(oldstations)
        logging.info('Found {:} stations: {:}'.format(nsta, ' '.join(oldstations)))
        ii = 0 # running variable, for progress
        for ii, station in enumerate(sorted(stationdict, key=stationdict.get, reverse=False)):
            
            dist = stationdict[station]
            ii += 1
            logging.info('Number {:}, station: {:}'.format(ii, station))
            for griddir in griddirs:
                if glob('%s*.P.%s*.angle.hdr' %(griddir,station)) != []:
                    grid = glob('%s*.P.%s*.angle.hdr' %(griddir,station))[0]
                    prefix, phase, station, _, _ = grid.split('/')[-1].split('.')
                    logging.info('Number {:}, station: {:}'.format(ii, station))
                    if station not in oldstations:
                        oldstations.append(station)
                        logging.info('\nThis is station {:}, {:} of {:}'.format(
                            station, ii, nsta))

                        pick = PickRequest(traces, station, event, dist)
                        pol, isvalid = MakePolarityRequest(traces, station, event, pick, dist)
                        if pol in ['+', '-']:
                            p, s, isvalid = MakeAmplitudeRequest(traces, station, pick, dist)

                        if isvalid:
                            azi, dip, qual = ReadNLLAngle(station, lon, lat, z, grid)
                            outline = '{:5s} {:1s} {:6.4f} {:03.0f} {:2.0f} {:03.0f} {:03.0f} {:08.3f}\n'.format(station, pol, dist, azi, 180-dip, 5, 5, s/p)
                            with open(hashfile, 'a') as hf:
                                hf.writelines(outline)

    # Run HASH #
    ############
    ans = []
    while ans not in ['y', 'n']:
        ans = input('Want to run HASH now? (y/n)')
    if ans == 'y':
        nobs = sum(1 for line in open(hashfile)) - 1
        strike, dip, rake, misfit = RunHASH(sID, resultdir)
        print(strike, dip, rake, misfit)
        WriteMechanismToHypfile(hypfile, strike, dip, rake, misfit, nobs, resultdir)
        PlotMechanism(sID, resultdir)
#    print(KillingTheProcessOnPurpose)

#This idiom means the below code only runs when executed from the command line
if __name__ == '__main__':
    print("This is HASHPY")
    if len(argv) != 2:
        raise SyntaxError('Parse Short ID of the event as only argument')
    sID = argv[1]

    main(sID)
#
#    # the terminal will be open/activated the whole time
#    # has to be killed manually with the except clause!! 
#    try:
#        cmd = 'ID=$(xdotool getactivewindow); while true; do sleep 0.2; xdotool windowfocus $ID; done;'
#        p   = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
#        main(sID)   
#    except Exception as e: 
#        print(e)
#        cmd = 'pkill -f xdotool'
#        p   = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
#        exit()

