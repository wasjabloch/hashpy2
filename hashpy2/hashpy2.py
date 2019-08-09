#! /usr/bin/env python3

"""
hashpy2 -- an interactive HASH wrapper for meassuring first motion P wave
polarities, P to S amplitude ratios and invert for focal mechanisms
"""

import logging
import yaml
import argparse
from glob import glob
from sys import argv
from subprocess import Popen, PIPE
from warnings import warn
from os.path import isfile
from numpy import cos, pi
import hashpy_utilities as hu



logging.basicConfig(level=logging.INFO)

par = argparse.ArgumentParser(prog='hashpy2',
                              description='An interactive HASH wrapper')
par.add_argument('ID', type=str, help="Origin time of event in format: YYYYMMDDhhmmss")
par.add_argument('--config', type=str, help="Configuration file", default="config.yaml")
#par.add_argument('--review', action="store_true", help="Review existing polarity amplitude data")
#par.add_argument('--add', action="store_true", help="Add more observations to existing data")

print("This is hashpy2")

args = par.parse_args()
ID = args.ID
configfile = args.config
with open(configfile, 'r') as stream:
    params = yaml.safe_load(stream)

resultdir = params['RESULTS'] + '/'
hashfile = resultdir + ID + '.pol.hash'
ctrlfile = resultdir + ID + '.inp'
wvdir = params['WAVEFORMS'] + '/'
stationfile = params['STATIONS']
hypfiles = glob(params['CATALOG'])

timestamp = hu.time_from_ID(ID)
event, hypfile = hu.get_event(timestamp, hypfiles)
traces = hu.make_traces(timestamp, wvdir)

stations = set([str(t.stats.station) for t in traces])
distance = hu.dist2stations(stations, event, stationfile)

logging.info('Found waveforms for %s stations', len(stations))
logging.info('and %s of them have available coordinates\n', len(distance))

review = False
dowaveforms = True
doadd = False
if isfile(hashfile):
    nobs = sum(1 for line in open(hashfile)) - 1
    print('Found existing HASH input file:' + hashfile)
    print('It has ' + str(nobs) + ' observations.')
    print('What should I do?')
    msg = ('(c)ompute focal mechanism\n' +
           '(o)verwrite\n' +
           '(r)eview, or\n' +
           '(a)dd missing stations?\n')
    ans = []
    while ans not in ['c', 'o', 'r', 'a']:
        ans = input(msg)
        if ans == 'n':
            dowaveforms = False
        elif ans == 'a':
            doadd = True
        elif ans == 'r':
            review = True
            dowaveforms = False

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
                pick = hu.PickRequest(traces, station, event, dist)
                ans = hu.AskIfPolarityIsGood(traces, station, event, pick, dist, pol)
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
                    p, s = hu.AskIfSPRatioIsGood(traces, station, pick, dist, soverp)
                    if p is 'x':
                        isvalid = False
                if isvalid:
                    filedata[l] = ('{:5s} {:1s} {:6.4f} {:03.0f} {:2.0f} ' +
                                   '{:03.0f} {:03.0f} {:08.3f}\n').format(
                                       station, pol, float(dist), float(azi),
                                       float(dip), 5, 5, s/p)

    except KeyboardInterrupt:
        with open(hashfile, 'w') as hf:
            hf.writelines(filedata)
        print('KillingTheProcessOnPurpose')

    with open(hashfile, 'w') as hf:
        hf.writelines(filedata)


# Do waveforms #
################
if dowaveforms:
    o = event.origins[0]
    lat = o.latitude
    lon = o.longitude
    slat = o.latitude_errors['uncertainty']*111.2
    slon = o.longitude_errors['uncertainty']*111.2*cos(lat*pi/180)
    sh = (slat+slon)/2
    z = o.depth/1000. # in km
    sz = o.depth_errors['uncertainty']/1000
    t = o.time

    if doadd:
        with open(hashfile, 'r') as hf:
            lines = hf.readlines()
        oldstations = [l.split()[0] for l in lines if l != 0]
    else:
        header = ('{:04d} {:02d} {:02d} {:02d} {:02d} {:05.2f} ' +
                  '{:010.6f} {:010.6f} {:09.6f} {:05.2f} ' +
                  '{:05.2f} {:}\n').format(
                      t.year, t.month, t.day, t.hour, t.minute,
                      t.second + t.microsecond*1e-6,
                      lat, lon, z, sh, sz, ID)
        with open(hashfile, 'w') as hf:
            hf.writelines(header)
        oldstations = []

    nsta = len(distance)
    onsta = len(oldstations)
    logging.info('Skipping %s previously picked stations: %s', onsta, ', '.join(oldstations))
    for ii, station in enumerate(sorted(distance, key=distance.get)):
        dist = distance[station]
        if station not in oldstations:
            logging.info('\nThis is station %s, %s of %s', station, ii, nsta)
            oldstations.append(station)
            pick = hu.PickRequest(traces, station, event, dist)
            pol, isvalid = hu.MakePolarityRequest(traces, station, event, pick, dist)
            if pol in ['+', '-']:
                p, s, isvalid = hu.MakeAmplitudeRequest(traces, station, pick, dist)
            if isvalid:
                outline = '{:5s} {:1s} {:08.3f}\n'.format(station, pol, s/p)
                with open(hashfile, 'a') as hf:
                    hf.writelines(outline)

# Run HASH #
############
#ans = []
#while ans not in ['y', 'n']:
#    ans = input('Want to run HASH now? (y/n)')
#if ans == 'y':
#    nobs = sum(1 for line in open(hashfile)) - 1
#    strike, dip, rake, misfit = RunHASH(sID, resultdir)
#    print(strike, dip, rake, misfit)
#    hu.WriteMechanismToHypfile(hypfile, strike, dip, rake, misfit, nobs, resultdir)
#    hu.PlotMechanism(sID, resultdir)

#    print(KillingTheProcessOnPurpose)

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

