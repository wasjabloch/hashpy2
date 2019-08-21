#! /usr/bin/env python3

"""
hashpy2 -- an interactive HASH wrapper for meassuring first motion P wave
polarities, P to S amplitude ratios and invert for focal mechanisms
"""

import logging
import yaml
import argparse
import sys
from glob import glob
from os.path import isfile
from numpy import cos, pi
import hashpy_utilities as hu


logging.basicConfig(level=logging.INFO)

par = argparse.ArgumentParser(prog='hashpy2',
                              description='An interactive HASH wrapper',
                              epilog=('Authors: Wasja Bloch ' +
                                      '(wasja@gfz-potsdam.de) ' +
                                      'and Lukas Lehmann ' +
                                      '(luklehma@uni-potsdam.de)'))

par.add_argument('ID', type=str,
                 nargs='?',
                 help="Origin time of event in format: YYYYMMDDhhmmss",
                 default=None)
par.add_argument('--config', type=str,
                 help="Configuration file [config.yaml]",
                 default="config.yaml")
par.add_argument('--setup', action="store_true",
                 help=('Print configuration tipps. ' +
                       'Produce default_config.yaml.'))

print("This is hashpy2")

args = par.parse_args()
ID = args.ID

if args.setup:
    defaultf = 'default_config.yaml'
    hu.print_config_tipps()
    hu.write_default_config(defaultf)
    sys.exit('Default configuration file written to: ' + defaultf)

if not ID:
    par.print_help()
    sys.exit()

#  Read config file
configfile = args.config
with open(configfile, 'r') as stream:
    params = yaml.safe_load(stream)
logging.basicConfig(level=params['LOGGING'])


#  Set directories
resultdir = params['RESULTS'] + '/'
hashfile = resultdir + ID + '.pol.hash'
ctrlfile = resultdir + ID + '.inp'
wvdir = params['WAVEFORMS'] + '/'
stationfile = params['STATIONS']
hypfiles = glob(params['CATALOG'])


#  Get event information and seismic traces
timestamp = hu.time_from_ID(ID)
event, hypfile = hu.get_event(timestamp, hypfiles)
traces = hu.make_traces(timestamp, wvdir)
stations = set([str(t.stats.station) for t in traces])
distance = hu.dist2stations(stations, event, stationfile)

logging.info('Found waveforms for %s stations', len(stations))
logging.info('and %s of them have available coordinates\n', len(distance))

#  Ask what to do if .pol.hash is present
dowaveforms = True
doadd = False
if isfile(hashfile):
    nobs = sum(1 for line in open(hashfile)) - 1
    print('Found existing HASH input file:' + hashfile)
    print('It has ' + str(nobs) + ' observations.')
    print('What should I do?')
    msg = ('(c)ompute focal mechanism\n' +
           '(a)dd missing stations?\n')
    ans = []
    while ans not in ['c', 'a']:
        ans = input(msg)
        if ans == 'c':
            dowaveforms = False
        elif ans == 'a':
            doadd = True

# Do waveforms #
################
if dowaveforms:
    o = event.origins[0]
    lat = o.latitude
    lon = o.longitude
    slat = o.latitude_errors['uncertainty']*111.2
    slon = o.longitude_errors['uncertainty']*111.2*cos(lat*pi/180)
    sh = (slat+slon)/2
    z = o.depth/1000.  # in km
    sz = o.depth_errors['uncertainty']/1000
    t = o.time

    if doadd:
        with open(hashfile, 'r') as hf:
            lines = hf.readlines()
        oldstations = [line.split()[0] for line in lines if line != 0]
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
    logging.info('Skipping %s previously picked stations: %s',
                 onsta, ', '.join(oldstations))
    for ii, station in enumerate(sorted(distance, key=distance.get)):
        dist = distance[station]
        if station not in oldstations:
            logging.info('\nThis is station %s, %s of %s', station, ii, nsta)
            oldstations.append(station)
            pick = hu.PickRequest(traces, station, event, dist)
            pol, qp, leave = hu.AskForPolarity(traces, station, event, pick,
                                               dist)
            if pol in ['+', '-']:
                sp = hu.get_spratio(traces, station, pick, dist)
                outline = '{:5s} {:1s} {:1d} {:08.3f}\n'.format(
                        station, pol, qp, sp)
                with open(hashfile, 'a') as hf:
                    hf.writelines(outline)
            if leave:
                break

# Run HASH #
############
ans = []
while ans not in ['y', 'n']:
    ans = input('Want to run HASH now? (y/n)')
if ans == 'y':
    hu.create_HASH_runfile(ctrlfile, hashfile, params, ID)
    nobs = sum(1 for line in open(hashfile)) - 1
    strike, dip, rake, misfit = hu.RunHASH(ctrlfile)
    print('Strike: Dip: Rake: RMS: (deg)')
    print('{: 6.0f}  {: 3.0f}  {: 4.0f}  {: 3.0f}'.format(
        strike, dip, rake, misfit))
    new_hypfile = hu.WriteMechanismToHypfile(hypfile, strike, dip, rake,
                                             misfit, nobs, resultdir)
    hu.PlotMechanism(resultdir, new_hypfile, ID)


# the terminal will be open/activated the whole time
# has to be killed manually with the except clause!!
#    try:
#        cmd = ('ID=$(xdotool getactivewindow); while true; do sleep 0.2; ' +
#               'xdotool windowfocus $ID; done;')
#        p   = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
#        main(sID)
#    except Exception as e:
#        print(e)
#        cmd = 'pkill -f xdotool'
#        p   = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
#        exit()
