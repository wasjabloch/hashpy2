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

loglevel = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
}


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    par = argparse.ArgumentParser(
        prog="hashpy2",
        description="An interactive HASH wrapper",
        epilog=(
            "Authors: Wasja Bloch "
            + "(wasja@gfz-potsdam.de) "
            + "and Lukas Lehmann "
            + "(luklehma@uni-potsdam.de)"
        ),
    )

    par.add_argument(
        "ID",
        type=str,
        nargs="?",
        help="Origin time of event in format: YYYYMMDDhhmmss",
        default=None,
    )
    par.add_argument(
        "--config",
        type=str,
        help="Configuration file [config.yaml]",
        default="config.yaml",
    )
    par.add_argument(
        "--setup",
        action="store_true",
        help=("Print configuration tipps. " + "Produce default_config.yaml."),
    )

    print("This is hashpy2")

    args = par.parse_args()
    ID = args.ID

    if args.setup:
        defaultf = "default_config.yaml"
        hu.print_config_tipps()
        hu.write_default_config(defaultf)
        sys.exit("Default configuration file written to: " + defaultf)

    if not ID:
        par.print_help()
        sys.exit()

    #  Read config file
    configfile = args.config
    with open(configfile, "r") as stream:
        params = yaml.safe_load(stream)
    logging.basicConfig(level=loglevel[params["LOGGING"]])

    #  Set directories
    resultdir = params["RESULTS"] + "/"
    hashfile = resultdir + ID + ".pol.hash"
    ctrlfile = resultdir + ID + ".inp"
    wvdir = params["WAVEFORMS"] + "/"
    stationfile = params["STATIONS"]
    hypfiles = glob(params["CATALOG"])

    #  Ask what to do if .pol.hash is present
    dowaveforms = True
    doadd = False

    if isfile(hashfile):
        nobs = sum(1 for line in open(hashfile)) - 1
        print("Found existing HASH input file:" + hashfile)
        print("It has " + str(nobs) + " observations.")
        print("What should I do?")
        msg = "(c)ompute focal mechanism\n" + "(a)dd missing stations?\n"
        ans = []
        while ans not in ["c", "a"]:
            ans = input(msg)
            if ans == "c":
                dowaveforms = False
            elif ans == "a":
                doadd = True


    #  Get event information and seismic traces
    locformat = None
    if hypfiles[0].endswith(".hyp"):
        logging.info("Detected NonLinLoc Hypfile format")
        locformat = "NLL"
        timestamp = hu.time_from_ID(ID)
        event, hypfile = hu.get_event_nll(timestamp, hypfiles)

    elif hypfiles[0].endswith(".pha") or hypfiles[0].endswith(".dat"):
        logging.info("Detected HypoDD Phasefile format")
        locformat = "HYPODD"
        event, timestamp = hu.get_event_hypodd(ID, hypfiles[0])
        
    else:
        msg = f"Unknown catalog file ending: " + hypfiles
        raise ValueError(msg)

    # Do waveforms #
    ################
    if dowaveforms:
        
        logging.info("Reading waveforms...")
        if params["FORMAT"] == "FUB2018":
            traces = hu.make_traces_nll(timestamp, wvdir)
        elif params["FORMAT"] == "UBC2022":
            traces = hu.make_traces_hypodd(timestamp, wvdir)
        else:
            msg = f"Unknown waveform database format: " + params["FORMAT"]
            raise ValueError(msg)
            

        stations = set([str(t.stats.station) for t in traces])
        if locformat == "NLL":
            distance = hu.dist2stations_nll(stations, event, stationfile)
        elif locformat == "HYPODD":
            distance = hu.dist2stations_hypodd(stations, event, stationfile)

        logging.info("Found waveforms for %s stations", len(stations))
        logging.info("and %s of them have available coordinates\n", len(distance))

        o = event.origins[0]
        lat = o.latitude
        lon = o.longitude
        try:
            slat = o.latitude_errors["uncertainty"] * 111.2
            slon = o.longitude_errors["uncertainty"] * 111.2 * cos(lat * pi / 180)
            sz = o.depth_errors["uncertainty"] / 1000
        except TypeError:
            logging.warning("No location uncertainty set. Setting 0")
            slat = 0
            slon = 0
            sz = 0
        sh = (slat + slon) / 2
        z = o.depth / 1000.0  # in km
        t = o.time

        if doadd:
            with open(hashfile, "r") as hf:
                lines = hf.readlines()
            oldstations = [line.split()[0] for line in lines if line != 0]
        else:
            header = (
                "{:04d} {:02d} {:02d} {:02d} {:02d} {:05.2f} "
                + "{:010.6f} {:010.6f} {:09.6f} {:05.2f} "
                + "{:05.2f} {:}\n"
            ).format(
                t.year,
                t.month,
                t.day,
                t.hour,
                t.minute,
                t.second + t.microsecond * 1e-6,
                lat,
                lon,
                z,
                sh,
                sz,
                ID,
            )
            with open(hashfile, "w") as hf:
                hf.writelines(header)
            oldstations = []

        nsta = len(distance)
        onsta = len(oldstations)
        logging.info(
            "Skipping %s previously picked stations: %s", onsta, ", ".join(oldstations)
        )
        for ii, station in enumerate(sorted(distance, key=distance.get)):
            dist = distance[station]
            if station not in oldstations:
                logging.info("\nThis is station %s, %s of %s", station, ii, nsta)
                oldstations.append(station)
                pick = hu.PickRequest(traces, station, event, dist)
                try:
                    pol, qp, leave = hu.AskForPolarity(traces, station, event, pick, dist)
                except IndexError:
                    continue
                if pol in ["+", "-"]:
                    sp = hu.get_spratio(traces, station, pick, dist)
                    outline = "{:5s} {:1s} {:1d} {:08.3f}\n".format(
                        station, pol, qp, sp
                    )
                    with open(hashfile, "a") as hf:
                        hf.writelines(outline)
                if leave:
                    break

    # Run HASH #
    ############
    ans = []
    while ans not in ["y", "n"]:
        ans = input("Want to run HASH now? (y/n)")
    if ans == "y":
        hu.create_HASH_runfile(ctrlfile, hashfile, params, ID)
        nobs = sum(1 for line in open(hashfile)) - 1
        strike, dip, rake, misfit = hu.RunHASH(ctrlfile)
        print("Strike: Dip: Rake: RMS: (deg)")
        print(
            "{: 6.0f}  {: 3.0f}  {: 4.0f}  {: 3.0f}".format(strike, dip, rake, misfit)
        )

        if locformat == "NLL":
            new_hypfile = hu.WriteMechanismToHypfile(
                hypfile, strike, dip, rake, misfit, nobs, resultdir
            )
            hu.PlotMechanismHyp(resultdir, new_hypfile, ID)
            
        if locformat == "HYPODD":
            rfn = resultdir + "/" + ID + ".focmec.dat"
            hu.WriteMechanismToResultfile(
                event, strike, dip, rake, misfit, nobs, rfn
            )
            hu.PlotMechanismRes(rfn, resultdir, ID)
            

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
