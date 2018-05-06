#!/usr/bin/env python

import os
import sys
import yaml
import time
import numpy
import datetime
import argparse
from matplotlib import pyplot

source_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(
    "/".join(source_path.split("/")[:-1]), "src"
))

from spectral_popp import PeriodicPoissonProcess, SpectralPoissonProcess


# Detection counter / process takes sensor data (leg, upper_body, or scene) to build a
# Poisson process. As the data are clustered into regions where the detections were made,
# the DetectionCounter has each process for each region.
class DetectionCounter(object):

    global source_path

    # @counter_type is the file name (leg, upper_body, scene) without ".yaml"
    # @increment is smallest interval (in seconds) where only one event happening is
    # allowed. For each minute increment, it is expected to have a rate (class Rate)
    # @periodic_cycle is the periodicity imposed to the Poisson process
    # @spectral_model is the option to go with periodic Poisson model or
    # Spectral-Poisson model.
    def __init__(
        self, counter_type, increment=60, periodic_cycle=1440*60, spectral_model=False
    ):
        print "%s process has %d second increment with %d periodic cycle..." % (
            counter_type, increment, periodic_cycle
        )
        self.counter_type = counter_type
        self.time_increment = increment
        path = os.path.join("/".join(source_path.split("/")[:-1]), "data")
        regions = yaml.load(open(path+"/regions.yaml", "r"))
        print "Process are built for %d regions" % len(regions)
        if spectral_model:
            self.process = {
                str(region): SpectralPoissonProcess(
                    increment, periodic_cycle,
                    path_to_db=os.path.join(
                        "/".join(source_path.split("/")[:-1]),
                        "db/%s/%s" % (region, counter_type)
                    )
                ) for region in regions
            }
        else:
            self.process = {
                str(region): PeriodicPoissonProcess(
                    increment, periodic_cycle,
                    path_to_db=os.path.join(
                        "/".join(source_path.split("/")[:-1]),
                        "db/%s/%s" % (region, counter_type)
                    )
                ) for region in regions
            }

    # retrieve stored rates to construct the poisson process
    def retrieve_from_db(self):
        print(
            "Retrieving %s process from db folder. It may take a while..." % self.counter_type
        )
        for region in self.process:
            print "Region", region
            self.process[region].retrieve_from_db()

    # get detection data from files
    def get_detection_data(self, start_time, end_time):
        path = os.path.join("/".join(source_path.split("/")[:-1]), "data")
        present_activity = yaml.load(open(path+"/present_activity.yaml", "r"))
        detections = yaml.load(open(path+"/%s.yaml" % self.counter_type, "r"))
        region_detection = {region: dict() for region in present_activity.keys()}
        # getting all detected and undetected activity from @self.counter_type
        # detectors given the activity was present
        for region in present_activity.keys():
            for timestamp in present_activity[region]:
                if timestamp >= start_time and timestamp < end_time:
                    region_detection[region][timestamp] = int(timestamp in detections[region])
        # getting all detected and undetected activity from @self.counter_type
        # detectors given the activity was absent
        absent_activity = yaml.load(open(path+"/absent_activity.yaml", "r"))
        for region in absent_activity.keys():
            if region not in region_detection:
                region_detection[region] = dict()
            for timestamp in absent_activity[region]:
                if timestamp >= start_time and timestamp < end_time:
                    region_detection[region][timestamp] = int(timestamp in detections[region])
        return region_detection

    # Estimate the rate function using detection count data from @start_time to @end_time
    # For all available regions
    def learn(self, start_time, end_time):
        print(
            "Estimating %s processes for each region from %s to %s" % (
                self.counter_type,
                datetime.datetime.fromtimestamp(start_time),
                datetime.datetime.fromtimestamp(end_time)
            )
        )
        detection_count = self.get_detection_data(start_time, end_time)
        for region in self.process:
            print "Region", region
            self.process[region].update(detection_count[region])
            self.process[region].store_to_db()

    # Plot arrival rate of the Poisson as a function of time. Point estimates
    # are used. Upper bound is shown
    def plot_per_region(self, region):
        if isinstance(self.process[region], PeriodicPoissonProcess().__class__):
            name = "Periodic Poisson Process for %s detections" % self.counter_type
        else:
            name = "Spectral-Poisson Process for %s detections" % self.counter_type
        mean = self.process[region].retrieve_full_periodic(mean=True)
        upper = self.process[region].retrieve_full_periodic(upper_bound=True)
        ordered_mean = list()
        ordered_upper = list()
        for key in sorted(mean.keys()):
            ordered_mean.append(mean[key])
            ordered_upper.append(upper[key])

        x = numpy.arange(len(mean))
        line = pyplot.plot(x, ordered_mean, "-", color="r", label="Mean")
        pyplot.setp(line, linewidth=3)
        line = pyplot.plot(x, ordered_upper, "--", color="b", label="Upper Bound")
        pyplot.setp(line, linewidth=2)
        pyplot.title("%s for Region %s" % (name, region), fontsize=30)
        pyplot.xlabel(
            "%d second time increment with %d periodic cycle" % (
                self.time_increment, self.process[region].periodic_cycle
            ), fontsize=15
        )
        pyplot.xticks(rotation="horizontal", fontsize=15)
        pyplot.xlim(xmax=len(x))
        pyplot.ylabel("Arrival Rate", fontsize=25)
        pyplot.yticks(fontsize=15)
        pyplot.ylim(ymax=max(upper.values())+0.2, ymin=0)
        pyplot.legend(prop={'size': 25}, loc='best')
        pyplot.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="detection_counter")
    parser.add_argument(
        "-i", dest="time_increment", default="60",
        help="Incremental time (in seconds). Default is 60 seconds"
    )
    parser.add_argument(
        "-c", dest="periodic_cycle", default="86400",
        help="Desired periodic cycle (in seconds). Default is daily (86400 seconds)"
    )
    parser.add_argument(
        "-l", dest="learn", default="0",
        help="Learning (1) the process from data or plotting (0) the available process from db"
    )
    parser.add_argument(
        "-m", dest="model", default="0",
        help="Periodic Poisson process (0) or Spectral-Poisson process (1)"
    )
    args = parser.parse_args()
    scene = DetectionCounter(
        "scene", int(args.time_increment),
        int(args.periodic_cycle), spectral_model=bool(int(args.model))
    )
    scene.retrieve_from_db()
    leg = DetectionCounter(
        "leg", int(args.time_increment),
        int(args.periodic_cycle), spectral_model=bool(int(args.model))
    )
    leg.retrieve_from_db()
    upper_body = DetectionCounter(
        "upper_body", int(args.time_increment),
        int(args.periodic_cycle), spectral_model=bool(int(args.model))
    )
    upper_body.retrieve_from_db()
    if int(args.learn):
        start_time = raw_input("Start learning time (format: 'YYYY MM DD hh mm', earliest available data: 2017 08 14 0 0): ")
        start_time = start_time.split(" ")
        start_time = datetime.datetime(
            int(start_time[0]), int(start_time[1]), int(start_time[2]),
            int(start_time[3]), int(start_time[4])
        )
        start_time = int(time.mktime(start_time.timetuple()))
        end_time = raw_input("End learning time (format: 'YYYY MM DD hh mm', latest available data: 2017 11 20 59 59): ")
        end_time = end_time.split(" ")
        end_time = datetime.datetime(
            int(end_time[0]), int(end_time[1]), int(end_time[2]),
            int(end_time[3]), int(end_time[4])
        )
        end_time = int(time.mktime(end_time.timetuple()))
        scene.learn(start_time, end_time)
        leg.learn(start_time, end_time)
        upper_body.learn(start_time, end_time)
    else:
        region = raw_input("Choose Regions %s: " % str(scene.process.keys()))
        scene.plot_per_region(region)
        leg.plot_per_region(region)
        upper_body.plot_per_region(region)
