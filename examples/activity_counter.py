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


# Activity counter / process takes ground truth activity to build a
# Poisson process. As the data are clustered into regions where the detections were made,
# the ActivityCounter has each process for each region.
class ActivityCounter(object):

    global source_path

    # @increment is smallest interval (in seconds) where only one event happening is
    # allowed. For each minute increment, it is expected to have a rate (class Rate)
    # @periodic_cycle is the periodicity imposed to the Poisson process
    # @spectral_model is the option to go with periodic Poisson model or
    # Spectral-Poisson model.
    def __init__(
        self, increment=60, periodic_cycle=1440*60, spectral_model=False
    ):
        print "Activity process has %d second increment with %d periodic cycle..." % (
            increment, periodic_cycle
        )
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
                        "db/%s/activity" % region
                    )
                ) for region in regions
            }
        else:
            self.process = {
                str(region): PeriodicPoissonProcess(
                    increment, periodic_cycle,
                    path_to_db=os.path.join(
                        "/".join(source_path.split("/")[:-1]),
                        "db/%s/activity" % region
                    )
                ) for region in regions
            }

    # retrieve stored rates to construct the activity poisson process
    def retrieve_from_db(self):
        print("Retrieving activity process from db folder. It may take a while...")
        for region in self.process:
            print "Region", region
            self.process[region].retrieve_from_db()

    # get ground truth activity (present and absent) from files
    def get_activity_data(self, start_time, end_time):
        path = os.path.join("/".join(source_path.split("/")[:-1]), "data")
        present_activity = yaml.load(open(path+"/present_activity.yaml", "r"))
        region_activity = {region: dict() for region in present_activity.keys()}
        # getting all present activity from activity yaml file
        for region in present_activity.keys():
            region_activity[region] = {
                timestamp: 1 for timestamp in present_activity[region] if timestamp >= start_time and timestamp < end_time
            }
        # getting all absent activity from absent activity yaml file
        absent_activity = yaml.load(open(path+"/absent_activity.yaml", "r"))
        for region in absent_activity.keys():
            if region not in region_activity:
                region_activity[region] = dict()
            region_activity[region].update({
                timestamp: 0 for timestamp in absent_activity[region] if timestamp >= start_time and timestamp < end_time
            })
        return region_activity

    # Estimate the rate function using activity count data from @start_time to @end_time
    # For all available regions
    def learn(self, start_time, end_time):
        print(
            "Estimating activity processes for each region from %s to %s" % (
                datetime.datetime.fromtimestamp(start_time),
                datetime.datetime.fromtimestamp(end_time)
            )
        )
        activity_count = self.get_activity_data(start_time, end_time)
        for region in self.process:
            print "Region", region
            self.process[region].update(activity_count[region])
            self.process[region].store_to_db()

    # Plot arrival rate of the Poisson as a function of time. Point estimates
    # are used. Upper bound is shown
    def plot_per_region(self, region):
        if isinstance(self.process[region], PeriodicPoissonProcess().__class__):
            name = "Periodic Poisson Process"
        else:
            name = "Spectral-Poisson Process"
        # @mean, @upper become time series data. However they are somehow not
        # timely ordered. Sorting is needed
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
    parser = argparse.ArgumentParser(prog="activity_counter")
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
    ac = ActivityCounter(
        int(args.time_increment), int(args.periodic_cycle), spectral_model=bool(int(args.model))
    )
    ac.retrieve_from_db()
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
        ac.learn(start_time, end_time)
    else:
        region = raw_input("Choose Regions %s: " % str(ac.process.keys()))
        ac.plot_per_region(region)
