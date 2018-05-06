#!/usr/bin/env python

import os
import yaml
import copy
from rate import Rate
from fourier import rectify_signal, reconstruct_signal


# (Nonhomogeneous) Poisson process representation with arrival rate represented as gamma distribution
class PoissonProcess(object):

    # @increment is smallest interval (in seconds) where only one event happening is
    # allowed. For each minute increment, it is expected to have a rate (class Rate)
    # @path_to_db is the path to database folder storing this class.
    def __init__(self, increment=1, path_to_db="", db_name="poisson_process"):
        self.poisson = dict()
        self.increment = increment
        if path_to_db == "":
            path_to_db = os.path.join(os.getcwd(), db_name)
        else:
            path_to_db = os.path.join(path_to_db, db_name)
        self._path_to_db = os.path.join(path_to_db, str(self.increment))
        if not os.path.exists(self._path_to_db):
            os.makedirs(self._path_to_db)

    # convert the start time of an event into an @increment step interval
    def _convert_time(self, start_time):
        return (start_time / self.increment) * self.increment

    def default_rate(self):
        return Rate()

    # manually set the rate of the smallest interval at @start_time with @rate
    # @rate is a class Rate
    def set_rate_at(self, start_time, rate):
        start_time = self._convert_time(start_time)
        self.poisson[start_time] = rate

    # get the rate (class Rate) at specific time
    def get_rate_at(self, start_time):
        start_time = self._convert_time(start_time)
        if start_time not in self.poisson:
            rate = self.default_rate()
        else:
            rate = self.poisson[start_time]
        return rate

    # @counts is in the form of {@timestamp: @counts}
    # @count, can be (the smallest) 0 or 1, or natural number.
    # @timestamp is the specific time corresponding to @count
    def update(self, counts):
        for start_time, count in counts.iteritems():
            start_time = self._convert_time(start_time)
            if start_time not in self.poisson:
                self.poisson[start_time] = self.default_rate()
            self.poisson[start_time].update_rate([count])

    # store specific rate of the poisson process to a file
    # KeyError is possible!
    def _store(self, start_time):
        start_time = self._convert_time(start_time)
        rate = self.poisson[start_time]
        path = os.path.join(self._path_to_db, "%s.yaml" % start_time)
        with open(path, "w") as f:
            f.write(yaml.dump({"alpha": rate.alpha, "beta": rate.beta}))

    # store the poisson process by storing its rates at each increment
    def store_to_db(self):
        print("Storing this Poisson process with %d data" % len(self.poisson))
        for start_time in self.poisson.iterkeys():
            self._store(start_time)

    # retrieve stored rates to construct the poisson process
    def retrieve_from_db(self):
        start_times = [
            f for f in os.listdir(self._path_to_db) if os.path.isfile(
                os.path.join(self._path_to_db, f)
            )
        ]
        for start_time in start_times:
            data = yaml.load(
                open(os.path.join(self._path_to_db, start_time), "r")
            )
            start_time = int(start_time.split(".")[0])
            self.poisson[start_time] = self.default_rate()
            self.poisson[start_time].set_rate(
                (data["alpha"] - 1) / float(data["beta"]),
                beta=float(data["beta"])
            )

        retrieved = len(start_times) > 0
        if retrieved:
            print("%d new poisson distributions are obtained from db..." % len(start_times))
        return retrieved

    # get point estimates of arrival rates from @start_time to @end_time
    # point estimates can be the MAP hypothesis (default), mean expectation,
    # the upper bound of the rate (Gamma) distribution, or the lower bound of
    # the rate.
    def retrieve(
        self, start_time, end_time, mean=False,
        upper_bound=False, lower_bound=False
    ):
        """ retrieve poisson distribution from specified start_time until specified end_time
        """
        # convert start_time and end_time to the closest time range (in minutes)
        start_time = self._convert_time(start_time)
        end_time = self._convert_time(end_time)
        print(
            "Retrieving arrival rate from %d to %d" % (start_time, end_time)
        )
        result = dict()
        while start_time < end_time:
            try:
                poisson = self.poisson[start_time]
            except:
                poisson = self.default_rate()
            # upper trumphs lower
            if upper_bound:
                result[start_time] = poisson.upper_end()
            elif lower_bound:
                result[start_time] = poisson.lower_end()
            elif mean:
                result[start_time] = poisson.get_rate(map_estimate=True)
            else:
                result[start_time] = poisson.get_rate(map_estimate=False)
            start_time = start_time + self.increment
        return result


# Periodic Poisson process is a nonhomogeneous Poisson process where the rate
# function repeats itself after some delta-time
class PeriodicPoissonProcess(PoissonProcess):

    # @periodic_cycle is the fixed periodicity imposed to the Poisson process
    def __init__(
        self, increment=1, periodic_cycle=3600,
        path_to_db="", db_name="poisson_process"
    ):
        self._pivot_time = None
        self.periodic_cycle = periodic_cycle
        super(
            PeriodicPoissonProcess, self
        ).__init__(increment, path_to_db, db_name+"_"+str(periodic_cycle))

    # converting start_time to relative time from the periodic cycle
    def _relative_start_time(self, start_time):
        if self._pivot_time is None:
            self._pivot_time = start_time
        delta = (start_time - self._pivot_time) % self.periodic_cycle
        return self._pivot_time + delta

    # get the rate (class Rate) at specific time
    def get_rate_at(self, start_time):
        start_time = self._relative_start_time(start_time)
        return super(PeriodicPoissonProcess, self).get_rate_at(start_time)

    # manually set the rate of the smallest interval at @start_time with @rate
    # @rate is a class Rate.
    def set_rate_at(self, start_time, rate):
        start_time = self._relative_start_time(start_time)
        return super(PeriodicPoissonProcess, self).set_rate_at(start_time, rate)

    # @counts is in the form of {@timestamp: @counts}
    # @count, can be (the smallest) 0 or 1, or natural number.
    # @timestamp is the specific time corresponding to @count
    def update(self, counts):
        for start_time, count in counts.iteritems():
            start_time = self._relative_start_time(start_time)
            start_time = self._convert_time(start_time)
            if start_time not in self.poisson:
                self.poisson[start_time] = self.default_rate()
            self.poisson[start_time].update_rate([count])

    # store the poisson process by storing its rates at each increment
    def store_to_db(self):
        print("Storing this Poisson process with %d data" % len(self.poisson))
        for start_time in self.poisson.iterkeys():
            start_time = self._relative_start_time(start_time)
            super(PeriodicPoissonProcess, self)._store(start_time)

    # retrieve stored rates to construct the poisson process
    def retrieve_from_db(self):
        retrieved = super(PeriodicPoissonProcess, self).retrieve_from_db()
        if retrieved:
            self._pivot_time = min(self.poisson.keys())
        return retrieved

    # get point estimates of arrival rates from @start_time to @end_time
    # point estimates can be the MAP hypothesis (default), mean expectation,
    # the upper bound of the rate (Gamma) distribution, or the lower bound of
    # the rate.
    def retrieve(
        self, start_time, end_time, mean=False,
        upper_bound=False, lower_bound=False
    ):
        """ retrieve poisson distribution from specified start_time until specified end_time
        """
        # convert start_time and end_time to the closest time range (in minutes)
        start_time = self._convert_time(start_time)
        end_time = self._convert_time(end_time)
        print(
            "Retrieving arrival rate from %d to %d" % (start_time, end_time)
        )
        result = dict()
        while start_time < end_time:
            relative_start = self._relative_start_time(start_time)
            try:
                poisson = self.poisson[relative_start]
            except:
                poisson = self.default_rate()
            # upper trumphs lower
            if upper_bound:
                result[start_time] = poisson.upper_end()
            elif lower_bound:
                result[start_time] = poisson.lower_end()
            else:
                result[start_time] = poisson.get_rate(map_estimate=(not mean))
            start_time = start_time + self.increment
        return result

    # get a complete cycle of point estimates of arrival rates
    # point estimates can be the MAP hypothesis (default), mean expectation,
    # the upper bound of the rate (Gamma) distribution, or the lower bound of
    # the rate.
    def retrieve_full_periodic(
        self, mean=False, upper_bound=False, lower_bound=False
    ):
        if self._pivot_time is not None:
            start_time = self._pivot_time
        else:
            start_time = 0
        end_time = start_time + self.periodic_cycle
        return self.retrieve(
            start_time, end_time, mean, upper_bound, lower_bound
        )


# Spectral-Poisson process is a nonhomogeneous Poisson process where the rate
# function is modelled with the help of Fourier transformation
class SpectralPoissonProcess(PeriodicPoissonProcess):

    # spectral-poisson has 3 additional processes apart from the @self.poisson
    # to represent the Fourier transformation of @self.poisson.
    # @self.poisson is the original rate function of Poisson, whereas
    # @self._spectral is where the point estimate of the rate function is
    # transformed via Fourier to obtain new smother rate function.
    # @self._upper and @self._lower are the representative of the upper bound
    # and the lower bound of the rate function transformed using Fourier.
    def __init__(
        self, increment=1, periodic_cycle=3600,
        path_to_db="", db_name="poisson_process"
    ):
        self._spectral = dict()
        super(
            SpectralPoissonProcess, self
        ).__init__(increment, periodic_cycle, path_to_db, db_name)

    # retrieve stored rates to construct the poisson and spectral processes
    def retrieve_from_db(self):
        retrieved = super(SpectralPoissonProcess, self).retrieve_from_db()
        if self._pivot_time is not None:
            self.fourier_transform()
        return retrieved

    # transform the rate function in @self.poisson using Fourier to create
    # the spectral model of the rate function
    def fourier_transform(self):
        if self._pivot_time is not None:
            start_time = self._pivot_time
        else:
            start_time = 0
        end_time = start_time + self.periodic_cycle
        stamped_rates = super(
            SpectralPoissonProcess, self
        ).retrieve(start_time, end_time, True, False, False)
        rates, _ = reconstruct_signal(
            [stamped_rates[timestamp] for timestamp in sorted(stamped_rates.keys())]
        )
        rates = rectify_signal(rates, low_thres=0.001)
        self._spectral = dict()
        for index, timestamp in enumerate(sorted(stamped_rates.keys())):
            self._spectral[timestamp] = copy.deepcopy(self.get_rate_at(timestamp))
            self._spectral[timestamp].set_rate(
                rates[index], beta=self._spectral[timestamp].beta,
                map_estimate=False
            )

    # @counts is in the form of {@timestamp: @counts}
    # @count, can be (the smallest) 0 or 1, or natural number.
    # @timestamp is the specific time corresponding to @count
    def update(self, counts):
        super(SpectralPoissonProcess, self).update(counts)
        self.fourier_transform()

    # get point estimates of arrival rates from @start_time to @end_time
    # point estimates can be the MAP hypothesis (default), mean expectation,
    # the upper bound of the rate (Gamma) distribution, or the lower bound of
    # the rate.
    # The retrieval of point estimates are based on the spectral model.
    def retrieve(
        self, start_time, end_time, mean=False,
        upper_bound=False, lower_bound=False
    ):
        """
            retrieve spectral-poisson distribution from
            specified start_time until specified end_time
        """
        # convert start_time and end_time to the closest time range (in minutes)
        start_time = self._convert_time(start_time)
        end_time = self._convert_time(end_time)
        print(
            "Retrieving arrival rate from %d to %d" % (start_time, end_time)
        )
        result = dict()
        while start_time < end_time:
            relative_start = self._relative_start_time(start_time)
            try:
                poisson = self._spectral[relative_start]
            except:
                poisson = self.default_rate()
            # upper trumphs lower
            if upper_bound:
                result[start_time] = poisson.upper_end()
            elif lower_bound:
                result[start_time] = poisson.lower_end()
            else:
                result[start_time] = poisson.get_rate(map_estimate=(not mean))
            start_time = start_time + self.increment
        return result
