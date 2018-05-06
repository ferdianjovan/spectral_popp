#!/usr/bin/env python

from scipy.stats import gamma


# The parameter of Poisson distribution represented as a gamma distribution.
# The gamma distribution is represented by shape/alpha and rate/beta parameters (opposing
# the standard gamma distribution which is represented with shape and scale).
class Rate(object):

    def __init__(self, interval=1):
        self.reset()
        self.interval = interval

    def reset(self):
        self.beta = 1.1
        self.alpha = 1.1
        self.mode = self._mode(self.alpha, self.beta)
        self.mean = gamma.mean(self.alpha, scale=1/float(self.beta))

    # set the rate (gamma) distribution using point estimate rate and beta parameter
    # should not be used excessively
    def set_rate(self, rate, beta=-1, map_estimate=True):
        if rate < 0:
            return
        if beta != -1:
            self.beta = beta
        if map_estimate:
            self.mode = rate
            self.alpha = (rate * self.beta) + 1
            self.mean = gamma.mean(self.alpha, scale=1/float(self.beta))
        else:
            self.mean = rate
            self.alpha = rate * self.beta
            self.mode = self._mode(self.alpha, self.beta)

    # get the point estimate of the rate
    def get_rate(self, map_estimate=True):
        if map_estimate:
            return self.mode
        else:
            return self.mean

    # get the mode of gamma distribution
    def _mode(self, alpha, beta):
        if alpha >= 1:
            return (alpha - 1) / float(beta)
        else:
            return -1.0

    # get the specified percentile of the rate
    def get_rate_percentile(self, percentile):
        return gamma.ppf(percentile, self.alpha, scale=1/float(self.beta))

    # get the upper bound of the percentile of the rate (default = 0.95)
    def upper_end(self, percentile=0.95):
        return self.get_rate_percentile(percentile)

    # get the lower bound of the percentile of the rate (default = 0.05)
    def lower_end(self, percentile=0.05):
        return self.get_rate_percentile(percentile)

    # posterior distribution of the rate using conjugacy between Poisson-Gamma
    def update_rate(self, data):
        self.alpha += sum(data)
        self.beta += len(data) * self.interval
        self.mode = self._mode(self.alpha, self.beta)
        self.mean = gamma.mean(self.alpha, scale=1/float(self.beta))
