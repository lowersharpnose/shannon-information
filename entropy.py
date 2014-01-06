    #  -----------------------------------------------------------------------
    #  Name:        entropy
    #  Purpose:
    #  Version 1.01
    #
    #  Author:      JdP
    #  -----------------------------------------------------------------------

import collections
import math

    # HX = HX + Px * math.log(1.0/Px, 2)
    # is equivalent to
    # HX = HX - Px * math.log(Px, 2)
    # The latter, more efficient form is used.
    # Where bins are used in the production of histograms from discrete data
    # that need to be approximated to probability density function, then the
    # binSize used is passed in.
    # This results in the term log(binSize) being added to the summed entropy
    # value (base 2 log).
    # http://www2.warwick.ac.uk/fac/soc/economics/staff/academic/wallis/publications/entropy.pdf


def EntropyFromProbabilityDistribution(distribution):

    """Calculate the entropy of the passed list of probabilities.

    distribution is an iterable list of probabilities.
    All values should be positive.
    The sum of the probabilities should be equal to 1.0 to 2 decimal places
    This function returns the entropy (HX) as a number of bits to 3 sig figs.

    The error return is -1"""

    if not isinstance(distribution, collections.Iterable):
        HX = -1

    # Check the sum of probabilities is 1.0
    elif (round(sum(distribution), 2) != 1.0):
        HX = -1

    # Check for negative probabilities
    elif (min(distribution) < 0.0):
        HX = -1

    else:
        HX = 0

        for i in range(0, len(distribution)):
            Px = distribution[i]

            if (Px > 0):
                HX = HX - Px * math.log(Px, 2)

    return (HX)


def DifferentialEntropyFromProbabilityDistribution(distribution, binSize):

    """Calculate the entropy of the passed list of probabilities.

    distribution is an iterable list of probabilities.
    All values should be positive.
    The sum of the probabilities should be equal to 1.0 to 2 decimal places
    This function returns the entropy (HX) as a number of bits to 3 sig figs.

    The error return is -1"""

    if not isinstance(distribution, collections.Iterable):
        HX = -1

    # Check the sum of probabilities is 1.0
    elif (round(sum(distribution) * binSize, 2) != 1.0):
        HX = -1

    # Check for negative probabilities
    elif (min(distribution) < 0.0):
        HX = -1

    elif (binSize < 0.0):
        HX = -1

    else:
        HX = 0

        if (len(distribution) > 1):
            for i in range(0, len(distribution)):
                Px = distribution[i] * binSize

                if (Px > 0):
                    HX = HX - Px * math.log(Px, 2)

            HX = HX + math.log(binSize, 2)

    return (HX)


def EntropyFromFrequencyDistribution(distribution):

    """Calculate the entropy of the passed list of frequencies.

    distribution is an iterable list of frequencies.
    This function returns the entropy (HX) as a number of bits to 3 sig figs.

    The error return is -1"""

    if not isinstance(distribution, collections.Iterable):
        HX = -1

    elif len(distribution) < 1:
        HX = -1

    elif (min(distribution) < 0.0):
        HX = -1

    else:

        HX = 0

        totalFrequency = sum(distribution)

        for i in range(0, len(distribution)):
            Px = float(distribution[i]) / totalFrequency

            if (Px > 0.0):
                HX = HX - Px * math.log(Px, 2)

    return (HX)


def DifferentialEntropyFromFrequencyDistribution(distribution, binSize):

    """Calculate the entropy of the passed list of frequencies.

    distribution is an iterable list of frequencies.
    This function returns the entropy (HX) as a number of bits to 3 sig figs.

    The error return is -1"""

    if not isinstance(distribution, collections.Iterable):
        HX = -1

    elif len(distribution) < 1:
        HX = -1

    elif (min(distribution) < 0.0):
        HX = -1

    elif (binSize <= 0.0):
        HX = -1

    else:

        HX = 0

        if (len(distribution) > 1):

            totalFrequency = sum(distribution)

            for i in range(0, len(distribution)):
                Px = float(distribution[i]) / totalFrequency

                if (Px > 0.0):
                    HX = HX - Px * math.log(Px, 2)

            HX = HX + math.log(binSize, 2)

    return (HX)


def EntropyFromSampleDistribution(distribution):

    """Calculate the entropy of the passed list of sample outcomes.

    distribution is an iterable list of outcomes/values/items/states/etc.
    e.g     [1,2,5,2,2,1,4,4,2,2,2,2]
            ["A","B","C","B","B","A","D","D","B","B","B","B"]
            ["H","T","H","T",,"H","H","H","H","H"]
    This function returns the entropy (HX) as a number of bits to 3 sig figs.

    The error return is -1"""

    if not isinstance(distribution, collections.Iterable):
        HX = -1

    elif len(distribution) < 1:
        HX = -1

    else:
        sampleSize = len(distribution)

        HX = 0

        for frequency in collections.Counter(sorted(distribution)).values():
            Px = float(frequency) / sampleSize

            if (Px > 0.0):
                HX = HX - Px * math.log(Px, 2)

    return (HX)


def EntropyHXFromFrequencyDistribution(distribution, rows):

    """returns HX = entropy of the row totals

    distribution is an iterable list of frequencies and rows is the
    number of rows in the data.  If the list comprises N values, then
    N = rows x columns.
    The error return is both set to -1"""

    HX = -1

    rows = max(rows, 1)
    columns = int(len(distribution) / rows)

    # check that rows is a factor of the size of the distribution - ensures
    # we have a proper AxB matrix with integer numbers of rows and columns
    if (columns * rows == len(distribution)):

        listx = []

        for row in range(0, rows):
            listx.append(0)
            for col in range(0, columns):
                listx[row] = listx[row] + distribution[row * columns + col]

        HX = EntropyFromFrequencyDistribution(listx)

    return (HX)


def EntropyHYFromFrequencyDistribution(distribution, rows):

    """returns HY = entropy of the column totals

    distribution is an iterable list of frequencies and rows is the
    number of rows in the data.  If the list comprises N values, then
    N = rows x columns.
    The error return is both set to -1"""

    HY = -1

    rows = max(rows, 1)
    columns = int(len(distribution) / rows)

    # check that rows is a factor of the size of the distribution - ensures
    # we have a proper AxB matrix with integer numbers of rows and columns
    if (columns * rows == len(distribution)):

        listy = []

        for col in range(0, columns):
            listy.append(0)
            for row in range(0, rows):
                listy[col] = listy[col] + distribution[row * columns + col]

            HY = EntropyFromFrequencyDistribution(listy)

    return (HY)


def EntropyValuesFromFrequencyDistribution(distribution, rows):

    """returns a six element DICTIONARY with keys HX, HY, HXY, HXgY, HYgX, IXY
    HXY = joint entropy
    HX = entropy of the row totals
    HY = entropy of the column totals
    HX|Y = conditional entropy of X given Y
    HY|X = conditional entropy of Y given X
    IXY = mutual information

    distribution is an iterable list of frequencies and rows is the
    number of rows in the data.  If the list comprises N values, then
    N = rows x columns.
    The error return is both set to -1"""

    HXgY = -1
    HYgX = -1
    IXY = -1

    HXY = EntropyFromFrequencyDistribution(distribution)

    HX = EntropyHXFromFrequencyDistribution(distribution, rows)
    if (HX != -1):
        HYgX = HXY - HX

    HY = EntropyHYFromFrequencyDistribution(distribution, rows)
    if (HY != -1):
        HXgY = HXY - HY

    if (HX != -1) and (HY != -1) and (HXY != -1):
        IXY = HX + HY - HXY

    HXY = RoundSF(HXY, 3)
    HX = RoundSF(HX, 3)
    HY = RoundSF(HY, 3)
    HXgY = RoundSF(HXgY, 3)
    HYgX = RoundSF(HYgX, 3)
    IXY = RoundSF(IXY, 3)

    rc = {'HXY': HXY, 'HX': HX, 'HY': HY, 'HX|Y': HXgY,
          'HY|X': HYgX, 'IXY': IXY}
    return (rc)


def RoundSF(num, sigfigs):

    """Round to specified number of sigfigs."""

    if num == 0:
        return (0)

    rc = round(num, -int(math.floor(math.log(abs(num), 10)) - (sigfigs - 1)))

    return (rc)


def CheckResult(result, target, sigfigs=3):

    rc = RoundSF(result, sigfigs)

    if rc == target:
        print ("OK", rc)
    else:
        print ("FAIL->", rc, " target->", target)

