""" ha1utils.pyx
This code is adapted for 80 columns
|------------------------------------------------------------------------------|
"""
import numpy
cimport numpy
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
cdef numpy.float64_t \
    max_in(numpy.ndarray[numpy.float64_t, ndim=1] arr):
    """
    Find the maximum value in arr.

    Parameters
    ----------
    arr : numpy.ndarray[numpy.float64_t, ndim=1]
        The array whose maximum value will be computed.

    Returns
    -------
    numpy.float64_t
        The maximum value in arr.
    """
    cdef numpy.float64_t maxval = arr[0]
    cdef int i, arrLen = len(arr)
    for i in range(arrLen):
        if maxval < arr[i]:
            maxval = arr[i]
    return maxval


@cython.boundscheck(False)
@cython.wraparound(False)
cdef numpy.uint32_t \
    avg(numpy.ndarray[numpy.uint32_t, ndim=1] arr):
    """
    find the average value of arr.

    Parameters
    ----------
    arr : numpy.ndarray[numpy.float64_t, ndim=1]
        The array whose average value will be computed.

    Returns
    -------
    numpy.uint32_t
        The average value of arr.
    """
    cdef numpy.uint32_t av = 0
    cdef int i, arrLen = len(arr)
    for i in range(arrLen):
        av += arr[i]
    return av/arrLen


@cython.boundscheck(False)
@cython.wraparound(False)
cdef numpy.uint32_t \
    occurence_of(numpy.float64_t val,
                 numpy.ndarray[numpy.float64_t, ndim=1] arr,
                 numpy.ndarray[numpy.uint32_t, ndim=1] out,
                 numpy.float64_t tol):
    """
    Compute number of occurences of val in arr, store the
    indices in out and return the number of occurences.

    Parameters
    ----------
    val : numpy.float64_t
        The value which the number of occurences is to be counted.
    arr : numpy.ndarray[numpy.float64_t, ndim=1]
        The input array, whose values will be compared to val.
    out : numpy.ndarray[numpy.uint32_t, ndim=1]
        The output array where the occurences will be stored.
    tol : numpy.float64_t
        The tolerance relative to the maximum value where the signal
        level is to be considered for peak values, i.e.
        abs(maxval - val) < tol.

    Returns
    -------
    numpy.uint32_t
        The number of occurences of val in arr.
    """
    cdef int i, arrLen = len(arr)
    cdef numpy.uint32_t nVals = 0, offset = 0
    for i in range(arrLen):
        if ((val - arr[i] < tol and val - arr[i] >= 0)
                or (val - arr[i] < -tol and arr[i] - val <= 0)):
            out[nVals] = offset
            nVals += 1
        offset += 1
    return nVals


@cython.boundscheck(False)
@cython.wraparound(False)
cdef numpy.uint32_t \
    peak_sep_cond(numpy.ndarray[numpy.uint32_t, ndim=1] arr):
    """
    Finds the condition for separating points into belonging to
    different peaks.

    This method is probably not general enough but it have worked
    so far.

    Parameters
    ----------
    index_darr : numpy.ndarray[numpy.uint32_t, ndim=1]
        The indices of all points that will be considered for peak values.

    Returns
    -------
    numpy.float64_t
        The contdition for separating points into belonging to
        different peaks. The contision is on the form:
        If the points separation is less than the condition, the points
        are considered to belong to the same peak.
    """
    return avg(arr)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef numpy.ndarray[numpy.uint32_t, ndim=2] \
    sort_groups(numpy.ndarray[numpy.uint32_t, ndim=1] index_darr,
                numpy.ndarray[numpy.uint32_t, ndim=1] dist,
                numpy.uint32_t peak_sep_cond):
    """
    Sorts index_darr into groups if values in dist are smaller than
    peak_sep_cond.

    Parameters
    ----------
    index_darr : numpy.ndarray[numpy.uint32_t, ndim=1]
        The indices of all points that will be considered for peak values.
    dist : numpy.ndarray[numpy.uint32_t, ndim=1]
        the distance (in terms of indices) betwee the points in index_darr.

    Returns
    -------
    numpy.ndarray[numpy.uint32_t, ndim=2]
        The indices of all points that will be considered for peak values,
        sorted into groups where the indices in each group belongs to the
        same peak.
    """
    cdef int arrLen = len(dist)
    cdef numpy.ndarray[numpy.uint32_t, ndim=2] pairs =\
        numpy.zeros((arrLen, arrLen), dtype=numpy.uint32)
    pairs[0][0] = index_darr[0]
    cdef int i, j = 1, pair = 0
    for i in range(arrLen):
        if dist[i] > peak_sep_cond:
            pair += 1
            j = 0
        pairs[pair][j] = index_darr[i+1]
        j += 1
    return pairs


@cython.boundscheck(False)
@cython.wraparound(False)
cdef numpy.ndarray[numpy.uint32_t, ndim=1] \
    find_peaks(numpy.ndarray[numpy.uint32_t, ndim=2] groups,
               numpy.ndarray[numpy.float64_t, ndim=1] signal):
    """
    Finds the index of the highest signal value in the arrays
    of groups.

    Parameters
    ----------
    groups : numpy.ndarray[numpy.uint32_t, ndim=2]
        The potential peak point indices in the signal data which belongs.
        to the same maximum value.
        Axis 0 contains each peak and axis 1 contains the data for each
        peak.
    signal : numpy.ndarray[numpy.float64_t, ndim=1]
        The signal data.

    Returns
    -------
    numpy.ndarray[numpy.uint32_t, ndim=1]
        The peak point indices in the signal data from lowest to highest.
    """
    cdef int arrLen = len(groups)
    cdef numpy.ndarray[numpy.uint32_t, ndim=1] peaks =\
        numpy.empty(arrLen, dtype=numpy.uint32)
    cdef numpy.float64_t arrMax
    cdef int i = 0, j, g_i0, g_ij, arrMaxIdx
    for i in range(arrLen):
        # look for max signal values and store index of those
        g_i0 = int(groups[i][0])
        arrMax = signal[g_i0]
        arrMaxIdx = g_i0
        for j in range(1, arrLen):
            g_ij = int(groups[i][j])
            if signal[g_ij] > arrMax:
                arrMax = signal[g_ij]
                arrMaxIdx = g_ij
            if g_ij == 0:
                break;  # skip trailing zeros
        if g_i0 == 0:
            break;  # skip empty pairs
        peaks[i] = arrMaxIdx  # index of max value found

    # remove trailing zeros
    cdef numpy.ndarray[numpy.uint32_t, ndim=1] out =\
        numpy.empty(i, dtype=numpy.uint32)
    for j in range(i):
        out[j] = peaks[j]
    return out


@cython.boundscheck(False)
@cython.wraparound(False)
cdef numpy.ndarray[numpy.uint32_t, ndim=1] \
    center_freqs(numpy.ndarray[numpy.uint32_t, ndim=1] peaks,
                 numpy.ndarray[numpy.float64_t, ndim=1] signal,
                 numpy.float64_t bwcond):
    """
    Find center frequencies.

    Parameters
    ----------
    peaks : numpy.ndarray[numpy.uint32_t, ndim=1]
        The peak point indices in the signal data from lowest to highest
        index.
    signal : numpy.ndarray[numpy.float64_t, ndim=1]
        The signal data.
    tol : numpy.float64_t
        The condition relative to the *maximum* value where the signal
        level is considered to be outside the band.

    Returns
    -------
    numpy.ndarray[numpy.uint32_t, ndim=1]
        The center frequency points for each band from lowest to highest
        index.
    """
    cdef int arrLen = len(peaks)
    cdef numpy.ndarray[numpy.uint32_t, ndim=2] pairs =\
        numpy.zeros((arrLen, 2), dtype=numpy.uint32)
    cdef int i, j = 0, lb, ub, pair = 0, start = 0
    cdef numpy.uint32_t upper, lower
    while (start < arrLen):
        lower = peaks[start]
        pairs[pair][0] = lower
        # set the upper peak of this band in case it is a single peak
        pairs[pair][1] = lower
        lb = int(lower)
        for i in range(start+1, arrLen):
            # find the next peak inside this band
            upper = peaks[i]
            ub = int(peaks[i])
            for j in range(lb, ub):
                # peaks belong to different bands if below bwcond
                if signal[j] < bwcond:
                    break
            if j == ub-1:
                # peaks belong to same band
                pairs[pair][1] = upper
                start += 1
            else:
                break
        start += 1  # already checked
        pair += 1
    cdef numpy.ndarray[numpy.uint32_t, ndim=1] av =\
        numpy.empty(pair, dtype=numpy.uint32)
    for i in range(pair):
        av[i] = avg(pairs[i])
    return av


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef numpy.ndarray[numpy.uint32_t, ndim=1] \
    find_centers(numpy.ndarray[numpy.float64_t, ndim=1] signal,
                 numpy.float64_t tol):
    """
    Compute center frequencies.

    Parameters
    ----------
    signal : numpy.ndarray[numpy.float64_t, ndim=1]
        The signal data
    tol : numpy.float64_t
        The tolerance relative to the maximum value where the signal level
        is to be considered for peak values, i.e. abs(maxval - val) < tol.

    Returns
    -------
    numpy.ndarray[numpy.uint32_t, ndim=1]
        The center frequency points for each band from lowest to highest index.
    """
    # Find maximum value
    cdef numpy.float64_t maxval = max_in(signal)

    # Compute number of occurences
    cdef numpy.ndarray[numpy.uint32_t, ndim=1] maxvalIndex =\
        numpy.empty(1000, dtype=numpy.uint32)
    cdef numpy.uint32_t nMaxvals = occurence_of(
        maxval, signal, maxvalIndex, tol)

    # computed distance to pairs
    cdef numpy.ndarray[numpy.uint32_t, ndim=1] indexDist =\
        numpy.empty(nMaxvals - 1, dtype=numpy.uint32)
    cdef int i
    for i in range(nMaxvals - 1):
        indexDist[i] = maxvalIndex[i+1] - maxvalIndex[i]

    # group into groups of potential peaks
    cdef numpy.ndarray[numpy.uint32_t, ndim=2] groups =\
        sort_groups(maxvalIndex, indexDist, peak_sep_cond(indexDist))

    # find the indices for the peaks
    cdef numpy.ndarray[numpy.uint32_t, ndim=1] peaks =\
        find_peaks(groups, signal)

    # find center frequencies given the peaks
    # 0.6 is lower limit relative to maximum value for separating bands
    cdef numpy.ndarray[numpy.uint32_t, ndim=1] center =\
        center_freqs(peaks, signal, 0.6)
    return center


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef numpy.ndarray[numpy.uint32_t, ndim=2] \
    find_bandwidths(numpy.ndarray[numpy.float64_t, ndim=1] signal,
           numpy.ndarray[numpy.uint32_t, ndim=1] centers,
           numpy.float64_t bwcond):
    """
    Compute bandwidths given center point(s) and a lower-bound condition.

    Parameters
    ----------
    signal : numpy.ndarray[numpy.float64_t, ndim=1]
        The signal data
    centers : numpy.ndarray[numpy.uint32_t, ndim=1]
        The indices of the band's center frequency points in the signal data.
    bwcond : numpy.float64_t
        The condition where the signal level is considered to be outside
        the band. This value is relative to center signal level and not
        necessarily the peak value inside the band.

    Returns
    -------
    numpy.ndarray[numpy.uint32_t, ndim=1]
        The bandwidth expressed in indices. It is up to the caller of this
        function to determine what this corresponds to.
    """
    cdef int centerLen = len(centers), arrLen = len(signal)
    cdef numpy.ndarray[numpy.uint32_t, ndim=2] bandwidths =\
        numpy.empty((centerLen, 2), numpy.uint32)
    cdef int i, j, upper, lower
    cdef numpy.float64_t cond
    for i in range(centerLen):
        upper = lower = int(centers[i])
        cond = signal[int(centers[i])]*bwcond
        for j in range(centers[i], arrLen-1):
            if signal[j+1] < cond:
                upper = j
                break
        for j in range(centers[i], 0, -1):
            if signal[j-1] < cond:
                lower = j
                break
        bandwidths[i][0] = lower
        bandwidths[i][1] = upper
    return bandwidths


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef numpy.ndarray[numpy.uint32_t, ndim=2] \
    findbw(numpy.ndarray[numpy.float64_t, ndim=1] signal,
           numpy.float64_t bwcond,
           numpy.float64_t tol):
    """
    Find bandwidth of signal.

    Parameters
    ----------
    signal : numpy.ndarray[numpy.float64_t, ndim=1]
        The signal data
    bwcond : numpy.float64_t
        The condition where the signal level is considered to be outside
        the band. This value is relative to center signal level and not
        necessarily the peak value inside the band.
    tol : numpy.float64_t
        The tolerance relative to the maximum value where the signal level
        is to be considered for peak values, i.e. abs(maxval - val) < tol.

    Returns
    -------
    numpy.ndarray[numpy.uint32_t, ndim=1]
        The bandwidth expressed in indices. It is up to the caller of this
        function to determine what this corresponds to.
    """
    return find_bandwidths(signal, find_centers(signal, tol), bwcond)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef int find_ordered_zero(numpy.ndarray[numpy.float64_t, ndim=1] arr,
                           int upper, int lower):
    """
    Finds the index of the element in arrLen with value 0.

    Parameters
    ----------
    arr : numpy.ndarray[numpy.float64_t, ndim=1]
        The ordered array.
    upper : int
        The upper limit index of the search.
    lower : int
        The lower limit index of the search.

    Returns
    -------
    int
        The index of the 0-element in arr, or -1 if no 0-element was found.
    """
    if ((arr[upper] > 0 and arr[lower] > 0)
            or (arr[upper] < 0 and arr[lower] < 0)):
        return -1
    cpdef int center = (upper-lower)/2
    if arr[center] < 0:
        return find_ordered_zero(arr, upper, center)
    elif arr[center] > 0:
        return find_ordered_zero(arr, center, lower)
    else:
        return center
