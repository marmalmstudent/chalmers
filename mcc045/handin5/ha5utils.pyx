"""
This code is adapted for 80 columns
|------------------------------------------------------------------------------|
"""
import numpy
cimport numpy
cimport cython

# turn off bounds-checking for entire function
@cython.boundscheck(False)
# turn off negative index wrapping for entire function
@cython.wraparound(False)
cpdef numpy.ndarray[numpy.double_t, ndim=3] \
    vecMatDot(numpy.ndarray[numpy.int32_t, ndim=1] dims,
              numpy.ndarray[numpy.double_t, ndim=3] mat1,
              numpy.ndarray[numpy.double_t, ndim=3] mat2):
    """
    Performs the dot product for two matrices whose elements are vectors.
    Axis 0 are treated as the rows, axis 1 are treated as the columns and
    axis 2 are the matrix elements. The dot product is performed on axis 0
    and axis 1 while multiplication with axis 2 will be elementwise.

    Parameters
    ----------
    dims : numpy.ndarray[numpy.int32_t, ndim=1]
        The dimensions of the returned array; [row, col, arrLen].
    mat1 : numpy.ndarray[numpy.double_t, ndim=3]
        The left matrix.
    mat2 : numpy.ndarray[numpy.double_t, ndim=3]
        The right matrix.

    Returns
    -------
    numpy.ndarray[numpy.double_t, ndim=3]
        A matrix whose elements are vectors and is the result of the dot
        product.

    Example
    -------
    >>> import ha5utils
    >>> import numpy as np
    >>> mat=np.arange(0, 12, dtype=np.float64).reshape((2,2,3))
    >>> dim=np.array([mat.shape[0], mat.shape[1], mat.shape[2]], dtype=np.int32)
    >>> mat
    array([[[  0.,   1.,   2.],
            [  3.,   4.,   5.]],

           [[  6.,   7.,   8.],
            [  9.,  10.,  11.]]])
    >>> ha5utils.vecMatDot(dim, mat, mat)
    array([[[  18.,   29.,   44.],
            [  27.,   44.,   65.]],

           [[  54.,   77.,  104.],
            [  99.,  128.,  161.]]])
    """
    cdef numpy.ndarray[numpy.double_t, ndim=3] out
    cdef int i,j, k, l
    out = numpy.zeros((dims[0], dims[1], dims[2]),
		      dtype=numpy.double)
    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[0]):
                for l in range(dims[2]):
                    out[i,j, l] += mat1[i, k, l]*mat2[k, j, l]
    return out
