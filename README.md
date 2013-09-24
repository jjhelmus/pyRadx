pyRadx is the beginning of a Python wrapper around 
[Radx](http://www.ral.ucar.edu/projects/titan/docs/radial_formats/radx.html)
using [Cython](http://cython.org/). Currently it is protype/testing code to
check the feasibility of such a wrapper.


Building pyRadx
---------------

* Download the RadX library and unpack.

* Download, build and install RadX with -fPIC flag. This can be accomplished
  using ./configure CXXFLAGS="-fPIC", make and [sudo] make install.

* Edit the `include_dirs` and `library_dirs` lines in setup.py to point
  to the location of the RadX libraries (.a files) and include files (.hh).

* Build the Cython wrapper around RadX with 'python setup.py build\_ext -i'.

Trying it out
-------------

The `example.py` script shows a bit of what pyRadX can do.  It reads an 
[example uf](http://www.ral.ucar.edu/projects/titan/example_data/uf/20080604/)
file and prints out some parameters.  The expected output is below:

    $ ./example.py 

    Volume parameters
    -----------------
    Title: TIMREX
    Institution: 
    NSweeps: 9
    NRays: 4343
    NFields: 0

    Sweep parameters
    ----------------
    SweepNumber: 0
    StartRayIndex: 0
    EndRayIndex: 482

    Ray parameters
    --------------
    AzimuthDeg: 121.5
    NFields: 2

    Field parameters
    ----------------
    Name: DZ
    LongName: reflectivity
    StandardName: 
    NRays: 1
    NPoints: 996
    NBytes: 1992
    DataType: 1

    Data
    ----
    Point 0: -1694
    Point 1: 8997
    Point 2: 20816
    Point 3: 20628
    Point 4: 20722
    Point 5: 2054
    Point 6: 3334
    Point 7: 3947
    Point 8: 3108
    Point 9: 16157
