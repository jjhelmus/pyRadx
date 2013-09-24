from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from cython.operator cimport dereference as deref
from cython cimport view

import numpy as np
cimport numpy as np

from _RadxVol_hh cimport RadxVol, RadxSweep, RadxRay, RadxField
# TODO cimport from pxd files for other classes.

cdef extern from "RadxFile.hh":
    cdef cppclass RadxFile:
        RadxFile()
        char* PATH_SEPARATOR
        int readFromPath(string &path, RadxVol &vol)
        string getErrStr()


cdef class PyRadxField:
    """
    Python wrapper around RadxField C++ object 
    """
    
    cdef RadxField *_p

    def __cinit__(self):
        pass

    def __dealloc__(self):
        pass 

    cdef __foo__(self, RadxField* p):
        """ method for loading RadxField object pointer. """
        self._p = p

    def getname(self):
        return self._p.getName()

    def getlongname(self):
        return self._p.getLongName()

    def getstandardname(self):
        return self._p.getStandardName()

    def getnrays(self):
        return self._p.getNRays()

    def getnpoints(self):
        return self._p.getNPoints()

    def getnbytes(self):
        return self._p.getNBytes()

    def getdatatype(self):
        return self._p.getDataType()

    def getdata(self):
        cdef const void *p
        cdef view.array my_array
        p = self._p.getData()
        n_points = self._p.getNPoints()
        data_type = self._p.getDataType()
        if data_type == 1:
            # this is a hack, no care is taken to track references to
            # p, the pointer to the data, also it is cast from const *void
            # to a *void which is a bad idea... but it works for the time 
            # being
            my_array = <short[:n_points]> p
            return np.array(my_array, copy=False)
        else:
            raise NotImplemented
        #cdef np.npy_intp shape[1]
        #shape[0] = <np.npy_intp> 801
        #arr = np.PyArray_SimpleNewFromData(1, shape, np.NPY_INT, p)
        #return arr

cdef class PyRadxRay:
    """
    Python wrapper around RadxRay C++ object 
    """
    
    cdef RadxRay *_p

    def __cinit__(self):
        pass

    def __dealloc__(self):
        pass 

    cdef __foo__(self, RadxRay* p):
        self._p = p

    def getazimuthdeg(self):
        return self._p.getAzimuthDeg()
    
    def getnfields(self):
        return self._p.getNFields()

    def getfield(self, index):
        cdef RadxField *radxfield = self._p.getField(index)
        field = PyRadxField()
        field.__foo__(radxfield)
        return field
    
    def getfields(self):
        cdef RadxField *radxfield
        cdef vector[RadxField *] radxfields

        radxfields = self._p.getFields()
        fields = []
        for radxfield in radxfields:
            field = PyRadxField()
            field.__foo__(radxfield)
            fields.append(field)
        return fields

 

cdef class PyRadxSweep:
    """
    Python wrapper around RadxSweep C++ object 
    """
    
    cdef RadxSweep *_p

    def __cinit__(self):
        pass

    def __dealloc__(self):
        pass 
        #del self._p

    cdef __foo__(self, RadxSweep* p):
        self._p = p

    def info(self):
        a = self._p.getSweepNumber()
        b = self._p.getStartRayIndex()
        c = self._p.getEndRayIndex()
        return a, b, c


cdef class PyRadxVol:
    """
    Python wrapper around RadxVol C++ object 
    """
    
    cdef RadxVol *_p

    def __cinit__(self):
        self._p = new RadxVol()
    
    def __dealloc__(self):
        del self._p

    def getsweepbynumber(self, int sweepnum):
        cdef RadxSweep *radxsweep
        radxsweep = self._p.getSweepByNumber(sweepnum)
        sweep = PyRadxSweep()
        sweep.__foo__(radxsweep)
        return sweep

    def getsweeps(self):
        cdef RadxSweep *radxsweep
        cdef vector[RadxSweep *] radxsweeps
        
        radxsweeps = self._p.getSweeps()
        
        #sweep = PyRadxSweep()
        #sweep.__foo__(radxsweeps[0])
        #return sweep
        
        sweeps = []
        for radxsweep in radxsweeps:
            sweep = PyRadxSweep()
            sweep.__foo__(radxsweep)
            sweeps.append(sweep)
        return sweeps

    def getrays(self):
        cdef RadxRay *radxray
        cdef vector[RadxRay *] radxrays

        radxrays = self._p.getRays()
        rays = []
        for radxray in radxrays:
            ray = PyRadxRay()
            ray.__foo__(radxray)
            rays.append(ray)
        return rays

    def setrayfieldpointer(self):
        self._p.setRayFieldPointers()

    # properties mapped to get/set routines in C++ class
    """ yank 5 lines
    property XXX:
        def __get__(self):
            return self._p.getXXX()
        def __set__(self, val):
            self._p.setXXX(val)
    """
    property debug:
        #def __get__(self):
        #    return self._p._debug
        def __set__(self, bool val):
            self._p.setDebug(val)
    
    property title:
        def __get__(self):
            return self._p.getTitle()
        
        def __set__(self, string val):
            self._p.setTitle(val)
    
    property institution:
        def __get__(self):
            return self._p.getInstitution()
        def __set__(self, const string val):
            self._p.setInstitution(val)
    
    property nsweeps:
        def __get__(self):
            return self._p.getNSweeps()

    property nrays:
        def __get__(self):
            return self._p.getNRays()

    property nfields:
        def __get__(self):
            return self._p.getNFields()

    # TODO add the following properties
    # references
    # source
    # history
    # comment
    # statusxml
    # scanname
    # scanid
    # volumenumber
    # starttime ??? need to wrap time_t
    # endtime   ??? need to wrap time_t
    # targetscanratedegpersec (can be overloaded)
    # ngates
    # ngatesconstant ???
    # maxrangekm
    # packingfromrays ???
    # instrumentname
    # sitename
    # latitude
    # longitude
    # altitudekm
    # sensorhtaglm
    # frequencyhz
    # wavelengthm       also add
    # wavelengthcm      also add
    # radarbeamwidthdegh
    # radarbeamwidthdegv
    # radarantennagaindbh
    # radarantennagaindbv
    # radarrecieverbandwidthmhz
    # islidar
    # lidarconstant
    # lidarpulseenergyj
    # lidarpeakpowerw
    # lidaraperturedaimcm
    # lidarapertureefficiency
    # lidarfieldofviewmrad
    # lidardeamdivergencemrad
    # nsweeps   No setter
    # nrays     No setter
    # nfields   No setter

cdef class PyRadxFile:
    """
    Python wrapper around RadxFile C++ object 
    """
    
    cdef RadxFile *_p
    cdef public PyRadxVol vol
    
    def __cinit__(self):
        self._p = new RadxFile()
    
    def __dealloc__(self):
        del self._p

    def load(self, fname):
        self.vol = PyRadxVol()
        return self._p.readFromPath(fname, deref(self.vol._p))

    def getErrStr(self):
        return self._p.getErrStr()
