#! /usr/bin/env python

import radx

rfile = radx.PyRadxFile()
rfile.load('20080604_002217_SPOLRVP8_v036_SUR.uf')

# examine the RadxVol
vol = rfile.vol
print ""
print "Volume parameters"
print "-----------------"
print "Title:", vol.title
print "Institution:", vol.institution
print "NSweeps:", vol.nsweeps
print "NRays:", vol.nrays
print "NFields:", vol.nfields

# examine the first RadxSweep
sweep = vol.getsweepbynumber(0)
sweep_num, start_idx, end_idx = sweep.info()
print ""
print "Sweep parameters"
print "----------------"
print "SweepNumber:", sweep_num
print "StartRayIndex:", start_idx
print "EndRayIndex:", end_idx

# examine the first ray
ray = vol.getrays()[0]
print ""
print "Ray parameters"
print "--------------"
print "AzimuthDeg:", ray.getazimuthdeg()
print "NFields:", ray.getnfields()

# examine the first field
field = ray.getfields()[0]
print ""
print "Field parameters"
print "----------------"
print "Name:", field.getname()
print "LongName:", field.getlongname()
print "StandardName:", field.getstandardname()
print "NRays:", field.getnrays()
print "NPoints:", field.getnpoints()
print "NBytes:", field.getnbytes()
print "DataType:", field.getdatatype()

# examine the data
data = field.getdata()
print ""
print "Data"
print "----"
for i in range(10):
    print "Point %i: %i" % (i, data[i])


