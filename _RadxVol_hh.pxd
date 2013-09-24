# Cython pxd version of RadxVol.hh
# Author: Jonathan J. Helmus
"""
/*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
** Copyright UCAR (c) 1992 - 1999
** University Corporation for Atmospheric Research(UCAR)
** National Center for Atmospheric Research(NCAR)
** Research Applications Program(RAP)
** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
** All rights reserved. Licenced use only.
** Do not copy or distribute without authorization
** 1999/03/14 14:18:54
*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*/
/////////////////////////////////////////////////////////////
// RadxVol.hh
//
// RadxVol object
//
// NetCDF data for radar radial data in CF-compliant file
//
// Mike Dixon, RAP, NCAR
// P.O.Box 3000, Boulder, CO, 80307-3000, USA
//
// Dec 2009
//
///////////////////////////////////////////////////////////////
"""

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
#include <string>
#include <vector>
#include <Radx/Radx.hh>
#include <Radx/RadxPlatform.hh>
#include <Radx/RadxField.hh>
#class RadxSweep;
#class RadxRay;
#class RadxRcalib;
#class RadxFile;
#class RadxCfactors;
#using namespace std;

# XXX this should be in the _RadxSweep_hh.pxd or similar file...
cdef extern from "RadxSweep.hh":
    cdef cppclass RadxSweep:
        RadxSweep()
        int getSweepNumber()
        int getStartRayIndex()
        int getEndRayIndex()

# XXX this should be in the _RadxRay_hh.pxd or similar file...
cdef extern from "RadxRay.hh":
    cdef cppclass RadxRay:
        RadxRay()
        double getAzimuthDeg()
        size_t getNFields()
        RadxField * getField(int index)
        vector[RadxField *] getFields()

# XXX this should be in the _RadxRay_hh.pxd or similar file...
cdef extern from "RadxField.hh":
    cdef cppclass RadxField:
        RadxField()
        string getName()
        string getLongName()
        string getStandardName()
        size_t getNRays()
        size_t getNPoints()
        size_t getNBytes()
        int getDataType()
        const void * getData()
"""
//////////////////////////////////////////////////////////////////////
/// CLASS FOR STORING VOLUME (and SWEEP) data
/// 
/// This is the primary class for storing data from RADAR and LIDAR
/// volumes.
///
/// A volume is made up of a number of sweeps (RadxSweep), each of
/// which contain a number of rays (RadxRay).
///
/// A volume may contain one or more sweeps.
///
/// The information held by the volume object is references in a
/// number of ways. There are 3 primary vectors on a volume, and
/// all are inter-related:
///
/// \code
///  (a) rays
///  (b) sweeps
///  (c) fields 
/// \endcode
///
/// The RAYS are stored as a vector of RadxRay objects. Each RadxRay
/// owns a vector of field objects - see RadxRay for more information
/// on this.
///
/// A SWEEP is formed from a sequence of RAYS. RadxSweep is a
/// light-weight object which merely keeps information about how the
/// rays are grouped into sweeps. The sweep objects do not hold data,
/// they just keep track of the ray indices, and sweep properties.
///
/// As each ray is created, it manages the field data for that ray.
/// When this process is complete, the volume is a collection of
/// rays, which autonomously manage their own data and geometry.
///
/// A call to loadVolFieldsFromRays() will create fields on the
/// volume, copy the data from the rays into contiguous arrays on
/// these new fields, and then change the ray fields so that they
/// point to the volume fields rather than manage their own data.
/// After this step, the volume fields hold contiguous arrays for the
/// field data. These volume-owned fields manage the field data
/// memory, and the rays now just hold pointers into these main
/// fields.
///
/// NOTE ON CONVERTING DATA TYPES IN FIELDS:
///
/// If you convert the data type for individual fields, rather than
/// calling the methods on the RadxVol class, you must call
/// setRayFieldPointers() before using the field data in the RadxRay
/// objects. The call will ensure that the pointers in the field
/// objects on the rays points correctly to the data in the field
/// objects on the volume.
"""
cdef extern from "RadxVol.hh":
    cdef cppclass RadxVol:

        #public RadxRangeGeom, public RadxPacking {

        #public:

        #  /// Constructor
  
        RadxVol()

#  /// Copy constructor
  
#  RadxVol(const RadxVol &rhs);

#  /// Copy constructor, one sweep only
  
#  RadxVol(const RadxVol &rhs, int sweepNum);

#  /// Destructor
  
#  virtual ~RadxVol();

#  /// Assignment
  
#  RadxVol& operator=(const RadxVol &rhs);
  
#  /// Copy one sweep only into existing object,
#  /// clearing other sweeps as needed
  
#  void copy(const RadxVol &rhs, int sweepNum);

#  /// copy following meta data only:
#  ///   main members, calibs, sweeps as in file
#  ///
#  /// does not copy sweeps, rays and fields
  
#  void copyMeta(const RadxVol &rhs);
  
#  //////////////////////////////////////////////////////////////////
#  /// \name Set methods - except for platform parameters
#  //@{
#
#  /// Set debugging on/off. Off by default.

        void setDebug(bool val)

#  /// Set the volume title, if available. Use this for the project name.

        void setTitle(string &val)

#  /// Set the institution responsible for gathering the data.

        void setInstitution(const string &val)

#  /// Set a references string. Use this for the flight number, if
#  /// appropriate.

        void setReferences(const string &val)
  
#  /// Set source. Use this for the name of the facility generating the
#  /// data.

        void setSource(const string &val)

#  /// Set history.
#  ///
#  /// This string should be appended to as different
#  /// operations are performed on the data set, to provide a history
#  /// of those operations.

        void setHistory(const string &val)

#  /// Set comment as applicable.

        void setComment(const string &val)

#  /// Set status XML as applicable. For general-purpose status information.

        void setStatusXml(const string &val)

#  /// Set the scan strategy name, if available.

        inline void setScanName(const string &val)

#  /// Set the scan strategy id (VCP), if available.

        inline void setScanId(int val)

#  /// Set the volume number.
#  ///
#  /// This increments with every volume, and may wrap.

        inline void setVolumeNumber(int val)
  
#  /// Set the volume start time.
#  ///
#  /// The time is split into two parts, (a) the seconds UTC since Jan
#  /// 1 1970, and the fractional seconds converted to nanosecons.

#        inline void setStartTime(time_t secs, double nanoSecs)

#  /// Set the volume end time.
#  ///
#  /// The time is split into two parts, (a) the seconds UTC since Jan
#  /// 1 1970, and the fractional seconds converted to nanosecons.

#        inline void setEndTime(time_t secs, double nanoSecs)

#  /// Set the correction factors on the volume, if applicable.
#  ///
#  /// If not set, all corrections will be assumed to be 0.0.

#  void setCfactors(RadxCfactors &cfac);
  
#  /// set the target scan rate for all rays in a volume

        void setTargetScanRateDegPerSec(double rate)
  
#  /// set the target scan rate for all rays in a sweep

        void setTargetScanRateDegPerSec(int sweepNum, double rate)

#  /// Compute the max number of gates, by searching through the rays.
#  /// Also determines if number of gates vary by ray.
#  /// After this call, use the following to get the results:
#  ///   size_t getMaxNGates() - see RadxPacking.
#  ///   bool nGatesVary() - see RadxPacking

#  void computeMaxNGates() const;

#  /// Set the number of gates.
#  ///
#  /// If more gates are needed, extend the field data out to a set number of
#  /// gates. The data for extra gates are set to missing values.
#  ///
#  /// If fewer gates are needed, the data is truncated.
  
        void setNGates(size_t nGates)

#  /// Set to constant number of gates per ray.
#  /// 
#  /// First we determine the max number of gates, and also check
#  /// for a variable number of gates. If the number of gates does
#  /// vary, the shorter rays are padded out with missing data.
  
        void setNGatesConstant()

#  /////////////////////////////////////////////////////////////////////////
#  /// Set the maximum range.
#  /// Removes excess gates as needed.
#  /// Does nothing if the current max range is less than that specified.
  
        void setMaxRangeKm(double maxRangeKm)
  
#  /// set the packing from the rays

        void setPackingFromRays()

#  //@}
#
#  //////////////////////////////////////////////////////////////////
#  /// \name Set methods specific to platform
#  //@{
#
#  /// Set the platform, using a filled-out object

#  void setPlatform(const RadxPlatform &val) { _platform = val; }

#  /// Set the instrument name, if available.

        inline void setInstrumentName(const string &val)

#  /// Set the site name, if available.

        inline void setSiteName(const string &val)

#  /// Set the instrument type. Default is RADAR.

#  inline void setInstrumentType(Radx::InstrumentType_t val) {
#    _platform.setInstrumentType(val);
#  }

#  /// Set the platform type. Default is FIXED.

#  inline void setPlatformType(Radx::PlatformType_t val) {
#    _platform.setPlatformType(val);
#  }

#  /// Set the primary rotation axis. Default is AXIS_Z.

#  inline void setPrimaryAxis(Radx::PrimaryAxis_t val) {
#    _platform.setPrimaryAxis(val);
#  }

#  /// Set the latitude of the platform in degrees.
#  ///
#  /// Used for non-mobile platforms.

        void setLatitudeDeg(double val)

#  /// Set the longitude of the platform in degrees.
#  ///
#  /// Used for non-mobile platforms.

        void setLongitudeDeg(double val)

#  /// Set the altitude of the platform in km.
#  ///
#  /// Used for non-mobile platforms.

        void setAltitudeKm(double val)

#  /// Set the sensor ht above the surface

        void setSensorHtAglM(double val)

#  /// Set up the list of frequencies, adding them one at a time.
#  /// Normally there is only one frequency, but multiple are supported.

#  /// The set methods clear the list first, and then add the value
#  /// The add methods do not clear the list first
  
        void setFrequencyHz(double val)
        void setWavelengthM(double val)
        void setWavelengthCm(double val)
  
        void addFrequencyHz(double val)
        void addWavelengthM(double val)
        void addWavelengthCm(double val)
  
#  /// Set the RADAR beam width, horizontal, in degrees.

        void setRadarBeamWidthDegH(double val)

#  /// Set the RADAR beam width, vertical, in degrees.

        void setRadarBeamWidthDegV(double val)

#  /// Set the RADAR antenna gain, horizontal, in dB.

        void setRadarAntennaGainDbH(double val)

#  /// Set the RADAR antenna gain, vertical, in dB.

        void setRadarAntennaGainDbV(double val)

#  /// Set the RADAR receiver bandwidth, in MHz.

        void setRadarReceiverBandwidthMhz(double val)

#  /// Set the instrument type to LIDAR. Default is RADAR.

        void setIsLidar(bool val)

#  /// Set the LIDAR constant.

        void setLidarConstant(double val)

#  /// Set the LIDAR pulse energy, in Joules.

        void setLidarPulseEnergyJ(double val)

#  /// Set the LIDAR peak power, in Watts.

        void setLidarPeakPowerW(double val)

#  /// Set the LIDAR aperture diameter, in cm.

        void setLidarApertureDiamCm(double val)

#  /// Set the LIDAR aperture efficiency, in percent.

        void setLidarApertureEfficiency(double val)

#  /// Set the LIDAR field of view, in milli-radians.

        void setLidarFieldOfViewMrad(double val)
  
#  /// Set the LIDAR beam divergence, in milli-radians.

        void setLidarBeamDivergenceMrad(double val)

#  //@}
#
#  //////////////////////////////////////////////////////////////////
#  /// \name Add methods:
#  //@{
#
#  /// Add a ray to the volume.
#  ///
#  /// The ray must be created with new() before calling this method.
#  /// RadxVol takes responsibility for freeing the object.
  
#  void addRay(RadxRay *ray);

#  /// Add a sweep object.
#  ///
#  /// The sweep must be created with new() before calling this method.
#  /// RadxVol takes responsibility for freeing the sweep object.
  
#  void addSweep(RadxSweep *sweep);

#  /// Add sweep to 'as-in-file' vector
#  /// These are the sweeps as originally in the file before selected
#  /// sweeps were requested.
#  /// A copy of the object is made, and is managed by this object.

#  void addSweepAsInFile(RadxSweep *sweep);
  
#  /// Add a calibration to the volume.
#  ///
#  /// The calibration must be created with new() before calling this method.
#  /// RadxVol takes responsibility for freeing the calibration object.
  
#  void addCalib(RadxRcalib *calib);
  
#  /// Add a field to the volume.
#  ///
#  /// The field must be created with new() before calling this method.
#  /// RadxVol takes responsibility for freeing the field object.

#  void addField(RadxField *field);

#  //@}
#
#  //////////////////////////////////////////////////////////////////
#  /// \name Change data representation in volume:
#  //@{
#
#  /// After you you have added rays with field data local to each ray,
#  /// you may wish to convert those to contiguous fields and
#  /// initialize other information on the volume before using it.
#  ///
#  /// The method loads up contiguous fields in the volume from data in
#  /// the rays. The fields in the rays are then set to point to the
#  /// contiguous fields.
#  ///
#  /// See also loadRaysFromFields().

#  void loadFieldsFromRays();

#  /// Load up the ray fields from the contiguous fields in the volume.
#  /// This is the inverse of loadFieldsFromRays()

#  void loadRaysFromFields();
  
#  /// Set field data pointers in the rays to point into the
#  /// main contiguous fields on the volume.
#  ///
#  /// If the rays that have been added to the volume do not hold the
#  /// data themselves, call this method to set the pointers in the ray
#  /// fields so that they refer to the main fields on the volume.
  
        void setRayFieldPointers()

#  /// Make a copy of the field with the specified name.
#  /// This forms a contiguous field from the ray data.
#  /// Returns a pointer to the field on success, NULL on failure.

#  RadxField *copyField(const string &fieldName) const;

#  /// Rename a field
#  /// returns 0 on success, -1 if field does not exist in any ray

#  int renameField(const string &oldName, const string &newName);
  
#  /// Load volume information from the rays.
#  ///
#  /// This sets the volume number and the start and end times.
  
#  void loadVolumeInfoFromRays();
  
#  /// Load the sweep information from the rays.
#  ///
#  /// This loops through all of the rays, and determines the sweep
#  /// information from them. The resulting information is stored
#  /// in the sweeps array on the volume.
#  ///
#  /// Also sets the start/end of sweep/volume flags
  
#  void loadSweepInfoFromRays();
  
#  /// Set the modes on the rays, by copying from the sweeps.
#  ///
#  /// If the sweep mode, prf mode etc have been set correctly on the 
#  /// sweep objects, you can call this method to copy the modes over to
#  /// the ray objects.
  
#  void loadModesFromSweepsToRays();

#  /// Set the fixed angle on the rays, by copying from the sweeps.
#  ///
#  /// If the fixed angles have been set correctly on the sweep
#  /// objects, you can call this method to copy the modes over to the
#  /// ray objects.
  
#  void loadFixedAnglesFromSweepsToRays();

#  /// load the calbration index on the rays, using the pulse width to
#  /// determine which calibration is relevant.
#  ///
#  /// This method checks the pulse witdth on each ray, and compares
#  /// these with the pulse width for each calibration object. It then
#  /// sets the calibration index to point to the calibration data with
#  /// the pulse width closest to that on the ray.
  
#  void loadCalibIndexOnRays();

#  /// Constrain the data by specifying fixedAngle limits.
#  ///
#  /// This operation will remove unwanted rays from the data set,
#  /// remap the field arrays for the remaining rays and set the field
#  /// pointers in the rays to the remapped fields on the volume.
  
#  void constrainByFixedAngle(double minFixedAngleDeg,
#                             double maxFixedAngleDeg);

#  /// Constrain the data by specifying sweep number limits.
#  ///
#  /// This operation will remove unwanted rays from the data set,
#  /// remap the field arrays for the remaining rays and set the field
#  /// pointers in the rays to the remapped fields on the volume.

#  void constrainBySweepNum(int minSweepNum,
#                           int maxSweepNum);

#  /// remove rays with all missing data
  
#  void removeRaysWithDataAllMissing();

#  /////////////////////////////////////////////////////////////////
#  /// Ensure ray times are monotonically increasing by
#  /// interpolating the times if there are duplicates

#  void interpRayTimes();

#  /////////////////////////////////////////////////////////////////
#  /// Remap all fields and rays onto the specified geometry.
#  ///
#  /// This leaves the memory managed by the rays.
#  /// Call loadFieldsFromRays() if you need the field data
#  /// to be managed by the volume.
#  ///
#  /// If interp is true, uses interpolation if appropriate.
#  /// Otherwise uses nearest neighbor.
  
#  virtual void remapRangeGeom(double startRangeKm,
#                              double gateSpacingKm,
#                              bool interp = false);

#  /// Remap data in all rays to the predominant range geometry.
#  ///
#  /// A search is made through the rays, to identify which is the
#  /// predominant range geometry.  All rays which do not match this
#  /// are then remapped to this predominant geometry.
#  ///
#  /// This leaves the memory managed by the rays.
#  /// Call loadFieldsFromRays() if you need the field data
#  /// to be managed by the volume.
  
#  void remapToPredomGeom();
  
#  /// Remove rays which do not match the predominant range geometry.
#  ///
#  /// A search is made through the rays, to identify which is the
#  /// predominant range geometry.  All rays which do not match this
#  /// are then removed from the volume.
#  ///
#  /// This leaves the memory managed by the rays.
#  /// Call loadFieldsFromRays() if you need the field data
#  /// to be managed by the volume.
  
#  void filterOnPredomGeom();

#  /// Copy the range geom from the fields to the rays, provided
#  /// the fields have consistent in geometry
  
#  void copyRangeGeomFromFieldsToRays();

#  ////////////////////////////////////////////////////
#  /// Copy the range geom to the fields from the rays
  
#  void copyRangeGeomFromRaysToFields();
  
#  /// filter based on ray vectors
#  /// Keep the good rays, remove the bad rays
  
#  void removeBadRays(vector<RadxRay *> &goodRays,
#                     vector<RadxRay *> &badRays);
  
#  /// Remove rays with the antenna transition flag set.
  
#  void removeTransitionRays();

#  /// Remove rays with transitions, with the specified margin.
#  ///
#  /// Sometimes the transition flag is turned on too early in
#  /// a transition, on not turned off quickly enough after a transition.
#  /// If you set this to a number greater than 0, that number of rays
#  /// will be included at each end of the transition, i.e. the
#  /// transition will effectively be shorter at each end by this
#  /// number of rays
  
#  void removeTransitionRays(int nRaysMargin);

#  /// Trim surveillance sweeps to 360 deg
#  ///
#  /// Remove extra rays in each surveillance sweep
  
#  void trimSurveillanceSweepsTo360Deg();

#  /// Reorder the sweeps into ascending angle order
#  ///
#  /// If the sweeps are reordered, this means that the rays times
#  /// will no longer be monotonically increasing

#  void reorderSweepsAscendingAngle();

#  /// Reorder the sweeps as in file into ascending angle order
  
#  void reorderSweepsAsInFileAscendingAngle();

#  /// Apply an azimuth offset to all rays in the volume
#  /// This applies to the rays currently in the volume, not to
#  /// any future reads
  
#  void applyAzimuthOffset(double offset);

#  /// Apply an elevation offset to all rays in the volume
#  /// This applies to the rays currently in the volume, not to
#  /// any future reads
  
#  void applyElevationOffset(double offset);

#  /// Set the fixed angle for a sweep.
#  /// Also sets the fixed angle for the rays in the sweep
  
#  void setFixedAngleDeg(int sweepNum, double fixedAngle);

#  /// combine rays from sweeps with common fixed angle and
#  /// gate geometry, but with different fields
  
#  void combineSweepsAtSameFixedAngleAndGeom(bool keepLongRange = false);

#  /// Make fields uniform in the volume.
#  /// This ensures that all rays in the volume have the same fields
#  /// and that they are in the same order in each ray.
#  /// If fields a missing from a ray, a suitable field is added
#  /// containing missing data.

#  void makeFieldsUniform();
  
#  /// Make fields uniform for each sweep.
#  /// This ensures that all rays in a sweep have the same fields.
#  /// and that they are in the same order in each ray.
#  /// If fields a missing from a ray, a suitable field is added
#  /// containing missing data.
  
#  void makeFieldsUniformPerSweep();
  
#  /// Reorder the fields, by name, removing any extra fields.
  
#  void reorderFieldsByName(const vector<string> &names);
  
#  /// apply the georeference corrections for moving platforms
#  /// to compute the earth-relative azimuth and elevation

#  void applyGeorefs();

#  //@}
#
#  //////////////////////////////////////////////////////////////////
#  /// \name Convert the data type for all fields
#  //@{

#  /// Convert all fields to 64-bit floats.
  
#  void convertToFl64();

#  /// Convert all fields to 32-bit floats.
  
#  void convertToFl32();

#  /// Convert all fields to scaled 32-bit signed integers,
#  /// using the specified scale and offset.
#  /// See RadxField for correct use of the scale and offset.
  
#  void convertToSi32(double scale, double offset);

#  /// Convert all fields to scaled 32-bit signed integers,
#  /// dynamically computing the scale and offset.
  
#  void convertToSi32();

#  /// Convert all fields to scaled 16-bit signed integers,
#  /// using the specified scale and offset.
#  /// See RadxField for correct use of the scale and offset.

#  void convertToSi16(double scale, double offset);

#  /// Convert all fields to scaled 16-bit signed integers,
#  /// dynamically computing the scale and offset.
  
#  void convertToSi16();

#  /// Convert all fields to scaled 8-bit signed integers,
#  /// using the specified scale and offset.
#  /// See RadxField for correct use of the scale and offset.
  
#  void convertToSi08(double scale, double offset);
  
#  /// Convert all fields to scaled 8-bit signed integers,
#  /// dynamically computing the scale and offset.
  
#  void convertToSi08();

#  /// Convert all fields to the specified data type.
  
#  void convertToType(Radx::DataType_t targetType);

#  /// Converts field type, and optionally changes the
#  /// names.
#  ///
#  /// If the data type is an integer type, dynamic scaling
#  /// is used - i.e. the min and max value is computed and
#  /// the scale and offset are set to values which maximize the
#  /// dynamic range.
#  ///
#  /// If targetType is Radx::ASIS, no conversion is performed.
#  ///
#  /// If a string is empty, the value on the field will
#  /// be left unchanged.
  
#  void convertField(const string &name,
#                    Radx::DataType_t type,
#                    const string &newName,
#                    const string &units,
#                    const string &standardName,
#                    const string &longName);

 # /// Converts field type, and optionally changes the
 # /// names.
 # ///
 # /// For integer types, the specified scale and offset
 # /// are used.
 # ///
 # /// If targetType is Radx::ASIS, no conversion is performed.
 # ///
 # /// If a string is empty, the value on the field will
 # /// be left unchanged.
  
#  void convertField(const string &name,
#                    Radx::DataType_t type,
#                    double scale,
#                    double offset,
#                    const string &newName,
#                    const string &units,
#                    const string &standardName,
#                    const string &longName);

#  //@}
#
#  //////////////////////////////////////////////////////////////////
#  /// \name Get methods - except for platform parameters
#  //@{

#  /// Get title. May be used for project name.

        string getTitle()

#  /// Get institution responsible for gathering the data.

        const string getInstitution()

#  /// Get references. May be used for flight number.

        inline const string &getReferences()

#  /// Get source. Should be used for generating facility.

        inline const string &getSource()

#  /// Get history.

        inline const string &getHistory()

#  /// Get comments.

        inline const string &getComment()

#  /// Get status XML. For general-purpose status information.

        inline const string &getStatusXml()

#  /// Get scan name.

        inline const string &getScanName()

#  /// Get scan id (VCP).

        inline int getScanId()

#  /// Get volume number.

        inline int getVolumeNumber()

#  /// Get start time in seconds.
#  /// Combine with getNanoSecs() for high-precision time.
#  ///
#  /// \code
#  /// double time = (time in secs) + (nano seconds) / 1.0e9.
#  /// \endcode

#  inline double getStartNanoSecs() const { return _startNanoSecs; } 

#  /// Get nano-seconds for start time.

#  inline time_t getStartTimeSecs() const { return _startTimeSecs; }

#  /// Get end time in seconds.
#  /// Combine with getNanoSecs() for high-precision time.
#  ///
#  /// \code
#  /// double time = (time in secs) + (nano seconds) / 1.0e9.
#  /// \endcode

#  inline time_t getEndTimeSecs() const { return _endTimeSecs; }

#  /// Get nano-seconds for end time.

#  inline double getEndNanoSecs() const { return _endNanoSecs; }

#  /// get flag to indicate that the ray times are in order -
#  /// i.e. they are increasing

#  inline bool getRayTimesIncrease() const { return _rayTimesIncrease; }

#  /// Get number of sweeps in volume.

        size_t getNSweeps()

#  /// Get number of rays in volume.
  
        inline size_t getNRays()

#  /// Get number of fields in volume.

        inline size_t getNFields()
  
#  /// Get the list of unique field names, compiled by
#  /// searching through all rays.
#  ///
#  /// The order of the field names found is preserved

#  vector<string> getUniqueFieldNameList() const;

#  /// Get field, based on index.
#  ///
#  /// Returns NULL on failure.

#  inline RadxField *getField(int index) const {
#    if (index < (int) _fields.size()) {
#      return _fields[index];
#    } else {
#      return NULL;
#    }
#  }

#  /// Get field on the volume, based on name.
#  /// Returns NULL on failure.

#  RadxField *getField(const string &name) const {
#    for (size_t ii = 0; ii < _fields.size(); ii++) {
#      if (_fields[ii]->getName() == name) {
#        return _fields[ii];
#      }
#    }
#    return NULL;
#  }

#  /// Get a field from a ray, given the name.
#  /// Find the first available field on a suitable ray.
#  /// Returns field pointer on success, NULL on failure.
  
#  const RadxField *getFieldFromRay(const string &name) const;

#  /// Get vector of fields.

#  inline const vector<RadxField *> &getFields() const { return _fields; }
#  inline vector<RadxField *> &getFields() { return _fields; }

#  /// Get sweep by sweep number (not the index).
#  /// 
#  /// Returns NULL on failure.

        const RadxSweep *getSweepByNumber(int sweepNum)
#  RadxSweep *getSweepByNumber(int sweepNum);
  
#  /// Get vector of sweeps.

        const vector[RadxSweep *] &getSweeps()
#  const vector<RadxSweep *> &getSweeps() const { return _sweeps; }
#  vector<RadxSweep *> &getSweeps() { return _sweeps; }

#  /// check if all rays in a sweep are in an antenna transition

#  bool checkAllSweepRaysInTransition(const RadxSweep &sweep) const;
#  bool checkAllSweepRaysInTransition(int sweepNum) const;
 
#  /// Get vector of sweeps as they appeared in original file.

#  inline const vector<RadxSweep *> &getSweepsAsInFile() const {
#    return _sweepsAsInFile;
#  }
  
#  /// Get vector of rays.

        inline const vector[RadxRay *] &getRays()
#  inline const vector<RadxRay *> &getRays() const { return _rays; }
#  inline vector<RadxRay *> &getRays() { return _rays; }
  
#  /// Get vector of radar calibrations.

#  inline size_t getNRcalibs() const { return _rcalibs.size(); }
#  inline const vector<RadxRcalib *> &getRcalibs() const { return _rcalibs; }
#  inline vector<RadxRcalib *> &getRcalibs() { return _rcalibs; }
  
#  /// Get pointer to correction factors.
#  ///
#  /// Returns NULL if there are no corrections available.
  
#  inline const RadxCfactors *getCfactors() const { return _cfactors; }
#  inline RadxCfactors *getCfactors() { return _cfactors; }
  
#  /// Determine whether rays are indexed in angle, and what
#  /// the predominant angular resolution is.
#  //
#  /// Returns true if all rays are indexed, false otherwise.
  
#  bool checkForIndexedRays() const;
  
#  /// check whether volume is predominantly in RHI mode
#  ///
#  /// Returns true if RHI, false otherwise

#  bool checkIsRhi() const;

#  //@}
#
#  //////////////////////////////////////////////////////////////////
#  /// \name Get methods for platform parameters
#  //@{

#  /// Get platform object

#  inline const RadxPlatform &getPlatform() const {
#    return _platform;
#  }

#  /// Get instrument name.

        inline const string &getInstrumentName()

#  /// Get site name.

        inline const string &getSiteName()

#  /// Get instrument type

#  inline Radx::InstrumentType_t getInstrumentType() const {
#    return _platform.getInstrumentType();
#  }

#  /// Get platform type

#  inline Radx::PlatformType_t getPlatformType() const {
#    return _platform.getPlatformType();
#  }

#  /// Get primary axis

#  inline Radx::PrimaryAxis_t getPrimaryAxis() const {
#    return _platform.getPrimaryAxis();
#  }

#  /// Get latitude in degrees. Applies to FIXED platform.

        inline double getLatitudeDeg()

#  /// Get longitude in degrees. Applies to FIXED platform.

        inline double getLongitudeDeg()

#  /// Get altitude in km. Applies to FIXED platform.

        inline double getAltitudeKm()

#  /// Get the sensor ht above the surface in meters

        inline double getSensorHtAglM()

#  /// Get vector of frequencies - normally only 1 entry

#  inline const vector<double> &getFrequencyHz() const {
#    return _platform.getFrequencyHz();
#  }
#  double getWavelengthM() const;
#  double getWavelengthCm() const;

#  /// For RADAR, get horizontal beam width, in degrees.

        inline double getRadarBeamWidthDegH()

#  /// For RADAR, get vertical beam width, in degrees.

        inline double getRadarBeamWidthDegV()

#  /// For RADAR, get horizontal antenna gain, in dB.

        inline double getRadarAntennaGainDbH()

#  /// For RADAR, get horizontal antrenna gain, in dB.

        inline double getRadarAntennaGainDbV()

#  /// For RADAR, get receiver band width, in Mhz.

        inline double getRadarReceiverBandwidthMhz()

#  /// For LIDAR, get lidar constant.

        inline double getLidarConstant()

#  /// For LIDAR, get pulse energy, in Joules.

        inline double getLidarPulseEnergyJ()

#  /// For LIDAR, get peak power, in watts.

        inline double getLidarPeakPowerW()

#  /// For LIDAR, get aperture diameter, in cm.

        inline double getLidarApertureDiamCm()

#  /// For LIDAR, get aperture efficiency, in percent.

        inline double getLidarApertureEfficiency()

#  /// For LIDAR, get field of view, in milli-radians.

        inline double getLidarFieldOfViewMrad()

#  /// For LIDAR, get beam divergence, in milli-radians.

        inline double getLidarBeamDivergenceMrad()

#  //@}
#
#  //////////////////////////////////////////////////////////////////
#  /// \name Path used in read/write
#  //@{

#  void setPathInUse(const string &val) const { _pathInUse = val; }

#  /// Get the full file path actually used for reading or writing.
  
#  const string &getPathInUse() const { return _pathInUse; }

#  //@}
#
#  //////////////////////////////////////////////////////////////////
#  /// \name Clearing data:
#  //@{

#  /// Clear all of the data in the object.
  
#  void clear();
  
#  /// Remove all ray on object.

#  void clearRays();

#  /// Clear sweep information on object.

#  void clearSweeps();

#  /// Clear sweep information as they appear in original file.

#  void clearSweepsAsInFile();

#  /// Clear radar calibration info on object.

#  void clearRcalibs();

#  /// Clear the field vector, deleting the RadxField objects first.
  
#  void clearFields();

#  /// Clear the correction factors.
  
#  void clearCfactors();
  
#  /// Clear fields from rays, but retain the rays with their metadata.
  
#  void clearRayFields();
  
#  /// Clear the frequency list.
  
#  void clearFrequency();

#  //@}

#  //////////////////////////////////////////////////////////////////
#  /// \name Printing:
#  //@{
#
#  /// Print metadata on volume object.
  
#  void print(ostream &out) const;
  
#  /// Print ray metadata

#  void printWithRayMetaData(ostream &out) const;

#  /// Print summary of each ray

#  void printRaySummary(ostream &out) const;

#  /// Print full metadata, and actual field data.

#  void printWithFieldData(ostream &out) const;

#  //@}
#
#protected:
#  
#private:

#  // debug state

        bool _debug

#  // class for keeping track of the geometry of the rays and
#  // remapping data onto a common geometry

#  class RayGeom {
#  public:
#    RayGeom();
#    RayGeom(double start_range,
#            double gate_spacing);
#    void print(ostream &out) const;
#    double startRange;
#    double gateSpacing;
#    int rayCount;
#  };

#  // class for combining sweeps with same fixed angle but different fields

#  class Combo {
#  public:
#    size_t target;
#    vector<size_t> sources; 
#    Combo() {
#      target = 0;
#    }
#    Combo(size_t index) {
#      target = index;
#    }
#  };

#  // meta strings

#  string _title;
#  string _institution;
#  string _references;
#  string _source;
#  string _history;
#  string _comment;
#  string _statusXml;

#  // scan details

#  string _scanName;
#  int _scanId; // VCP

#  // platform parameters

#  RadxPlatform _platform;

#  // volume number

#  int _volNum;
  
#  // times

#  time_t _startTimeSecs;
#  time_t _endTimeSecs;
#  double _startNanoSecs;
#  double _endNanoSecs;

#  // flag to indicate that the ray times are in order -
#  // i.e. they are increasing

#  bool _rayTimesIncrease;

#  // transitions

#  vector<bool> _transitionFlags;
  
#  // path in use - for reading/writing
  
#  mutable string _pathInUse; ///< path in use

#  // sweeps
  
#  vector<RadxSweep *> _sweeps;

#  // sweeps as they were in the file
  
#  vector<RadxSweep *> _sweepsAsInFile;

#  // rays

#  vector<RadxRay *> _rays;

#  // calibrations

#  vector<RadxRcalib *> _rcalibs;

#  // fields
  
#  vector<RadxField *> _fields;

#  // correction factors

#  RadxCfactors *_cfactors;

#  // searching for angle match between sweeps

#  static const int _searchAngleN = 36000;
#  static const double _searchAngleRes;
#  int _searchMaxWidth;
#  vector<const RadxRay *> _searchRays;
  
#  // private methods
  
#  void _init();
#  RadxVol & _copy(const RadxVol &rhs);

#  RayGeom _getPredomGeom() const;
#  void _constrainBySweepIndex(vector<int> &sweepIndexes);
#  void _checkForIndexedRays(RadxSweep &sweep) const;
#  double _computeRoundedAngleRes(double res) const;
#  void _findTransitions(int nRaysMargin);
#  void _checkRayTimesIncrease();
#  void _augmentSweepFields(size_t target, size_t source);

#  int _setupAngleSearch(size_t sweepNum);
#  int _getSearchAngleIndex(double angle);
#  double _getSearchAngle(int index);
#  void _populateSearchRays(int start, int end);
#  void _populateSearchAcross360(int first, int last);

#  void _makeFieldsUniform(size_t startIndex, size_t endIndex);

#};

#endif
