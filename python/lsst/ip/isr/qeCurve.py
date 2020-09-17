# This file is part of ip_isr.
#
# LSST Data Management System
# Copyright 2008-2017 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
"""
Something to do with transmission curves.
"""
from astropy.table import Table
from collections import defaultdict

import lsst.afw.image as afwImage
import lsst.afw.cameraGeom.utils as cgUtils
from lsst.geom import Point2I
from lsst.ip.isr import IsrCalib


__all__ = ["Curve"]


class QeCurve(IsrCalib):
    """Parameter set for quantum efficiency data.

    Parameters
    ----------
    curveType : `str`, optional
    level : `str`, optional
    listDict : `dict` [`str` `dict`]

    camera : `lsst.afw.cameraGeom.Camera`, optional
    detector : `lsst.afw.cameraGeom.Detector`, optional
        Detector object.  Passed to self.fromDetector() on init.
    log : `lsst.log.Log`, optional
        Logger to handle messages.
    kwargs : `dict`, optional
        Other keyword arguments to pass to the parent init.

    Raises
    ------
    RuntimeError :
        Raised if

    Notes
    -----
    The curve attributes stored are:

    curveType : `str`
        Type of curve stored.
    level : `str`
        Level at which the curve is defined.  Should be one of 'AMP',
        'DETECTOR', 'FP'.
    radiusScale : `float`
        Radius scale factor to apply.  Should be equal to the pixel size.
    units : `str`
        Units the wavelengths are tabulated in.
    curves : `dict` [`str`, `dict`]
        Dictionary of curve data, indexed by the name of the level
        segments.  Each contains the following keys:

        dimension : `int`
            Dimension of throughputs array.
        radii : `numpy.array`
            Array of radii samples.
        wavelengths : `numpy.array`
            Array of wavelength samples.
        throughputs : `numpy.array`
            Array of throughputs.
    """

    _OBSTYPE = 'qe_curve'
    _SCHEMA = 'Gen3 qe curves'
    _VERSION = 1.0

    def __init__(self, curveType=None, level=None, radiusScale=1.0, units='nm', **kwargs):
        self.curveType = curveType if curveType else 'QE'
        self.level = level if level else None
        self.radiusScale = radiusScale
        self.units = units
        self.curves = {}

        super().__init__(**kwargs)
        self.requiredAttributes.update(['curveType', 'level', 'radiusScale', 'units', 'curves'])

    def updateMetadata(self, setDate=False, **kwargs):
        """Update calibration metadata.

        This calls the base class's method after ensuring the required
        calibration keywords will be saved.

        Parameters
        ----------
        setDate : `bool`, optional
            Update the CALIBDATE fields in the metadata to the current
            time.  Defaults to False.
        kwargs :
            Other keyword parameters to set in the metadata.
        """
        kwargs['TYPE'] = self.curveType
        kwargs['LEVEL'] = self.level
        kwargs['RAD_SCL'] = self.radiusScale
        kwargs['UNITS'] = self.units

        super().updateMetadata(setDate=setDate, **kwargs)

    def addCurveData(key, dimension, radii=None, wavelengths=None, throughputs=None, length=(0, )):
        """Add curve data to the calibration dictionary.

        Parameters
        ----------
        key : `str`
            Index for the curve.
        dimension : `int`
            Dimension of the throughputs array.
        radii : `numpy.array`, optional
            Array of radii samples.
        wavelengths : `numpy.array`, optional
            Array of wavelength samples.
        throughputs : `numpy.array`, optional
            Array of throughputs.
        length : `tuple` [`int`], optional
            Array lengths to initialize to zero if no arrays are
            supplied.

        Raises
        ------
        RuntimeError :
            Raised if array lengths do not match.
        KeyError :
            Raised if the curve to create already exists.
        """
        if key in self.curves:
            raise KeyError(f"Key {key} already exists.")

        if dimension == 1:
            if radii:
                raise RuntimeError("Radii suppled for 1d curve.")

            if wavelengths and througputs:
                if len(wavelengths) != len(throughputs):
                    raise RuntimeError("Wavelength and throughput arrays differ in length.")
            elif len(length) == 1:
                radii = []
                wavelengths = np.zeros(shape=length)
                throughputs = np.zeros(shape=length)
            else:
                raise RuntimeError("No data supplied.")
        elif dimension == 2:
            if radii and wavelengths and throughputs:
                if througputs.shape != (len(wavelengths), len(radii)):
                    raise RuntimeError("Throughputs not equal to radii and wavelengths.")
            elif len(length) == 2:
                radii = np.zeros(length[1])
                wavelengths = np.zeros(length[0])
                throughputs = np.zeros(length)
            else:
                raise RuntimeError("No data supplied.")
        else:
            raise RuntimeError(f"Unsupported dimension {dimension}.")

        # Inputs should be validated now.
        self.curve[key] = {'dimension': dimension,
                           'radii' : radii,
                           'wavelengths' : wavelengths,
                           'throughputs' : througputs,
                           }

    @classmethod
    def fromDict(cls, dictionary):
        """Construct a calibration from a dictionary of properties.

        Parameters
        ----------
        dictionary : `dict`
            Dictionary of properties.

        Returns
        -------
        calib : `lsst.ip.isr.Curve`
            Constructed calibration.

        Raises
        ------
        RuntimeError :
            Raised if the supplied dictionary is for a different
            calibration.
        """
        calib = cls()

        metadata = dictionary['metadata']
        calib.setMetadata(metadata)

        calib.calibInfoFromDict(dictionary)

        calib.curveType = metadata.get('TYPE', 'QE')
        calib.level = metadata.get('LEVEL', metadata.get('MODE', 'DETECTOR'))
        calib.radiusScale = metadata.get('RAD_SCL', 1.0)
        calib.units = metadata.get('UNITS', 'nm')

        for key in dictionary['curves']:
            calib.curve[key] = dictionary['curves'][key]

        return calib

    def toDict(self):
        """Return a dictionary containing the calibration properties.

        Returns
        ----------
        outDict : `dict`
            Dictionary of properties.
        """
        self.updateMetadata()

        outDict = {}
        outDict['metadata'] = self.getMetadata()
        outDict['type'] = self.curveType
        outDict['level'] = self.level
        outDict['radiusScale'] = self.radiusScale
        outDict['units'] = self.units

        outDict.curve = {}
        for key in self.curves:
            outDict.curve[key] = self.curves[key]

        return outDict

    @classmethod
    def fromTable(cls, tableList):
        """Construct calibration from a list of tables.

        This method uses the `fromDict` method to create the
        calibration, after constructing an appropriate dictionary from
        the input tables.

        Parameters
        ----------
        tableList : `list` [`astropy.table.Table`]
            List of tables to use to construct the qe_curve
            calibration.

        Returns
        -------
        calib : `lsst.ip.isr.Curve`
            The calibration defined in the tables.
        """
        qeTable = tableList[0]
        metadata = qeTable.meta

        inDict = dict()
        inDict['metadata'] = metadata
        inDict['curves'] = defaultdict(dict)

        for row in qeTable:
            key = row['key']
            for column in ['dimension', 'radii', 'wavelengths' 'throughputs']:
                if column in row:
                    inDict['curves'][key][column] = row[column]

        return cls().fromDict(inDict)

    def toTable(self):
        """Construct a list of tables containing the information in this calibration.

        Thse list of tables should create an identical calibration
        after being passed to this class's fromTable method.

        Returns
        -------
        tableList : `list` [`astropy.table.Table`]
            List of tables containing the qe_curve calibration
            information.
        """
        self.updateMetadata()

        (keys, dimension, wavelengths, throughputs,
         unit, radii, radiusScale) = ([], [], [], [],
                                      [], [], [])
        keys = [key for key in self.curves]
        dimensions = [self.curve[key]['dimension'] for key in self.curves]
        radii = [self.curve[key]['radii'] for key in self.curves]
        wavelengths = [self.curve[key]['wavelengths'] for key in self.curves]
        throughputs = [self.curve[key]['throughputs'] for key in self.curves]

        catalog = Table([{'key': keys,
                          'radii': radii,
                          'dimension': dimension,
                          'wavelengths': wavelengths,
                          'throughputs': throughputs}])
        catalog.meta = self.getMetadata().toDict()

        return([catalog])

    # Science methods
    @staticmethod
    def standardCurveFromLists(wavelengths, throughputs, wavelengthUnit="nm", atMin=None, atMax=None):
        """Generate a curve from wavelength/throughput list.

        Parameters
        ----------
        wavelengths : `numpy.array`
            Wavelengths the curve is sampled at.
        throughputs : `numpy.array`
            Throughput at each wavelength sample.
        wavelengthUnit : `str`, optional
            Unit the wavelengths are tabulated in.  Defaults to 'nm'.
        atMin : `float`, optional
            Value to use at and below the minimum wavelength value.
        atMax : `float`, optional
            Value to use at and above the maximum wavelength value.

        Returns
        -------
        curve : `lsst.afw.image.TransmissionCurve`
            The transmission curve constructed.

        Raises
        ------
        RuntimeError :
            Raised if the wavelength and througput arrays have
            different lengths.
        """
        if len(wavelengths) != len(throughputs):
            raise RuntimeError("Wavelength and throughput lists are different length!")

        if wavelengthUnit == 'nm':
            wavelengthFactor = 10.0

        atMin = atMin if atMin else throughputs[0]
        atMax = atMax if atMax else throughputs[-1]
        return afwImage.TransmissionCurve.makeSpatiallyConstant(throughput=throughputs,
                                                                wavelengths=wavelengths * wavelengthFactor,,
                                                                throughputAtMin=atMin,
                                                                throughputAtMax=atMax)

    @staticmethod
    def radialCurveFromLists(radii, wavelengths, throughputs, radiusScale=1.0, wavelengthUnit="nm",
                             atMin=None, atMax=None):
        """Generate a curve from a list of wavelengths and radius points, and
        a grid of throughput values.

        Parameters
        ----------
        radii : `numpy.array`
            Radii the throughputs are sampled at.
        wavelengths : `numpy.array`
            Wavelengths the curve is sampled at.
        throughputs : `numpy.array`
            Throughput at each wavelength sample.
        wavelengthUnit : `str`, optional
            Unit the wavelengths are tabulated in.  Defaults to 'nm'.
        atMin : `float`, optional
            Value to use at and below the minimum wavelength value.
        atMax : `float`, optional
            Value to use at and above the maximum wavelength value.

        Returns
        -------
        curve : `lsst.afw.image.TransmissionCurve`
            The transmission curve constructed.

        Raises
        ------
        RuntimeError :
            Raised if the radii, wavelength, and througput arrays have
            different lengths.
        """
        if len(radii) != throughputs.shape[0]:
            raise RuntimeError(f"Radii dimension mismatch: {len(radii)} != {throughputs.shape[0]}")
        if len(wavelengths) != throughputs.shape[1]:
            raise RuntimeError(f"Wavelength dimension mismatch: {len(wavelengths)} != {throughputs.shape}")

        if wavelengthUnit == 'nm':
            wavelengthFactor = 10.0

        atMin = atMin if atMin else throughputs[0]
        atMax = atMax if atMax else throughputs[-1]
        return afwImage.TransmissionCurve.makeRadial(throughput=throughputs,
                                                     wavelengths=wavelengths * wavelengthFactor,
                                                     radii=radii / radiusScale,
                                                     throughputAtMin=0.0,
                                                     throughputAtMax=0.0)

    def evaluate(self, key):
        """Generate the afwImage.TransmissionCurve from the calibration
        dictionary.

        Parameters
        ----------
        key : `str`
            Name of the curve to construct.

        Returns
        -------
        curve : `lsst.afw.image.TransmissionCurve`
            The constructed curve.

        Raises
        ------
        KeyError
            Raised if the requested key does not exist.
        """
        if key not in self.curves:
            raise KeyError(f"No entry {key} found in {self}")

        curveData = self.curves[key]
        if curveData['dimension'] == 1:
            return self.standardCurveFromLists(curveData['wavelengths'], curveData['throughputs'],
                                               wavelengthUnit=self.units,
                                               atMin=curveData['atMin'],
                                               atMax=curveData['atMax'])
        elif curveData['dimension'] == 2:
            return self.radialCurveFromLists(curveData['radii'],
                                             curveData['wavelengths'], curveData['throughputs'],
                                             wavelengthUnit=self.units,
                                             atMin=curveData['atMin'],
                                             atMax=curveData['atMax'])
        else:
            raise RuntimeError(f"Unsupported curve dimension {curveData['dimension']}")
