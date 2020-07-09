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


class Curve(IsrCalib):
    """An arbitrary curve with interpolation.

    Parameters
    ----------
    log : `lsst.log.Log`, optional
        Log to write messages to.
    **kwargs :
        Parameters to pass to parent constructor.
    """

    _OBSTYPE = 'qe_curve'
    _SCHEMA = 'Gen3 qe curves'
    _VERSION = 1.0

    def __init__(self, curveType=None, level=None, filterName=None, listDict=None,
                 **kwargs):
        self.curveType = curveType if curveType else 'QE'
        self.level = level if level else None
        self.filter = filterName if filterName else None
        self.curves = {}

        super().__init__(**kwargs)
        self.requiredAttributes.update(['curveType', 'level', 'filter'
                                        'curves'])

        if listDict:
            for key in listDict:
                wavelengthUnit = listDict[key].get('unit', "nm")
                wavelengths = listDict[key].get('wavelengths')
                throughputs = listDict[key].get('throughputs')
                self.curve[key] = self.fromLists(wavelengths, throughputs,
                                                 wavelengthUnit=wavelengthUnit)

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
        kwargs['FILTER'] = self.filter

        super().updateMetadata(setDate=setDate, **kwargs)

    @staticmethod
    def standardCurveFromLists(wavelengths, throughputs,
                               wavelengthUnit="nm"):
        """Method to generate a curve from wavelength/throughput list
        """
        if len(wavelengths) != len(throughputs):
            raise RuntimeError("Wavelength and throughput lists are different length!")

        if wavelengthUnit == 'nm':
            wavelengths *= 10.0
        atMin = throughputs[0]
        atMax = throughputs[-1]
        return afwImage.TransmissionCurve.makeSpatiallyConstant(throughput=throughputs,
                                                                wavelengths=wavelengths,
                                                                throughputAtMin=atMin,
                                                                throughputAtMax=atMax)

    @staticmethod
    def radialCurveFromLists(radii, wavelengths, throughputs,
                             radiusScale=1.0, wavelengthUnit="nm"):
        """Method to generate a curve from a set of throughputs defined on a
        grid of wavelength/radius.

        """
        if len(radii) != throughputs.shape[0]:
            raise RuntimeError(f"Radii dimension mismatch: {len(radii)} != {throughputs.shape[0]}")
        if len(wavelengths) != throughputs.shape[1]:
            raise RuntimeError(f"Wavelength dimension mismatch: {len(wavelengths)} != {throughputs.shape}")

        radii *= radiusScale
        if wavelengthUnit == 'nm':
            wavelengths *= 10.0

        return afwImage.TransmissionCurve.makeRadial(throughput=throughputs,
                                                     wavelengths=wavelengths,
                                                     radii=radii,
                                                     throughputAtMin=0.0,
                                                     throughputAtMax=0.0)

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
        if calib._OBSTYPE != metadata['OBSTYPE']:
            raise RuntimeError(f"Incorrect curve supplied.  Expected {calib._OBSTYPE}, "
                               f"found {metadata['OBSTYPE']}")

        calib.setMetadata(metadata)
        calib.curveType = metadata.get('TYPE', 'QE')
        calib.level = metadata.get('LEVEL', metadata.get('MODE', 'DETECTOR'))
        calib.filter = metadata.get('FILTER', None)

        calib._detectorName = metadata.get('DETECTOR', None)
        calib._detectorSerial = metadata.get('DETECTOR_SERIAL', None)

        for key in dictionary['curves']:
            curveDict = dictionary['curves']['key']
            if curveDict['dimension'] == 1:
                calib.curve[key] = cls.standardCurvefromLists(curveDict['wavelengths'],
                                                              curveDict['throughputs'],
                                                              wavelengthUnit=curveDict.get('unit', 'nm'))
            elif curveDict['dimension'] == 2:
                calib.curve[key] = cls.radialCurveFromLists(curveDict['radii'],
                                                            curveDict['wavelengths'],
                                                            curveDict['throughputs'],
                                                            curveDict.get('radiusScale', 1.0),
                                                            curveDict.get('unit', 'nm'))
            else:
                raise RuntimeError(f"Unknown  dimension: {curveDict['dimension']}")

        return calib

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
            for column in ['dimension', 'wavelengths' 'throughputs',
                           'unit', 'radii', 'radiusScale']:
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
        for key in self.curves:
            keys.append(key)
            dimension.append(self.curves[key]['dimension'])
            wavelengths.append(self.curves[key]['wavelengths'])
            throughputs.append(self.curves[key]['throughputs'])
            unit.append(self.curves[key].get('unit', 'nm'))
            radii.append(self.curves[key].get('radii', []))
            radiusScale.append(self.curves[key].get('radiusScale', 0.0))

        catalog = Table([{'key': keys,
                          'dimension': dimension,
                          'wavelengths': wavelengths,
                          'throughputs': throughputs,
                          'unit': unit,
                          'radii': radii,
                          'radiusScale': radiusScale}])
        catalog.meta = self.getMetadata().toDict()

        return([catalog])

    # Actual science methods
    def evaluate(self, wavelengths, detector=None, position=None,
                 kind='linear', bounds_error=False, fill_value=0):
        """Interpolate the curve at the specified position and wavelength.

        Parameters
        ----------
        wavelength : `astropy.units.Quantity`
            The wavelength(s) at which to make the interpolation.
        position : `lsst.geom.Point2D`, optional
            The position on the detector at which to evaluate the curve.
        detector : `lsst.afw.cameraGeom.Detector`, optional
            Is used to find the appropriate curve given the position for
            curves that vary over the detector.  Ignored in the case where
            there is only a single curve per detector.
        kind : `str`, optional
            The type of interpolation to do (default is 'linear').
            See documentation for `scipy.interpolate.interp1d` for
            accepted values.
        bounds_error : `bool`, optional
            Raise error if interpolating outside the range of x?
            (default is False)
        fill_value : `float`, optional
            Fill values outside the range of x with this value
            (default is 0).
        Returns
        -------
        value : `astropy.units.Quantity`
            Interpolated value(s).  Number of values returned will match the
            length of `wavelength`.

        Raises
        ------
        ValueError
            If the ``bounds_error`` is changed from the default, it will raise
            a `ValueError` if evaluating outside the bounds of the curve.
        """
        if self.mode == 'AMP':
            targetAmp = cgUtils.findAmp(detector, Point2I(position))
            return self.curves[targetAmp].sampleAt(position, wavelengths)
        elif self.mode in ('DETECTOR', 'IMAGE'):
            # Should IMAGE have subregions defined?
            keys = list(self.curves.keys())
            if len(keys) == 1:
                return self.curves[keys[0]].sampleAt(position, wavelengths)
            else:
                raise RuntimeError("unexpected multi-key value!")
