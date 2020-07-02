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
Arbitrary curves with interpolation.
"""
import astropy.units as u
import numpy as np
import itertools

from astropy.table import Table
from collections import defaultdict
from scipy.interpolate import interp1d

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
    _SCHEMA = 'Simple Curve'
    _VERSION = 1.0

    def __init__(self, mode="UNKNOWN", **kwargs):
        self.type = 'QE'
        self.mode = mode

        super().__init__(**kwargs)

        self.requiredAttributes.update(['mode', 'type',
                                        'wavelength', 'efficiency'])

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
        kwargs['MODE'] = self.mode or "UNKNOWN"
        kwargs['TYPE'] = self.type

        super().updateMetadata(setDate=setDate, **kwargs)

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

        if calib._OBSTYPE != dictionary['metadata']['OBSTYPE']:
            raise RuntimeError(f"Incorrect curve supplied.  Expected {calib._OBSTYPE}, "
                               f"found {dictionary['metadata']['OBSTYPE']}")

        calib.setMetadata(dictionary['metadata'])
        calib._detectorName = dictionary['metadata']['DETECTOR']
        calib._detectorSerial = dictionary['metadata'].get('DETECTOR_SERIAL', None)
        calib.mode = dictionary['metadata']['MODE']

        calib.data = dictionary['data']

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
        inDict['mode'] = metadata['MODE']
        inDict['type'] = metadata['TYPE']
        inDict['data'] = defaultdict(dict)

        wavelength = qeTable['wavelength']
        efficiency = qeTable['efficiency']

        if inDict['mode'] == 'AMP':
            ampNameList = qeTable['amp_name']
            ampNames = set(ampNameList)
            for ampName in ampNames:
                idx = np.where(ampNameList == ampName)[0]
                inDict['data'][ampName]['wavelength'] = wavelength[idx]
                inDict['data'][ampName]['efficiency'] = efficiency[idx]

        elif inDict['mode'] in ('DETECTOR', 'IMAGE'):
            inDict['data'][inDict['mode']]['wavelength'] = wavelength
            inDict['data'][inDict['mode']]['efficiency'] = efficiency

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
        if self.mode == 'AMP':
            amp_names = [[ampName] * len(self.data[ampName]) for ampName in self.data]
            wavelength = list(itertools.chain.from_iterable([self.data[ampName]['wavelength']
                                                             for ampName in self.data]))
            efficiency = list(itertools.chain.from_iterable([self.data[ampName]['efficiency']
                                                             for ampName in self.data]))
            catalog = Table([{'amp_names': amp_names, 'wavelength': wavelength,
                              'efficiency': efficiency}])
        elif self.mode in ('DETECTOR', 'IMAGE'):
            catalog = Table([{'wavelength': self.data[self.mode]['wavelength'],
                              'efficiency': self.data[self.mode]['efficiency']}])

        catalog.meta = self.getMetadata().toDict()
        return([catalog])

    # Actual science methods
    def evaluate(self, wavelength, detector=None, position=None,
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
            wavelengths = self.data[targetAmp.getName()]['wavelength']
            efficiencies = self.data[targetAmp.getName()]['efficiency']
        elif self.mode in ('DETECTOR', 'IMAGE'):
            # Should IMAGE have subregions defined?
            wavelengths = self.data[self.mode]['wavelength']
            efficiencies = self.data[self.mode]['efficiency']

        return self.interpolate(wavelengths, efficiencies, wavelength,
                                kind=kind, bounds_error=bounds_error, fill_value=fill_value)

    def interpolate(self, wavelengths, values, wavelength, kind, bounds_error, fill_value):
        """Interplate the curve at the specified wavelength(s).
        Parameters
        ----------
        wavelengths : `astropy.units.Quantity`
            The wavelength values for the curve.
        values : `astropy.units.Quantity`
            The y-values for the curve.
        wavelength : `astropy.units.Quantity`
            The wavelength(s) at which to make the interpolation.
        kind : `str`
            The type of interpolation to do.  See documentation for
            `scipy.interpolate.interp1d` for accepted values.
        Returns
        -------
        value : `astropy.units.Quantity`
            Interpolated value(s)
        """
        if not isinstance(wavelength, u.Quantity):
            raise ValueError("Wavelengths at which to interpolate must be astropy quantities")
        if not (isinstance(wavelengths, u.Quantity) and isinstance(values, u.Quantity)):
            raise ValueError("Model to be interpreted must be astropy quantities")
        interp_wavelength = wavelength.to(wavelengths.unit)
        f = interp1d(wavelengths, values,
                     kind=kind, bounds_error=bounds_error, fill_value=fill_value)
        return f(interp_wavelength.value)*values.unit
