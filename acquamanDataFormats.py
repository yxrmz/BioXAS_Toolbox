# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 17:23:52 2020

@author: chernir
"""

format_BioXASMain = {'xaxis': 'AxisValues::',
                     'I0': 'I0Detector_darkCorrected',
                     'I1': 'I1Detector_darkCorrected',
                     'I2': 'I2Detector_darkCorrected',
                     'I3': 'PIPSDetector_darkCorrected',
                     'dwellTime': 'ScalerDwellTimeFeedback',
                     'mcaTotal': ('BioXASMainOutboardDetector',
                                  'BioXASMainInboardDetector'),
                     'mcaSingle': ('BioXASMainOutboardDetectorRawSpectrum{}',
                                   'BioXASMainInboardDetectorRawSpectrum{}')}

format_BioXASSide = {'xaxis': 'AxisValues::',
                     'I0': 'I0Detector_darkCorrected',
                     'I1': 'I1Detector_darkCorrected',
                     'I2': 'I2Detector_darkCorrected',
                     'I3': 'PIPSDetector_darkCorrected',
                     'dwellTime': 'ScalerDwellTimeFeedback',
                     'mcaTotal': ('Ge32Element',),
                     'mcaSingle': ('Ge32ElementRawSpectrum{}',)}

format_IDEAS = {'xaxis': 'AxisValues::',
                'I0': 'I_0',
                'I1': 'Sample',
                'I2': 'Reference',
                'dwellTime': 'dwellTime',
                'mcaTotal': ('KETEK',)}

format_VESPERS = {'xaxis': 'AxisValues::',
                  'I0': 'SplitIonChamber',
                  'I1': 'MiniIonChamber',
                  'I2': 'PostIonChamber',
                  'I3': 'PreKBIonChamber',
                  'dwellTime': 'MasterDwellTime',
                  'mcaTotal': ('FourElementXMapVortex',),
                  'mcaSingle': ('FourElementVotexRawSpectrum{}',)}


dataFormats = {'BioXAS-Main': format_BioXASMain,
               'BioXAS-Side': format_BioXASSide,
               'IDEAS': format_IDEAS,
               'VESPERS': format_VESPERS}