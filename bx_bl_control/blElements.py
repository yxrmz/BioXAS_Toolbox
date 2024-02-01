# -*- coding: utf-8 -*-
"""
Created on Sun Jul 23 13:20:02 2017

@author: bioxasmain
"""


import epics
import time
import numpy as np
import inspect
from collections import OrderedDict, namedtuple
import sys
from datetime import datetime
import subprocess


myInput = raw_input if sys.version[0] == '2' else input

hc_ = 12398.42
crystal2D_ = 3.8403117

pCenter = namedtuple('pCenter', ['x', 'y', 'z'])
pOffset = namedtuple('pOffset', ['x', 'y', 'z'])
AngleOffset = namedtuple('AngleOffset', ['pitch', 'roll', 'yaw'])

configurationMain = {
    'pseudoM1pitch': '',
    'pseudoM1roll': '',
    'pseudoM1yaw': '',
    'pseudoM1height': '',
    'pseudoM2pitch': '',
    'pseudoM2roll': '',
    'pseudoM2yaw': '',
    'pseudoM2height': '',
    'monoBraggPosition': '',
    'monoHeight': '',
    'BeWindowHeight': '',
    'pseudoTableHeight': ''}

configurationSide = {
    'pseudoM1pitch': '',
    'pseudoM1roll': '',
    'pseudoM1yaw': '',
    'pseudoM1height': '',
    'pseudoM2pitch': '',
    'pseudoM2roll': '',
    'pseudoM2yaw': '',
    'pseudoM2height': '',
    'monoBraggPosition': '',
    'monoHeight': '',
    'BeWindowHeight': '',
    'pseudoTableHeight': ''}
# regionOffset_ = 180
#
# goniometerMotor_ = 0
# region_ = 0
# m1MirrorPitch_ = 0


class Motor(object):
    def __init__(self, pvStr=None, fieldStr=None, fbk=None, description=None):
        self.status = 0
        self.position = 0
        self.__pvReader = epics.PV(
            '{0}:{1}{2}'.format(pvStr, fieldStr,
                                '' if fbk is None else ':' + fbk),
            callback=self.monitor_position, auto_monitor=True)
        self.__pvWriter = epics.PV(pvStr+':'+fieldStr)
        self.__pvStatus = epics.PV(pvStr+':status',
                                   callback=self.monitor_status,
                                   auto_monitor=True)
        self.__pvStop = epics.PV(pvStr+':stop')

    def monitor_status(self, pvname=None, value=None, char_value=None, **kw):
        self.status = value

    def monitor_position(self, pvname=None, value=None, char_value=None, **kw):
        self.position = value

    def get_position(self):
        self.position = self.__pvReader.get()
        return self.position

    def get_status(self):
        self.status = self.__pvStatus.get()
        return self.status

    def move(self, position):
        self.status = 1
        self.__pvWriter.put(position)

    def stop(self):
        self.__pvStop.put(1)


class MotorK(object):
    def __init__(self, pvStr=None, fieldStr=None, fbk=None, description=None,
                 kStr=None, kOffset=0):
        self.status = 0
        self.position = 0
        self.positionK = 0
        self.kOffset = kOffset
        self.__pvReader = epics.PV(
            '{0}:{1}{2}'.format(pvStr, fieldStr,
                                '' if fbk is None else ':' + fbk),
            callback=self.monitor_position, auto_monitor=True)
        self.__pvReaderK = epics.PV(kStr,
            callback=self.monitor_positionK, auto_monitor=True)
        self.__pvWriter = epics.PV(pvStr+':'+fieldStr)
        self.__pvStatus = epics.PV(pvStr+':status',
                                   callback=self.monitor_status,
                                   auto_monitor=True)
        self.__pvStop = epics.PV(pvStr+':stop')

    def monitor_status(self, pvname=None, value=None, char_value=None, **kw):
        self.status = value

    def monitor_position(self, pvname=None, value=None, char_value=None, **kw):
        self.position = value

    def monitor_positionK(self, pvname=None, value=None, char_value=None, **kw):
        self.positionK = -value + self.kOffset

    def get_position(self):
        self.position = self.__pvReader.get()
        return self.position

    def get_positionK(self):
        self.positionK = self.__pvReaderK.get()
        return self.positionK

    def get_status(self):
        self.status = self.__pvStatus.get()
        return self.status

    def move(self, position):
        self.status = 1
        self.__pvWriter.put(position)

    def stop(self):
        self.__pvStop.put(1)


class DCMenergyMotor(object):
    def __init__(self, pvStr=None, fieldStr=None, fbk=None, description=None,
                 mirror1pitchMotor=None):
        self.status = 0
        self.position = 0
#        print pvStr, fieldStr, fbk
        self.__pvReader = epics.PV(
            '{0}:{1}{2}'.format(pvStr, fieldStr,
                                '' if fbk is None else ':' + fbk),
            callback=self.monitor_position, auto_monitor=True)
        self.__pvWriter = epics.PV(pvStr+':'+fieldStr)
        self.__pvStatus = epics.PV(pvStr+':status',
                                   callback=self.monitor_status,
                                   auto_monitor=True)
        self.__pvStop = epics.PV(pvStr+':stop')
        self.region = 'B'
        self.regionOffset = 0
        self.thetaBraggOffset = 180.0
        self.M1PV = mirror1pitchMotor
        self.positionMax = 33000
        self.positionMin = 4000
        self.positionGonio = 0
        self.positionBragg = 0

    def monitor_position(self, pvname=None, value=None, char_value=None, **kw):
        self.positionGonio = value
        self.positionBragg = self.bragg_from_position(self.positionGonio)
        self.position = self.energy_from_bragg(self.positionBragg)

    def monitor_status(self, pvname=None, value=None, char_value=None, **kw):
        self.status = value

    def get_position(self):
        self.positionGonio = self.__pvReader.get()
#        print self.positionGonio
        self.positionBragg = self.bragg_from_position(self.positionGonio)
        self.position = self.energy_from_bragg(self.positionBragg)
        return self.position

    def get_status(self):
        self.status = self.__pvStatus.get()
        return self.status

    def move(self, energy):
        if energy > self.positionMin and energy < self.positionMax:
            self.status = 1
    #        print self.position_from_bragg(self.bragg_from_energy(energy))
            self.__pvWriter.put(self.position_from_bragg(
                self.bragg_from_energy(energy)))

    def moveBragg(self, braggAngle):
        self.status = 1
        self.__pvWriter.put(self.position_from_bragg(braggAngle))

    def bragg_from_position(self, goniometerPosition):
        braggAngle = self.thetaBraggOffset - goniometerPosition -\
            (2*self.M1PV.position)
        if self.region == 'A':
            braggAngle -= self.regionOffset
        return braggAngle

    def position_from_bragg(self, braggAngle):
        braggPosition = self.thetaBraggOffset - (2*self.M1PV.position) -\
            braggAngle
        if self.region == 'A':
            braggPosition -= self.regionOffset
        return braggPosition

    def bragg_from_energy(self, energy):
        return np.degrees(np.arcsin(hc_ / (energy * crystal2D_)))

    def energy_from_bragg(self, braggAngle):
        return hc_ / (crystal2D_ * np.sin(np.radians(braggAngle)))

    def stop(self):
        self.__pvStop.put(1)


class DummyMotor(object):
    def __init__(self, pvStr=None, fieldStr=None, fbk=None, description=None):
        self.status = 0
        self.position = 0

    def monitor_status(self):
        self.status = 0

    def check_status(self):
        self.status = 0

    def monitor_position(self):
        self.status = 0

    def move(self, position):
        self.status = 0
        self.position = position
        self.status = 0

    def stop(self):
        self.status = 0


class Counter(object):
    def __init__(self, pvStr=None, channels=None):
        self.status = 0
        self._pvStr = pvStr
        self.__pvStarter = epics.PV(pvStr+':startScan')
        self.__pvTimes = epics.PV(pvStr+':delay')
        self.__pvMode = epics.PV(pvStr+':inputMode')
        self.__pvTrigger = epics.PV(pvStr+':triggerSource')
        self.chIndex = {}
        for iChannel, channel in enumerate(channels):
            self.__pvChannel = epics.PV(
                '{0}{1}:fbk'.format(pvStr, channel),
                callback=self.monitor_channel,
                auto_monitor=True)
            self.chIndex[str(channel)] = iChannel

        self.data = [0]*len(channels)
        self.updateMonitor = [1]*len(channels)
        self.__pvScanCount = epics.PV(pvStr+':scanCount')

    def monitor_channel(self, pvname=None, value=None, char_value=None, **kw):
        channelNum = str(pvname).split(':')[1][3:]
        self.data[self.chIndex[channelNum]] = value
        self.updateMonitor[self.chIndex[channelNum]] = 0
#        print pvname, value, time.ctime()

    def monitor_status(self, pvname=None, value=None, char_value=None, **kw):
        self.status = value

    def set_scan_count(self, value):
        self.__pvScanCount.put(value)

    def start(self):
        self.__pvStarter.put(1)

    def stop(self):
        self.__pvStarter.put(0)

    def clear(self):
        self.data = [0]*len(self.data)
        self.updateMonitor = [1]*len(self.data)

    def read(self):
        return self.data

    def set_time(self, timeMS):
        self.__pvTimes.put(timeMS)

    def set_mode(self, mode):
        self.__pvMode.put(mode)

    def set_trigger(self, trigger):
        self.__pvTrigger.put(trigger)


class GateGenerator(object):
    def __init__(self, pvStr):
        self.__pvGate = epics.PV(pvStr+':SOFT_IN:B0')
        self.__pvPulseWidth = [epics.PV(
            pvStr+':PULSE{}_WID'.format(ch+1)) for ch in range(4)]

    def set_pulse_width(self, channel, val):
        self.__pvPulseWidth[channel-1].put(val)

    def sendGate(self):
        self.__pvGate.put(1)


class ISEGVHQ(object):
    def __init__(self, pvStr=None, fieldStr=None, fbk=None, description=None):
        self.status = 0
        self.position = 0
        self.__pvReader = epics.PV(
            '{0}:{1}{2}'.format(pvStr, fieldStr,
                                '' if fbk is None else ':' + fbk),
            callback=self.monitor_position, auto_monitor=True)
#        self.__pvWriter = epics.PV(pvStr+':'+fieldStr)
        self.__pvStatus = epics.PV(pvStr+':VoutStatus',
                                   callback=self.monitor_status,
                                   auto_monitor=True)
        self.__pvStart = epics.PV(pvStr+':set_start')
        self.__pvStop = epics.PV(pvStr+':stop')

    def monitor_status(self, pvname=None, value=None, char_value=None, **kw):
        self.status = value

    def check_status(self):
        self.status = self.__pvStatus.get()

    def monitor_position(self, pvname=None, value=None, char_value=None, **kw):
        self.position = value

    def check_position(self):
        self.position = self.__pvReader.get()

    def move(self, position):
        self.status = 1
#        self.__pvWriter.put(position)
        self.__pvStart.put(position)

    def stop(self):
        self.__pvStop.put(1)


class DCM(object):
    def __init__(self, bl=None, braggPVStr=None, heightPVStr=None):
                # , yawPVStr=None,
        #             heightPVStr=None, lateralPVStr=None, bendPVStr=None):
        self.beamLine = bl
        paramList = inspect.getargspec(self.__init__)[0][2:]
        localVars = locals()
        initValuesList = [localVars[iP] for iP in paramList]
        paramDict = dict(zip(paramList, initValuesList))
        for pKey, pVal in paramDict.items():
            tmpParaName = '_{}'.format(pKey[:-3])
#            print tmpParaName[1:-2]
            if pVal is not None:
                pvStrList = pVal.split(':')
                if tmpParaName == '_braggPV':
                    tMotor = DCMenergyMotor(
                        pvStr=':'.join(pvStrList[:-2]),
                        fieldStr=pvStrList[-2],
                        fbk=pvStrList[-1],
                        mirror1pitchMotor=self.beamLine.M1._pitchPV)
                    tMotor.get_position()
                    setattr(self, tmpParaName, tMotor)
                else:
                    tMotor = Motor(pvStr=':'.join(pvStrList[:-2]),
                                   fieldStr=pvStrList[-2],
                                   fbk=pvStrList[-1])
                    tMotor.get_position()
                    setattr(self, tmpParaName, tMotor)
            else:
                setattr(self, tmpParaName, None)

    @property
    def bragg(self):
        return None if self._braggPV is None else self._braggPV.positionBragg

    @bragg.setter
    def bragg(self, pos):
        return None if self._braggPV is None else self._braggPV.moveBragg(pos)

    @property
    def height(self):
        return None if self._heightPV is None else self._heightPV.position

    @height.setter
    def height(self, pos):
        return None if self._heightPV is None else self._heightPV.move(pos)

    @property
    def energy(self):
        return None if self._braggPV is None else self._braggPV.position

    @energy.setter
    def energy(self, pos):
        return None if self._braggPV is None else self._braggPV.move(pos)


class Mirror(object):
    def __init__(self, pitchPVStr=None, rollPVStr=None, yawPVStr=None,
                 heightPVStr=None, lateralPVStr=None, bendPVStr=None,
                 maskPVStr=None):
        self.listOfPVs = []
        self.center = pCenter(0, 0, 0)
        self.positionOffset = pOffset(0, 0, 0)
        self.angularOffset = AngleOffset(0, 0, 0)
        paramList = inspect.getargspec(self.__init__)[0][1:]
        localVars = locals()
        initValuesList = [localVars[iP] for iP in paramList]
        paramDict = dict(zip(paramList, initValuesList))
        for pKey, pVal in paramDict.items():
            tmpParaName = '_{}'.format(pKey[:-3])
#            print tmpParaName[1:-2]
            if pVal is not None:
                pvStrList = pVal.split(':')
                tMotor = Motor(pvStr=':'.join(pvStrList[:-2]),
                               fieldStr=pvStrList[-2],
                               fbk=pvStrList[-1])
                tMotor.get_position()
                setattr(self, tmpParaName, tMotor)
                self.listOfPVs.append(tmpParaName)
            else:
                setattr(self, tmpParaName, None)

    @property
    def pitch(self):
        return None if self._pitchPV is None else\
            self._pitchPV.position + self.angularOffset.pitch

    @pitch.setter
    def pitch(self, pos):
        print("Actual position will be set to", pos - self.angularOffset.pitch)
        return None if self._pitchPV is None else\
            self._pitchPV.move(pos - self.angularOffset.pitch)

    @property
    def roll(self):
        return None if self._rollPV is None else self._rollPV.position

    @roll.setter
    def roll(self, pos):
        return None if self._rollPV is None else self._rollPV.move(pos)

    @property
    def yaw(self):
        return None if self._yawPV is None else self._yawPV.position

    @yaw.setter
    def yaw(self, pos):
        return None if self._yawPV is None else self._yawPV.move(pos)

    @property
    def height(self):
        return None if self._heightPV is None else\
            self._heightPV.position  # - self.positionOffset.z

    @height.setter
    def height(self, pos):
        return None if self._heightPV is None else\
            self._heightPV.move(pos)

    @property
    def lateral(self):
        return None if self._lateralPV is None else self._lateralPV.position

    @lateral.setter
    def lateral(self, pos):
        return None if self._lateralPV is None else self._lateralPV.move(pos)

    @property
    def bend(self):
        return None if self._bendPV is None else self._bendPV.position

    @bend.setter
    def bend(self, pos):
        return None if self._bendPV is None else self._bendPV.move(pos)

    @property
    def mask(self):
        return None if self._maskPV is None else self._maskPV.position

    @mask.setter
    def mask(self, pos):
        return None if self._maskPV is None else self._maskPV.move(pos)
#    @property
#    def status(self):
#        statusVar = 0
#        for motor in self.listOfPVs:
#            motorPV = getattr(self, motor)
#            statusVar += int(motorPV.status)
#        return statusVar


class Table(object):
    def __init__(self, heightPVStr=None, pitchPVStr=None, rollPVStr=None):
        paramList = inspect.getargspec(self.__init__)[0][1:]
        localVars = locals()
        initValuesList = [localVars[iP] for iP in paramList]
        paramDict = dict(zip(paramList, initValuesList))
        for pKey, pVal in paramDict.items():
            tmpParaName = '_{}'.format(pKey[:-3])
#            print tmpParaName[1:-2]
            if pVal is not None:
                pvStrList = pVal.split(':')
                tMotor = Motor(pvStr=':'.join(pvStrList[:-2]),
                               fieldStr=pvStrList[-2],
                               fbk=pvStrList[-1])
                tMotor.get_position()
                setattr(self, tmpParaName, tMotor)
            else:
                setattr(self, tmpParaName, None)

    @property
    def pitch(self):
        return None if self._pitchPV is None else self._pitchPV.position

    @pitch.setter
    def pitch(self, pos):
        return None if self._pitchPV is None else self._pitchPV.move(pos)

    @property
    def roll(self):
        return None if self._rollPV is None else self._rollPV.position

    @roll.setter
    def roll(self, pos):
        return None if self._rollPV is None else self._rollPV.move(pos)

    @property
    def height(self):
        return None if self._heightPV is None else self._heightPV.position

    @height.setter
    def height(self, pos):
        return None if self._heightPV is None else self._heightPV.move(pos)


class BeWindow(object):
    def __init__(self, heightPVStr=None):
        paramList = inspect.getargspec(self.__init__)[0][1:]
        localVars = locals()
        initValuesList = [localVars[iP] for iP in paramList]
        paramDict = dict(zip(paramList, initValuesList))
        for pKey, pVal in paramDict.items():
            tmpParaName = '_{}'.format(pKey[:-3])
#            print tmpParaName[1:-2]
            if pVal is not None:
                pvStrList = pVal.split(':')
                tMotor = Motor(pvStr=':'.join(pvStrList[:-2]),
                               fieldStr=pvStrList[-2],
                               fbk=pvStrList[-1])
                tMotor.get_position()
                setattr(self, tmpParaName, tMotor)
            else:
                setattr(self, tmpParaName, None)

    @property
    def height(self):
        return None if self._heightPV is None else self._heightPV.position

    @height.setter
    def height(self, pos):
        return None if self._heightPV is None else self._heightPV.move(pos)


class BeamLine(object):
    def __init__(self, configuration='Main'):
        self.config = configurationMain if configuration == 'Main' else\
            configurationSide
        self.M1 = Mirror(pitchPVStr=self.config['pseudoM1pitch'],
                         rollPVStr=self.config['pseudoM1roll'],
                         yawPVStr=self.config['pseudoM1yaw'],
                         heightPVStr=self.config['pseudoM1height'],
                         bendPVStr=self.config['pseudoM1bend'],
                         maskPVStr=self.config['M1Mask'])
        self.M1.center = self.config['M1Center']
        self.M1.offset = self.config['M1Offset']
        self.M1.angularOffset = self.config['M1AngularOffset']
        self.M1.bending_curve = self.config['M1BendingCurve']
        self.M1.pitch_curve = self.config['M1PitchCurve']
        self.M2 = Mirror(pitchPVStr=self.config['pseudoM2pitch'],
                         rollPVStr=self.config['pseudoM2roll'],
                         yawPVStr=self.config['pseudoM2yaw'],
                         heightPVStr=self.config['pseudoM2height'],
                         bendPVStr=self.config['pseudoM2bend'])
        self.M2.center = self.config['M2Center']
        self.M2.offset = self.config['M2Offset']
        self.M2.angularOffset = self.config['M2AngularOffset']
        self.M2.bending_curve = self.config['M2BendingCurve']
        self.M2.pitch_curve = self.config['M1PitchCurve']
        self.DCM = DCM(bl=self,
                       braggPVStr=self.config['monoBraggPosition'],
                       heightPVStr=self.config['monoHeight'])
        self.DCM.center = self.config['monoCenter']
        self.DCM.offset = self.config['monoOffset']
        self.BeWindow = BeWindow(heightPVStr=self.config['BeWindowHeight'])
        self.BeWindow.center = self.config['BeWindowCenter']
        self.BeWindow.offset = self.config['BeWindowOffset']
        self.Table = Table(heightPVStr=self.config['pseudoTableHeight'])
        self.Table.center = self.config['tableCenter']
        self.Table.offset = self.config['tableOffset']
        self.Table.angularOffset = self.config['tableAngularOffset']

#        try:
#            self.__pvMon = epics.PV(
#                'BIOXAS-SP-BL-DEV:AO',
#                callback=self.monitor_in,
#                auto_monitor=True)
#        except:
#            pass

    @property
    def energy(self):
        return self.DCM.energy

    @energy.setter
    def energy(self, value):
        self.DCM.energy = value

    def monitor_in(self, pvname=None, value=None, char_value=None, **kw):
        print(value)
        print("Are you sure you want to reposition beamline?")

    def calculate_mirror_bendmm(self, theta, mirror):
        """
        theta: desired pitch angle in degrees
        mirror: pointer to the instanse of the Mirror class
        focus: absolute longitudinal coordinate of the focal spot (usually the
        JJ slit). If None, considered to be the source.
        """
        def get_Rmer_from_Coddington(p, q, theta=None):  # from xrt.oes_base
            return 2 * p * q / (p+q) / np.sin(np.radians(theta))

        focus = 0 if mirror == self.M1 else 30350  # self.JJslit.center.y
        Rmer = get_Rmer_from_Coddington(np.abs(mirror.center.y - focus),
                                        1e12, theta) * 1e-3
        print("bending radius =", Rmer)
        return mirror.bending_curve(Rmer) if hasattr(mirror, 'bending_curve')\
            else None

    def calculatePositions(self, newEnergy):
        returnDict = OrderedDict()
        eStart = newEnergy - 200
        kyaw = 0.5882352941176472

        if eStart < 10500:
            m1Pitch = 0.16
            m1Bend = 75
            # DBHR is required
        elif eStart < 14000:
            m1Pitch = 0.205
            m1Bend = 50
#        elif eStart < 13200:
#            m1Pitch = 0.20
#        elif eStart < 14000:
#            m1Pitch = 0.19
#        elif eStart < 14400:
#            m1Pitch = 0.18
#        elif eStart < 15200:
#            m1Pitch = 0.18
#        elif eStart < 16300:
#            m1Pitch = 0.17
#        elif eStart < 17400:
        elif eStart < 22400:
            m1Pitch = 0.16
            m1Bend = 75
#        elif eStart < 18600:
#            m1Pitch = 0.14
#        elif eStart < 20400:
#            m1Pitch = 0.13
#        elif eStart < 22500:
#            m1Pitch = 0.12
        else:
            m1Pitch = 0.12  # M1 Pitch limit 0.1
            m1Bend = 77

        returnDict['mirrorPitch'] = m1Pitch
        returnDict['M1Bend'] = m1Bend
        returnDict['monoHeight'] = (self.DCM.center.y - self.M1.center.y) *\
            np.tan(np.radians(2 * m1Pitch)) - self.DCM.offset.z
#        returnDict['M1Bend'] = self.calculate_mirror_bendmm(m1Pitch, self.M1)
        monoVShift = 2 * 6.5 * np.cos(
            np.arcsin(hc_ / (newEnergy * crystal2D_)))

        M2AbsHeight = (self.M2.center.y - self.M1.center.y) *\
            np.tan(np.radians(2 * m1Pitch)) + monoVShift
        returnDict['M2Height'] = M2AbsHeight - self.M2.offset.z
#        returnDict['M2Bend'] = self.calculate_mirror_bendmm(m1Pitch, self.M2) - 0.1
        returnDict['BeWindowHeight'] = M2AbsHeight - self.BeWindow.offset.z
        if returnDict['BeWindowHeight'] < -21:
            returnDict['BeWindowHeight'] = -21
        elif returnDict['BeWindowHeight'] > 16.9:
            returnDict['BeWindowHeight'] = 16.9
        returnDict['tableHeight'] = M2AbsHeight - self.Table.offset.z
        returnDict['M1Yaw'] = 0.055
        returnDict['M2Yaw'] = 0.08 + kyaw*(m1Pitch-0.12)
        returnDict['M1Mask'] = 1.6 if m1Pitch < 0.15 else 3.6
        for k, v in returnDict.items():
            print(k, '->', v)
        return returnDict

    def repositionToEnergy(self, newEnergy, silent=True):
        def monitor_move(motor, target):
            retryCounts = 0
            time.sleep(0.5)
            while motor.status != 0:
                if hasattr(motor, 'check_status'):
                    motor.check_status()
                if motor.status in [3, 4]:  # ERROR, trying to repeat move
                    if retryCounts < 10:
                        motor.move(target)
                        retryCounts += 1
                    else:
                        print('Error while moving', motor)
                        return -1
                        break
            return 0
        newPositions = self.calculatePositions(newEnergy)
        myInput("Press Enter to continue...")
#        return
        # Step1 changing M1 pitch
        print('M1 mask will be closed!')
        if not silent:
            myInput("Press Enter to continue...")
        self.M1.mask = 0
        if monitor_move(self.M1._maskPV, 0):
            return

        print('M1.Pitch:', self.M1.pitch, '->', newPositions['mirrorPitch'])
        if not silent:
            myInput("Press Enter to continue...")
        self.M1.pitch = newPositions['mirrorPitch']
        if monitor_move(self.M1._pitchPV, newPositions['mirrorPitch']):
            return

        print('M1.Yaw:', self.M1.yaw, '->', newPositions['M1Yaw'])
        if not silent:
            myInput("Press Enter to continue...")
        self.M1.yaw = newPositions['M1Yaw']
        if monitor_move(self.M1._yawPV, newPositions['M1Yaw']):
            return

        print('M1.Bend:', self.M1.bend, '->', newPositions['M1Bend'])
        if not silent:
            myInput("Press Enter to continue...")
        self.M1.bend = newPositions['M1Bend']
        if monitor_move(self.M1._bendPV, newPositions['M1Bend']):
            return

        nsteps = int(np.ceil(np.abs(self.M2.height - newPositions['M2Height'])))
        print('DCM.Height:', self.DCM.height, '->' ,newPositions['monoHeight'])
        print('M2.Height:', self.M2.height, '->' ,newPositions['M2Height'])
        print(nsteps, 'steps to move')
        if not silent:
            myInput("Press Enter to continue...")
        for M2Pos, monoPos in zip(
                np.linspace(self.M2.height, newPositions['M2Height'], nsteps),
                np.linspace(self.DCM.height, newPositions['monoHeight'],
                            nsteps)):
#            print 'DCM.Height to', monoPos, 'M2.Height to', M2Pos
            dtime = datetime.now()
            filename = "{0}_M2_Height_{1:.2f}.txt".format(
                    dtime.strftime("%H-%M-%S"), self.M2.height)
            fobj = open(filename, "w")
            subprocess.call(['/home/bioxas-m/Documents/KeyenceTest/snapshotM2.sh'],
                            stdout=fobj)
            time.sleep(2)
            fobj.close()
            print('DCM.Height:', self.DCM.height, '->', monoPos)
            print('M2.Height:', self.M2.height, '->', M2Pos)
#            myInput("Press Enter to continue...")
            self.DCM.height = monoPos
            self.M2.height = M2Pos
            monitor_move(self.DCM._heightPV, monoPos)
            monitor_move(self.M2._heightPV, M2Pos)
            print('Movement completed')

        print('M2.Pitch:', self.M2.pitch, '->', newPositions['mirrorPitch'])
        if not silent:
            myInput("Press Enter to continue...")
        self.M2.pitch = newPositions['mirrorPitch']
        if monitor_move(self.M2._pitchPV, newPositions['mirrorPitch']):
            return

        print('M2.Yaw:', self.M2.yaw, '->', newPositions['M2Yaw'])
        if not silent:
            myInput("Press Enter to continue...")
        self.M2.yaw = newPositions['M2Yaw']
        if monitor_move(self.M1._yawPV, newPositions['M2Yaw']):
            return

        print('DCM.Energy:', self.DCM.energy, '->', newEnergy)
        if not silent:
            myInput("Press Enter to continue...")
        self.DCM.energy = newEnergy

        print('Be Window Height:', self.BeWindow.height, '->', newPositions['BeWindowHeight'])
        if not silent:
            myInput("Press Enter to continue...")
        self.BeWindow.height = newPositions['BeWindowHeight']
        if monitor_move(self.BeWindow._heightPV,
                        newPositions['BeWindowHeight']):
            return

        print('Table Height:', self.Table.height, '->', newPositions['tableHeight'])
        if not silent:
            myInput("Press Enter to continue...")
        self.Table.height = newPositions['tableHeight']
        if monitor_move(self.Table._heightPV, newPositions['tableHeight']):
            return

        print('M1 mask opening!')
        if not silent:
            myInput("Press Enter to continue...")
        self.M1.mask = newPositions['M1Mask']
        if monitor_move(self.M1._maskPV, newPositions['M1Mask']):
            return

#        print('M1.Bend:', self.M1.bend, '->', newPositions['M1Bend'])
#        if not silent:
#            myInput("Press Enter to continue...")
#        self.M1.bend = newPositions['M1Bend']
#        if monitor_move(self.M1._bendPV, newPositions['M1Bend']):
#            return
#
#
#        print('M2.Bend:', self.M2.bend, '->', newPositions['M2Bend'])
#        if not silent:
#            myInput("Press Enter to continue...")
#        self.M2.bend = newPositions['M2Bend']
#        if monitor_move(self.M2._bendPV, newPositions['M2Bend']):
#            return

        print('Beamline successfully repositioned to', newEnergy, 'eV.')
        print('Do not forget to calibrate energy!')
