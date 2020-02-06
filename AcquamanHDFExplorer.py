# -*- coding: utf-8 -*-
"""
Created on Apr 26. Modified by chernir on Jan 22 2020

@author: konkle
"""
import time
import os
import re
from functools import partial
import json
import itertools
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import LSQUnivariateSpline, interp2d
from scipy.ndimage import gaussian_filter1d
from scipy.integrate import simps
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import h5py
import matplotlib as mpl
from matplotlib.figure import Figure
mpl.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.widgets import SpanSelector
pg.setConfigOptions(antialias=False)
from acquamanDataFormats import dataFormats
from lmfit import Model, models
from matplotlib import pyplot as plt

isTest = True
needFY = True
allAtOnce = True

FRAME_CHUNK_GENERAL = 1
FRAME_CHUNK_MU = 20
FRAME_CHUNK_EXAFS = 50

DETECTOR_NCHANNELS = 1
#DETECTOR_SUM_STR = 'BioXASMainOutboardDetector'
#DETECTOR_PIX_STR = 'BioXASMainOutboardDetectorRawSpectrum{}'
#ROOTPATH = r"X:\bioxas-m\AcquamanMainData\users"
ROOTPATH = ""

FLUOBIN = 10






ch = 12398.419297617678  # c*h[eV*A]
eV2revA = 0.2624682843  # 2m_e(eV)/(hbar(eVs)c(Å/s))^2

LABELS_TR = ['aem', 'current']
LABELS_FY = ['lima', 'roi', 'fy', 'fluo']
SIGNAL_KIND_UNKNOWN, SIGNAL_KIND_TR, SIGNAL_KIND_FY = range(3)

COLOR_LIME = '#BFFF00'
COLOR_ORCHID = '#DA70D6'
COLOR_GRAY = '#888888'
COLOR_DEFAULT = COLOR_GRAY
colorCycleT = ['#ff0000', '#00ff00', '#00ffff', '#ffff00',
               '#0000ff', '#ff00ff', '#ffffff', '#aaaaaa']
colorCycleF = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
               '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
labelStyle = {'color': '#AAA', 'font-size': '10pt'}
PEN_WIDTH = 1.

"""
To pass configs:
door = taurus.Attribute('Balder/Door/01)
door.get_property("my_variable") after
door.put_property("my_variable")
"""


class DataReceiver(QtCore.QObject):
    msgClear = QtCore.pyqtSignal()
    msgMoveable = QtCore.pyqtSignal(list)
    msgCounters = QtCore.pyqtSignal(list)
    msgRepeats = QtCore.pyqtSignal(list)
    msgHintE0 = QtCore.pyqtSignal(list)

    msgRawData = QtCore.pyqtSignal(list)
    msgMuData = QtCore.pyqtSignal(list)

    def __init__(self):
        QtCore.QObject.__init__(self)
        if isTest:
#            self.testTimer = QtCore.QTimer()
#            self.testTimer.setSingleShot(True)
#            self.testTimer.timeout.connect(self.test_random_listener)
#            self.testTimer.timeout.connect(self.test_data_listener)
#            self.testTimer.timeout.connect(self.test_data_listener_FY)
#            self.testTimer.start(1000)
            self.door = None
        else:
            import taurus
            self.door = taurus.Device('Balder/Door/01')
            dp = taurus.Attribute('Balder/Door/01/RecordData')
            dp.addListener(self.listener)

    def listener(self, evt_src, evt_type, evt_value):
        if evt_type != 0:
            return
#        print(evt_src, evt_type, evt_value)
        try:
            dataJSON = json.loads(evt_value.rvalue[1])['data']
            if 'ref_moveables' in dataJSON:
#                print(dataJSON)
                self.msgClear.emit()

                moveable = dataJSON['ref_moveables'][0]
                column_desc = dataJSON['column_desc']
                # column_desc[0]: Pt No
                # column_desc[1]: moveable
                # column_desc[2]: timer
                # column_desc[3:-1]: counters
                # column_desc[-1]: dt
                movMin = column_desc[1]['min_value']
                movMax = column_desc[1]['max_value']
                if 'energy' in moveable:
                    unit = 'eV'
                elif moveable.endswith(('_rol', '_yaw')):  # mirrors
                    unit = 'mrad'
                elif moveable.endswith(('rol', 'yaw')):  # crystals
                    unit = u'µrad'
                elif moveable.endswith(('_x', '_y', '_z')):
                    unit = 'mm'
                else:
                    unit = ''
                self.msgMoveable.emit([moveable, unit, movMin, movMax])
                self.moveable = moveable

                allCounters = list(dataJSON['counters'])
                counters, labels, units = [], [], []
                kinds, locations, colors = [], [], []
                iTR, iFY = 0, 0
                for cd, c in zip(column_desc[3:], allCounters):
                    cond = cd['conditioning'].lower()
                    kind = SIGNAL_KIND_UNKNOWN
                    color = COLOR_DEFAULT
                    for label in LABELS_TR:
                        if label in c or label in cond:
                            kind = SIGNAL_KIND_TR
                            iTR += 1
                            color = colorCycleT[(iTR-1) % len(colorCycleT)]
                            break
                    for label in LABELS_FY:
#                        print(label, cond)
                        if label in c or label in cond:
                            kind = SIGNAL_KIND_FY
                            iFY += 1
                            color = colorCycleF[(iFY-1) % len(colorCycleF)]
#                            print("found fy")
                            break
                    if 'skip' in cond or 'hid' in cond:
                        continue  # can't be right after 'for' for right color
                    isRight = 'right' in cond
                    if isRight and kind == SIGNAL_KIND_TR:
                        loc = 2
                    elif not isRight and kind == SIGNAL_KIND_FY:
                        loc = 3
                    elif isRight and kind == SIGNAL_KIND_FY:
                        loc = 4
                    else:
                        loc = 1
                    unit = u'A' if 'aem' in c else 'counts'
                    label = cd['label']

                    counters.append(c)
                    labels.append(label)
                    units.append(unit)
                    kinds.append(kind)
                    locations.append(loc)
                    colors.append(color)
                self.counters = counters
                self.locations = locations
                self.msgCounters.emit([labels, units, locations, colors])

                serialNo = dataJSON['serialno']
                repeat, repeats = None, None
                try:
                    if self.door is not None:
                        repeat = self.door.get_property('repeat')['repeat']  # a list of str  # noqa
                        if repeat != []:
                            repeat = repeat[0]
                        repeats = self.door.get_property('repeats')['repeats']
                        if repeats != []:
                            repeats = repeats[0]
                except IndexError:
                    pass
                self.msgRepeats.emit([serialNo, repeat, repeats])
                try:
                    hintE0 = self.door.get_property('hintE0')['hintE0']  # a list of str  # noqa
                    if hintE0 != []:
                        self.msgHintE0.emit(
                            [float(hintE0[0]), float(hintE0[1])])
                except:  # noqa
                    pass
            else:
                data = [dataJSON[d] for d in [self.moveable] + self.counters]
                for i, c in enumerate(self.counters):
                    if 'aem' in c:
                        data[i+1] *= 1e-3  # ind 0 is moveable
                self.msgRawData.emit(data)
                if 'energy' in self.moveable:
                    i0 = data[1] + data[2]
                    try:
                        muTData = np.log(i0 / data[3])
                    except:  # noqa
                        muTData = 0.
                    muData = [data[0], muTData]
                    for d, loc in zip(data[1:], self.locations):
                        if loc in [2, 3]:
                            try:
                                muFData = d / i0
                            except:  # noqa
                                muFData = 0.
                            muData.append(muFData)
                    self.msgMuData.emit(muData)  # 1st is T
        except:  # noqa
            pass

    def test_random_listener(self):
        self.msgClear.emit()
#        moveable = 'mot_dummy'
        moveable = 'energy'
        movMin = 6000
        movMax = 7000
        self.msgMoveable.emit([moveable, 'eV', movMin, movMax])

        counters = ['aem01.ch1', 'aem01.ch2', 'aem01.ch3',
                    'ch1_roi1', 'ch2_roi1']
#        counters = ['aem01.ch1', 'aem01.ch2', 'aem01.ch3', 'aem01.ch4']
        # 1: T-left, 2: T-right, 3: F-left, 4: F-right ::
        locations = [1, 1, 2, 3, 3]
#        locations = [1, 1, 2, 2]
        units = [u'A' if 'aem' in c else 'counts' for c in counters]
        colors = colorCycleT[0:3] + colorCycleF[0:2]
        self.msgCounters.emit([counters, units, locations, colors])
        self.msgRepeats.emit([None, 2, 3])

        xs = np.linspace(movMin, movMax, 300)
        for x in xs:
            data = [x]
            data += [np.random.random()*1e-5 for d in counters[:3]]
            data += [1e4 + np.random.random()*1e3 for d in counters[3:]]
            data[3] *= 2
            I0 = data[1] + data[2]
            try:
                muDataT = np.log(I0 / data[3])
            except:  # noqa
                muDataT = 0
            try:
                muDataF1 = data[4] / I0
            except:  # noqa
                muDataF1 = 0
            try:
                muDataF2 = data[5] / I0
            except:  # noqa
                muDataF2 = 0
            self.msgRawData.emit(data)
            if 'energy' in moveable:
                self.msgMuData.emit([x, muDataT, muDataF1, muDataF2])
            time.sleep(0.01)

    def test_data_listener(self):
        fileNames = ['test_data/Cu_004.cur']
#        fileNames = ['test_data/CuO_fly1.cur']
        for iFile, fileName in enumerate(fileNames):
            self.msgClear.emit()
            e, i0, i1 = np.loadtxt(
                fileName, skiprows=1, unpack=True, usecols=(0, 1, 2))
            i0 *= 1e-3
            i1 *= 1e-3
            moveable = 'energy'
            movMin = e.min()
            movMax = e.max()
            self.msgMoveable.emit([moveable, 'eV', movMin, movMax])
            counters = ['i0', 'i1']
            # 1: T-left, 2: T-right, 3: F-left, 4: F-right ::
            locations = [1, 2]
            units = [u'A', u'A']
            colors = colorCycleT[0:2]
            self.msgCounters.emit([counters, units, locations, colors])
            self.msgRepeats.emit([None, iFile+1, len(fileNames)])

            print("started")
            t0 = time.time()
            for eRun, i0Run, i1Run in zip(e, i0, i1):
                data = [eRun, i0Run, i1Run]
                try:
                    muDataT = np.log(data[1]/data[2])
                except:  # noqa
                    muDataT = 0
                self.msgRawData.emit(data)
                self.msgMuData.emit([eRun, muDataT])
                time.sleep(0.002)
            print("finished in {0} s".format(time.time()-t0))

    def test_data_listener_FY(self):
#        fileNames = ['test_data/S90T3_spot6_010.dat']
#        fileNames = [r"G:\hdf5\JG_slac_Ask_exafs_13k_3_1.hdf5",
#                     r"G:\hdf5\JG_slac_Ask_exafs_13k_3_2.hdf5",
#                     r"G:\hdf5\JG_slac_Ask_exafs_13k_3_3.hdf5"]
        fileNames = [r"D:\BioXAS\George\JG_slac_Ask_exafs_13k_1.h5"]
        for iFile, fileName in enumerate(fileNames):
            self.msgClear.emit()
#            e, ch1, ch2, ch3, ch4 = np.loadtxt(
#                fileName, unpack=True, usecols=(1, 3, 4, 5, 6))
            f = h5py.File(fileName, 'r')
            scanGroup = f['scan']
            e = np.array(scanGroup['AxisValues::BioXASEnergyControl'])
#            e = np.array(scanGroup['EnergyFeedback'])  Feedback values not monotonic
            ch1 = np.array(scanGroup['I0Detector_darkCorrected'])
            ch2 = np.array(scanGroup['I1Detector_darkCorrected'])
            ch3 = np.array(scanGroup['I2Detector_darkCorrected'])
            ch4 = np.sum(scanGroup['BioXASMainOutboardDetector'][:,950:1200], axis=1)

            moveable = 'energy'
            movMin = e.min()
            movMax = e.max()
            self.msgMoveable.emit([moveable, 'eV', movMin, movMax])
            counters = ['ch1', 'ch2', 'ch3', 'ch4']
            locations = [1, 1, 1, 3]
            units = [u'counts', u'counts', u'counts', u'counts']
            colors = colorCycleT[0:3] + colorCycleF[0:1]
            self.msgCounters.emit([counters, units, locations, colors])
            self.msgRepeats.emit([None, iFile+1, len(fileNames)])
#            self.msgHintE0.emit([7112-50, 7112+50])
            self.msgHintE0.emit([11867-50, 11867+50])

            print("started")
            t0 = time.time()
            if allAtOnce:
                data = [e, ch1, ch2, ch3, ch4]
                i0 = ch1
                try:
                    muDataT = np.log(i0/ch2)
                except:  # noqa
                    muDataT = 0
                self.msgRawData.emit(data)
                try:
                    muDataF = ch4 / i0
                except:  # noqa
                    muDataF = 0
                self.msgMuData.emit([e, muDataT, muDataF])
            else:
                for eRun, ch1Run, ch2Run, ch3Run, ch4Run in zip(
                        e, ch1, ch2, ch3, ch4):
                    data = [eRun, ch1Run, ch2Run, ch3Run, ch4Run]
                    i0 = ch1Run + ch2Run
                    try:
                        muDataT = np.log(i0/ch3Run)
                    except:  # noqa
                        muDataT = 0
                    self.msgRawData.emit(data)
                    try:
                        muDataF = ch4Run / i0
                    except:  # noqa
                        muDataF = 0
                    self.msgMuData.emit([eRun, muDataT, muDataF])
                    time.sleep(0.002)

            print("finished in {0} s".format(time.time()-t0))

    def hdf_data_listener(self, inplist):
        fileName, edge, e, ch1, ch2, ch3, ch4list = inplist
        self.msgClear.emit()
        moveable = 'energy'
        movMin = e.min()
        movMax = e.max()
        self.msgMoveable.emit([moveable, 'eV', movMin, movMax])        
        counters = ['ch1', 'ch2', 'ch3', 'ch4']
        locations = [1, 1, 1, 3]
        units = [u'counts', u'counts', u'counts', u'counts']
        colors = colorCycleT[0:3] + colorCycleF[0:1]
        self.msgCounters.emit([counters, units, locations, colors])
        self.msgRepeats.emit([None, 1, len(fileName)])
#            self.msgHintE0.emit([7112-50, 7112+50])
        self.msgHintE0.emit([edge-50, edge+50])        

        data = [e, ch1, ch2, ch3, ch4list[0]]
        i0 = ch1
        try:
            muDataT = np.log(i0/ch2)
        except:  # noqa
            muDataT = 0
        self.msgRawData.emit(data)
        try:
            muDataF = ch4list[0] / i0
        except:  # noqa
            muDataF = 0
        self.msgMuData.emit([e, muDataT, muDataF])
        time.sleep(0.002)


class EXAFS_Extractor(QtCore.QObject):
    """With some modifications, largely taken from
    github.com/vadmu/EXAFS_Monitor
    """
    msgReady = QtCore.pyqtSignal(list)

    def __init__(self):
        super(EXAFS_Extractor, self).__init__()
        self.pre1 = -150
        self.pre2 = -50
        self.post1 = 50
        self.r = np.arange(0.0, 6.0, 0.02)
        self.hintE0 = []

    def prepare(self, lst):
        self.e, self.muT, self.muFs = lst

    def victoreenMod(self, e, c, d):
        f = ch / e
        return c*f**-3 + d*f

    def e0(self, e, mu=None, precalcDeriv=None):
        try:
            if precalcDeriv is not None:
                if self.hintE0:
                    cond = (e >= self.hintE0[0]) & (e <= self.hintE0[1])
                    indD = np.argmax(cond)
                    deriv = precalcDeriv[cond]
                    testE = e[cond]
                else:
                    deriv = precalcDeriv
                    testE = e
                    indD = 0
                ind = np.argmax(deriv)
                return ind + indD, testE[ind], deriv[ind]
            else:
                if self.hintE0:
                    cond = (e >= self.hintE0[0]) & (e <= self.hintE0[1])
                    indD = np.argmax(cond)
                    testMu = mu[cond]
                    testE = e[cond]
                else:
                    testMu = mu
                    testE = e
                    indD = 0
                deriv = np.gradient(testMu) / np.gradient(testE)
                deriv = gaussian_filter1d(deriv, 3)
                ind = np.argmax(deriv)
                return ind + indD, testE[ind], deriv[ind]
        except (IndexError, ValueError):
            return

    def pre_edge(self, e, ie0, mu):
        try:
            idx1, idx2 = \
                np.searchsorted(e, [e[ie0]+self.pre1, e[ie0]+self.pre2])
            popt, pcov = \
                curve_fit(self.victoreenMod, e[idx1:idx2+1], mu[idx1:idx2+1])
            return self.victoreenMod(e, popt[0], popt[1])
        except TypeError:
            return

    def k(self, e, mu):
        rese0 = self.e0(self.e, mu)
        if rese0 is None:
            return
        ie0 = rese0[0]
        k = np.sqrt(eV2revA * (e[ie0:]-e[ie0]))
        kmin_idx = np.searchsorted(e, e[ie0]+self.post1) - ie0
        if kmin_idx > len(k) - 2:
            return
        return ie0, k, kmin_idx

    def chi(self, e, ie0, mu, pre, k, kmin_idx, kw=2):
        # k starts from e0
        fun_fit = (mu[ie0:] - pre[ie0:]) * k**kw
        n_knots = int(2.0 * (k.max()-k.min()) / np.pi) + 1
        knots = np.linspace(k[kmin_idx], k[-2], n_knots)
#        print(k, fun_fit, knots)
#        print(np.diff(k))
        spl = LSQUnivariateSpline(k, fun_fit, knots)
        post = spl(k[kmin_idx:]) / k[kmin_idx:]**kw
        edge_step = post[0]
        chi = (mu - pre)[ie0:] - edge_step
        chi[kmin_idx:] += edge_step - post
        chi /= edge_step

# =============================================================================
#         import matplotlib.pylab as plt
# #        plt.plot(e, mu)
# #        plt.plot(e, pre)
# #        plt.plot(e, mu-pre)
# #        plt.axvline(e[ie0], color='g')
# #        plt.plot(e[ie0+kmin_idx:], post)
# #        plt.axhline(edge_step, color='g')
#         plt.plot(e[ie0:], chi)
#         plt.show()
#          # raise here
# =============================================================================

#        return gaussian_filter1d(chi*k*k, 3)
        return chi*k*k

    def hann_win(self, k, kmin, kmax, dk):
        win = np.ones_like(k)
        win[k <= kmin] = 0
        cond = (kmin < k) & (k < kmin+dk)
        win[cond] = 0.5*(1 - np.cos(np.pi*(k[cond]-kmin)/dk))
        cond = (kmax-dk < k) & (k < kmax)
        win[cond] = 0.5*(1 - np.cos(np.pi*(kmax-k[cond])/dk))
        win[k >= kmax] = 0
        return win

    def make_ft(self, k, chi, win):
        ft_re = simps(np.cos(2.0*k*self.r[:, np.newaxis])*chi*win, k)
        ft_im = simps(np.sin(2.0*k*self.r[:, np.newaxis])*chi*win, k)
        ft_mag = (1.0/np.pi * np.abs(ft_re**2 + ft_im**2))**0.5
        return ft_mag

    def calculate(self):

        xax = np.copy(self.e)
        diffTest = np.diff(xax) <= 0
        if any(diffTest):
            xax[np.where(diffTest)] -= 1e-3

        try:
            mu = self.muT
            kres = self.k(xax, mu)
            if kres is None:
                self.terminate()
                return
            ie0, k, kmin_idx = kres

            pre_edge = self.pre_edge(xax, ie0, mu)
            if pre_edge is None:
                self.terminate()
                return
            chi = self.chi(xax, ie0, mu, pre_edge, k, kmin_idx)
            winFT = self.hann_win(k, 2., max(k), 1.0)
            ft = self.make_ft(k, chi, winFT)
            ks = [k]
            chis = [chi]
            fts = [ft]
        except ValueError:
            ks = []
            chis = []
            fts = []
            
        try:
            for mu in self.muFs:
                kres = self.k(xax, mu)
                if kres is None:
                    self.terminate()
                    return
                ie0, k, kmin_idx = kres
#                print(ie0, k, kmin_idx)
                pre_edge = self.pre_edge(xax, ie0, mu)
                if pre_edge is None:
                    self.terminate()
                    return
                pre_edge[:] *= 0
#                print(pre_edge)
                chi = self.chi(xax, ie0, mu, pre_edge, k, kmin_idx)
                winFT = self.hann_win(k, 2., max(k), 1.0)
                ft = self.make_ft(k, chi, winFT)
                ks.append(k)
                chis.append(chi)
                fts.append(ft)
            res = [ks, chis, fts, winFT]
        except (TypeError, IndexError, ValueError):
            raise
            res = []
        self.terminate(res)

    def terminate(self, res=[]):
        self.thread().terminate()
        self.msgReady.emit(res)


def common_substring(sa, sb, isReversed=False):
    """finds the longest common substring of sa and sb"""
    def _iter():
        for a, b in zip(sa[::-1] if isReversed else sa,
                        sb[::-1] if isReversed else sb):
            if a == b:
                yield a
            else:
                return
    res = ''.join(_iter())
    return res[::-1] if isReversed else res


def intervals_extract(iterable):
    iterable = sorted(set(iterable))
    for key, gr in itertools.groupby(
            enumerate(iterable), lambda t: int(t[1])-int(t[0])):
        gr = list(gr)
        yield [gr[0][1], gr[-1][1]]


def make_int_ranges(iterable):
    """examples:

    a = ['2', '3', '4', '5', '7', '8', '9', '11', '15', '16', '17', '18']
    print(make_int_ranges(a))
     -> "[2..5, 7..9, 11, 15..18]"

    a = ['03', '02', '04', '05']
    print(make_int_ranges(a))
     -> "[02..05]"
    """
    ranges = list(intervals_extract(iterable))
    aslist = ["{0[0]}..{0[1]}".format(r) if r[0] < r[1] else
              "{0[0]}".format(r) for r in ranges]
    return "[{}]".format(', '.join(aslist))


class MonitorWidget(QtGui.QWidget):
    mesRequestEXAFS = QtCore.pyqtSignal(list)

    def __init__(self, parent=None):
        super(MonitorWidget, self).__init__(parent)
        self.app = QtCore.QCoreApplication.instance()
#        pg.setConfigOption('background',
#                           self.palette().color(QtGui.QPalette.Background))
#        pg.setConfigOption('background', '#444444')

        self.font = QtGui.QFont()
        self.font.setPixelSize(15)

        self.threadEXAFS = QtCore.QThread()
        self.exafs_extractor = EXAFS_Extractor()
        self.exafs_extractor.moveToThread(self.threadEXAFS)
        self.exafs_extractor.msgReady.connect(self.setEXAFSdata)
        self.threadEXAFS.started.connect(self.exafs_extractor.calculate)
        self.setClear()

        labelScanNo = QtGui.QLabel('Scan #')
        self.scanNo = QtGui.QLabel('not set')
        self.pointsReceived = QtGui.QLabel(': none points received')
        self.coordinates = QtGui.QLabel()
        self.coordinates.setFont(self.font)

        self.repeatsBox = QtGui.QWidget()
        labelRepeatNo = QtGui.QLabel('Repeat')
        self.repeatNo = QtGui.QLabel(r'<i> not set <\i>')
        labelRepeats = QtGui.QLabel('of')
        self.repeats = QtGui.QLabel(r'<i> not set <\i>')
        repeatsLayout = QtGui.QHBoxLayout()
        repeatsLayout.setContentsMargins(0, 0, 0, 0)
        repeatsLayout.addWidget(labelRepeatNo)
        repeatsLayout.addWidget(self.repeatNo)
        repeatsLayout.addWidget(labelRepeats)
        repeatsLayout.addWidget(self.repeats)
        self.repeatsBox.setLayout(repeatsLayout)
        self.repeatsBox.setVisible(False)

        headerLayout = QtGui.QHBoxLayout()
        headerLayout.setContentsMargins(0, 0, 0, 0)
        headerLayout.addWidget(labelScanNo)
        headerLayout.addWidget(self.scanNo)
        headerLayout.addWidget(self.pointsReceived)
        headerLayout.addStretch()
        headerLayout.addWidget(self.coordinates)
        headerLayout.addStretch()
        headerLayout.addWidget(self.repeatsBox)

        self.tabWidget = QtGui.QTabWidget(parent=self)
        self.tabWidget.setStyleSheet("font-weight: bold")

        self.tabRaw = pg.GraphicsLayoutWidget()
        self.plotI = self.tabRaw.addPlot()
        self.plotI.hideButtons()
        self.plotI.vb.setMouseEnabled(x=True, y=False)
#        self.plotI.showGrid(y=True)
        # create right axis
        # From pyqtgraph/examples/MultiplePlotAxes.py
        self.plotIR = pg.ViewBox()
        self.plotIR.setMouseEnabled(x=True, y=False)
        self.plotI.showAxis('right')
        self.plotI.scene().addItem(self.plotIR)
        self.plotI.getAxis('right').linkToView(self.plotIR)
        self.plotIR.setXLink(self.plotI)
        self._updateTwinViews(self.plotI, self.plotIR)
        self.plotI.vb.sigResized.connect(
            partial(self._updateTwinViews, self.plotI, self.plotIR))

        if needFY:
            self.tabRaw.nextRow()
            self.plotFY = self.tabRaw.addPlot()
            self.plotFY.hideButtons()
            self.plotFY.vb.setMouseEnabled(x=True, y=False)
            self.plotFY.setXLink(self.plotI)
            # create right axis
            # From pyqtgraph/examples/MultiplePlotAxes.py
            self.plotFYR = pg.ViewBox()
            self.plotFYR.setMouseEnabled(x=True, y=False)
            self.plotFY.showAxis('right')
            self.plotFY.scene().addItem(self.plotFYR)
            self.plotFY.getAxis('right').linkToView(self.plotFYR)
            self.plotFYR.setXLink(self.plotFY)
            self._updateTwinViews(self.plotFY, self.plotFYR)
            self.plotFY.vb.sigResized.connect(
                partial(self._updateTwinViews, self.plotFY, self.plotFYR))

        self.curvesI = []
        self.curvesIR = []
        self.curvesFY = []
        self.curvesFYR = []
        self.tabWidget.addTab(self.tabRaw, 'Raw signals')

        self.tabMu = pg.GraphicsLayoutWidget()
        self.plotMuT = self.tabMu.addPlot()
        self.plotMuT.hideButtons()
        self.plotMuT.getAxis('bottom').setLabel('energy (eV)', **labelStyle)
        labelStyleL = dict(labelStyle)
        labelStyleL['color'] = COLOR_LIME
        self.plotMuT.getAxis('left').setLabel(
            r'µd <i>in transmission<\i>', **labelStyleL)
        if needFY:
            self.tabMu.nextRow()
            self.plotMuF = self.tabMu.addPlot()
            self.plotMuF.hideButtons()
            self.plotMuF.setXLink(self.plotMuT)
            self.plotMuF.getAxis('bottom').setLabel(
                'energy (eV)', **labelStyle)
            labelMuF = r'µd <i>in FY<\i>'
            self.plotMuF.getAxis('left').setLabel(labelMuF, **labelStyle)
            self.plotMuF.getAxis('right').setLabel(' ', **labelStyle)
            self.plotMuF.getAxis('right').setVisible(True)
        # create right axis (derivative of mu).
        # From pyqtgraph/examples/MultiplePlotAxes.py
        self.plotMuTP = pg.ViewBox()
        self.plotMuT.showAxis('right')
        self.plotMuT.scene().addItem(self.plotMuTP)
        self.plotMuT.getAxis('right').linkToView(self.plotMuTP)
        self.plotMuTP.setXLink(self.plotMuT)
        labelStyleO = dict(labelStyle)
        labelStyleO['color'] = COLOR_ORCHID
        self.plotMuT.getAxis('right').setLabel(u"(µd)'", **labelStyleO)
        self._updateTwinViews(self.plotMuT, self.plotMuTP)
        self.plotMuT.vb.sigResized.connect(
            partial(self._updateTwinViews, self.plotMuT, self.plotMuTP))
        self.curvesMuT = []
        self.curvesMuF = []
        self.curvesMuTP = []
        self.tabWidget.addTab(self.tabMu, u'µd')

        self.tabChi = pg.GraphicsLayoutWidget()
        self.plotChi = self.tabChi.addPlot()
        self.plotChi.hideButtons()
        self.plotChi.getAxis('bottom').setLabel(
                u'k (Å'+u"\u207B"+u'¹)', **labelStyle)
        self.plotChi.getAxis('left').setLabel(
            u'χ·k² (Å'+u"\u207B"+u'²)', **labelStyle)
        self.tabChi.nextRow()
        self.plotFT = self.tabChi.addPlot()
        self.plotFT.hideButtons()
        self.plotFT.getAxis('bottom').setLabel(u'r (Å)', **labelStyle)
        self.plotFT.getAxis('left').setLabel(
            u'|FT[χ·k²]| (Å'+u"\u207B"+u'³)', **labelStyle)
        self.curvesChi = []
        self.curvesFT = []
        self.tabWidget.addTab(self.tabChi, 'EXAFS  ')
        self.cbEXAFSTr = QtGui.QCheckBox('Tr')
        self.cbEXAFSFY = QtGui.QCheckBox('FY')
        self.cbBoxEXAFS = QtGui.QWidget()
        cbLayout = QtGui.QHBoxLayout()
        cbLayout.setContentsMargins(0, 0, 0, 0)
        cbLayout.addWidget(self.cbEXAFSTr)
        cbLayout.addWidget(self.cbEXAFSFY)
        self.cbBoxEXAFS.setLayout(cbLayout)
        self.tabWidget.tabBar().setTabButton(
            self.tabWidget.tabBar().count()-1, QtGui.QTabBar.RightSide,
            self.cbBoxEXAFS)
        self.cbEXAFSTr.stateChanged.connect(self.updatePlots)
        self.cbEXAFSFY.stateChanged.connect(self.updatePlots)
        self.cbEXAFSTr.setChecked(True)
        if needFY:
            self.cbEXAFSFY.setChecked(True)

        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addLayout(headerLayout)
        mainLayout.addWidget(self.tabWidget)

        self.setLayout(mainLayout)
        self.setMinimumWidth(500)
        self.setMinimumHeight(500)
        self.resize(800, 800)

        plots = [self.plotI, self.plotFY, self.plotMuT, self.plotMuF,
                 self.plotChi, self.plotFT] if needFY else\
            [self.plotI, self.plotMuT, self.plotChi, self.plotFT]
        for plot in plots:
            pg.SignalProxy(plot.scene().sigMouseMoved, rateLimit=60,
                           slot=partial(self.mouseMoved, plot))
            plot.scene().sigMouseMoved.connect(partial(self.mouseMoved, plot))

        self.textE0 = pg.TextItem("", anchor=(1.1, 1), color=COLOR_ORCHID)
        self.arrowE0 = pg.ArrowItem(angle=-90, brush='#00FFFF', pen='#00FFFF')
        self.plotMuTP.addItem(self.textE0)
        self.plotMuTP.addItem(self.arrowE0)

        self.tabWidget.setCurrentIndex(0)
        self.tabWidget.currentChanged.connect(self.updatePlots)

        self.dataReceiver = DataReceiver()
        self.dataReceiver.msgClear.connect(self.setClear)
        self.dataReceiver.msgMoveable.connect(self.setMoveable)
        self.dataReceiver.msgCounters.connect(self.setCounters)
        self.dataReceiver.msgRepeats.connect(self.setRepeats)
        self.dataReceiver.msgHintE0.connect(self.setHintE0)
        self.dataReceiver.msgRawData.connect(self.setRawData)
        self.dataReceiver.msgMuData.connect(self.setMuData)

    def mouseMoved(self, plot, evt):
        vb = plot.vb
        try:
            if plot.sceneBoundingRect().contains(evt):
                mousePoint = vb.mapSceneToView(evt)
                self.coordinates.setText('x={0:.7g}, y={1:.7g}'.format(
                                         mousePoint.x(), mousePoint.y()))
        except:  # noqa
            pass

    def _updateTwinViews(self, plotLeft, plotRight):
        # view has resized; update auxiliary views to match
        plotRight.setGeometry(plotLeft.vb.sceneBoundingRect())
        plotRight.linkedViewChanged(plotLeft.vb, plotRight.XAxis)

    def setClear(self):
        self.isFYenabled = False
        self.nPointsReceived = 0
        self.counters, self.locations = [], []
        self.emu = []
        self.muT = []
        self.muFs = []
        self.moveableValue = []
        self.counterValues = []
        self.chi = []
        self.ft = []
        self.chiFs = []
        self.ftFs = []
        self.setHintE0()

    def setMoveable(self, lst):
        self.moveable, unit, movMin, movMax = lst
        self.isXAFSenabled = 'energy' in self.moveable.lower()
        self.tabWidget.setTabEnabled(1, self.isXAFSenabled)
        self.tabWidget.setTabEnabled(2, self.isXAFSenabled)
        if not self.isXAFSenabled:
            self.tabWidget.setCurrentIndex(0)
        label = self.moveable + '{0}'.format(' ('+unit+')' if unit else '')
        # if use 'units' kwarg of axis: pg adds a prefix to it and scales data
        self.plotI.getAxis('bottom').setLabel(label, **labelStyle)
        if needFY:
            self.plotFY.getAxis('bottom').setLabel(label, **labelStyle)

        self.plotI.setRange(xRange=[movMin, movMax])
        self.plotI.vb.setLimits(xMin=movMin, xMax=movMax, yMin=None, yMax=None)
        self.plotIR.setLimits(xMin=movMin, xMax=movMax)
        if needFY:
            self.plotFY.setRange(xRange=[movMin, movMax])
            self.plotFY.vb.setLimits(xMin=movMin, xMax=movMax)
            self.plotFYR.setLimits(xMin=movMin, xMax=movMax)
        if self.isXAFSenabled:
            self.plotMuT.setRange(xRange=[movMin, movMax])
            self.plotMuT.vb.setLimits(xMin=movMin, xMax=movMax)
            self.plotMuTP.setRange(xRange=[movMin, movMax])
            if needFY:
                self.plotMuF.setRange(xRange=[movMin, movMax])
                self.plotMuF.vb.setLimits(xMin=movMin, xMax=movMax)
            r = self.exafs_extractor.r
            self.plotChi.vb.setLimits(xMin=0)
            self.plotFT.setRange(xRange=[r.min(), r.max()])
            self.plotFT.vb.setLimits(xMin=r.min(), xMax=r.max(), yMin=0)

    def setCounters(self, lst):
        """ lst = counters, units, locations, colors """
        self.counters, units, self.locations, colors = lst
        csTL = [c for c, loc in zip(self.counters, self.locations) if loc == 1]
        clTL = [c for c, loc in zip(colors, self.locations) if loc == 1]
        uTL = [u for u, loc in zip(units, self.locations) if loc == 1]
        labelTL = self.makeAxisLabel(csTL)
        self.plotI.getAxis('left').setLabel(
            labelTL, units=uTL[0], **labelStyle)

        csTR = [c for c, loc in zip(self.counters, self.locations) if loc == 2]
        clTR = [c for c, loc in zip(colors, self.locations) if loc == 2]
        uTR = [u for u, loc in zip(units, self.locations) if loc == 2]
        if csTR:
            self.plotI.getAxis('right').setVisible(True)
            labelTR = self.makeAxisLabel(csTR)
            self.plotI.getAxis('right').setLabel(
                labelTR, units=uTR[0], **labelStyle)
        else:
            self.plotI.getAxis('right').setVisible(True)
            self.plotI.getAxis('right').setLabel(units=uTL[0], **labelStyle)
            self.plotIR.setYLink(self.plotI.vb)

        self.isFYenabled = \
            len([l for l in self.locations if l in [3, 4]]) > 0 and needFY
        if self.isFYenabled:
            csFL = [c for c, loc in zip(self.counters, self.locations) if
                    loc == 3]
            clFL = [c for c, loc in zip(colors, self.locations) if loc == 3]
            uFL = [u for u, loc in zip(units, self.locations) if loc == 3]
            labelFL = self.makeAxisLabel(csFL)
            self.plotFY.getAxis('left').setLabel(
                labelFL, units=uFL[0], **labelStyle)

            csFR = [c for c, loc in zip(self.counters, self.locations) if
                    loc == 4]
            clFR = [c for c, loc in zip(colors, self.locations) if loc == 4]
            uFR = [u for u, loc in zip(units, self.locations) if loc == 4]
            if csFR:
                self.plotFY.getAxis('right').setVisible(True)
                labelFR = self.makeAxisLabel(csFR)
                self.plotFY.getAxis('right').setLabel(
                    labelFR, units=uFR[0], **labelStyle)
            else:
                self.plotFY.getAxis('right').setVisible(True)
                self.plotFY.getAxis('right').setLabel(
                    ' ', units=uFL[0], **labelStyle)
                self.plotFYR.setYLink(self.plotFY.vb)
        else:
            csFL, clFL, csFR, clFR = [], [], [], []
            self.cbEXAFSFY.setChecked(False)
            self.cbEXAFSFY.setEnabled(False)

        self.setupFYplots()

        self._curves(self.plotI, self.curvesI, csTL, clTL)
        self._curves(self.plotIR, self.curvesIR, csTR, clTR)
        if self.isFYenabled:
            self._curves(self.plotFY, self.curvesFY, csFL, clFL)
            self._curves(self.plotFYR, self.curvesFYR, csFR, clFR)

        if self.isXAFSenabled:
            if len(self.curvesMuT) < 1:
                curveItem = pg.PlotCurveItem()
                self.plotMuT.addItem(curveItem)
                self.curvesMuT.append(curveItem)
            self.curvesMuT[0].setPen(COLOR_LIME, width=PEN_WIDTH)
            if len(self.curvesMuTP) < 1:
                curveItem = pg.PlotCurveItem()
                self.plotMuTP.addItem(curveItem)
                self.curvesMuTP.append(curveItem)
            self.curvesMuTP[0].setPen(COLOR_ORCHID, width=PEN_WIDTH)
            if self.isFYenabled:
                self._curves(
                    self.plotMuF, self.curvesMuF, csFL+csFR, clFL+clFR)

            self._curves(self.plotChi, self.curvesChi, csFL+csFR, clFL+clFR, 2)
            self._curves(self.plotFT, self.curvesFT, csFL+csFR, clFL+clFR, 1)
            self.curvesChi[0].setPen(COLOR_LIME, width=PEN_WIDTH)
            self.curvesChi[1].setPen(COLOR_GRAY, width=PEN_WIDTH)
            self.curvesFT[0].setPen(COLOR_LIME, width=PEN_WIDTH)

    def _curves(self, plot, curves, counters, colors, dc=0):
        while len(curves) < len(counters)+dc:
            curveItem = pg.PlotCurveItem()
            plot.addItem(curveItem)
            curves.append(curveItem)
        while len(curves) > len(counters)+dc:
            del curves[-1]
            plot.removeItem(curves[-1])
        for curve, color in zip(curves[dc:], colors):
            curve.setPen(color, width=PEN_WIDTH)

    def makeAxisLabel(self, counters):
        try:
            cs = ''
            for c in counters:
                if not cs:
                    cs = c
                    continue
                cs = common_substring(cs, c)
                if not cs:
                    break
            cr = ''
            for c in counters:
                if not cr:
                    cr = c
                    continue
                cr = common_substring(cr, c, True)
                if not cr:
                    break
            csN = len(cs)
            crN = len(cr)
            num = [c[csN:] if crN == 0 else c[csN:-crN] for c in counters]
            label = cs + make_int_ranges(num) + cr
        except:  # noqa
            label = ', '.join([c for c in counters])
        return label

    def setupFYplots(self):
        if needFY:
            self.plotFY.setVisible(self.isFYenabled)
        if self.isFYenabled:
            self.plotI.hideAxis('bottom')
        else:
            self.plotI.showAxis('bottom')

        if not self.isXAFSenabled:
            return

        if needFY:
            self.plotMuF.setVisible(self.isFYenabled)
        if self.isFYenabled:
            self.plotMuT.hideAxis('bottom')
        else:
            self.plotMuT.showAxis('bottom')

    def setRepeats(self, lst):
        try:
            visible = lst[2] > 1
        except:  # noqa
            visible = False

        self.repeatsBox.setVisible(visible)
        if lst[0]:
            self.scanNo.setText(str(lst[0]))
        if lst[1]:
            self.repeatNo.setText(str(lst[1]))
        if lst[2]:
            self.repeats.setText(str(lst[2]))

    def setHintE0(self, lst=[]):
        self.exafs_extractor.hintE0 = lst

    def setRawData(self, lst):
        """lst = [moveable] + counters"""
        incomingSignalShape = np.array(lst[0]).shape
        if len(incomingSignalShape) > 0:
            self.moveableValue = lst[0]
        else:
            self.moveableValue.append(lst[0])

        if not self.counterValues:
            self.counterValues = [[] for i in lst[1:]]

        if len(lst[1:]) != len(self.counterValues):
            raise ValueError("wrong length of data")

        if len(incomingSignalShape) > 0:
            for data, (ic, counterValue) in zip(lst[1:], enumerate(self.counterValues)):
                self.counterValues[ic] = data
        else:            
            for data, counterValue in zip(lst[1:], self.counterValues):
                counterValue.append(data)
            
        self.nPointsReceived += 1
        self.pointsReceived.setText(': {0} point{1}'.format(
            self.nPointsReceived, 's' if self.nPointsReceived > 1 else ''))

        if self.tabWidget.currentIndex() == 0:
            if self.nPointsReceived % FRAME_CHUNK_GENERAL == 0:
                self.updatePlots()
#        if not self.isXAFSenabled:
        self.app.processEvents()

    def setMuData(self, lst):
        """lst = [e, muT] + [muFY]*N"""
        if not self.isXAFSenabled:
            return
        incomingSignalShape = np.array(lst[0]).shape
        if len(incomingSignalShape) > 0:
            self.emu = lst[0]
            self.muT = lst[1]
        else:
            self.emu.append(lst[0])
            self.muT.append(lst[1])


        if self.isFYenabled:
            if not self.muFs:
                self.muFs = [[] for i in lst[2:]]
            if len(lst[2:]) != len(self.muFs):
                raise ValueError("wrong length of mu data")
            if len(incomingSignalShape) > 0:
                for data, (im, muf) in zip(lst[2:], enumerate(self.muFs)):
                    self.muFs[im] = data
            else:
                for data, muf in zip(lst[2:], self.muFs):
                    muf.append(data)

#        if self.tabWidget.currentIndex() == 1:
#            if self.nPointsReceived % FRAME_CHUNK_MU == 0:
#                self.updatePlots()
#        elif self.tabWidget.currentIndex() == 2:
#            if self.nPointsReceived % FRAME_CHUNK_EXAFS == 0:
#                self.updatePlots()
        self.updatePlots()
        self.app.processEvents()

    def setEXAFSdata(self, lst):
        if not lst:
            return
        ks, chis, fts, winFT = lst
        visTr = self.cbEXAFSTr.checkState()
        visFY = self.cbEXAFSFY.checkState()
        self.curvesChi[0].setData(ks[0], chis[0])
        self.curvesChi[0].setVisible(visTr)
        try:
            chimaxTr = chis[0].max() if visTr else 1
            chimaxFY = max([chi.max() for chi in chis[1:]]) if visFY else 1
            chimax = max(chimaxTr, chimaxFY)
        except ValueError:
            chimax = 1
        self.curvesChi[1].setData(ks[-1], winFT/winFT.max()*chimax)
        self.curvesFT[0].setData(self.exafs_extractor.r, fts[0])
        self.curvesFT[0].setVisible(visTr)

        for k, chi, ft, curveChi, curveFT in zip(
                ks[1:], chis[1:], fts[1:],
                self.curvesChi[2:], self.curvesFT[1:]):
            curveChi.setData(k, chi)
            curveFT.setData(self.exafs_extractor.r, ft)
            curveChi.setVisible(visFY)
            curveFT.setVisible(visFY)
        self.app.processEvents()

    def updatePlots(self):
        if self.tabWidget.currentIndex() == 0:  # raw signals
            x = np.array(self.moveableValue)
            d = [np.array(counterValue) for counterValue in self.counterValues]
            if not (len(d) == len(self.locations)):
                return
#                raise ValueError("wrong length of raw data")
            nloc1, nloc2, nloc3, nloc4 = 0, 0, 0, 0
            for a, loc in zip(d, self.locations):
                if loc == 1:
                    self.curvesI[nloc1].setData(x, a)
                    nloc1 += 1
                elif loc == 2:
                    self.curvesIR[nloc2].setData(x, a)
                    nloc2 += 1
                elif loc == 3 and self.isFYenabled:
                    self.curvesFY[nloc3].setData(x, a)
                    nloc3 += 1
                elif loc == 4 and self.isFYenabled:
                    self.curvesFYR[nloc4].setData(x, a)
                    nloc4 += 1
        elif self.tabWidget.currentIndex() == 1:  # mu
            if not self.curvesMuT:
                return
            e = np.array(self.emu)
#            try:
            muT = np.array(self.muT)
            self.curvesMuT[0].setData(e, muT)
            muTP = np.gradient(muT) / np.gradient(e)
            muTP = gaussian_filter1d(muTP, 3)
            self.curvesMuTP[0].setData(e, muTP)
            rese0 = self.exafs_extractor.e0(e, precalcDeriv=muTP)
            if rese0 is not None:
                _, e0, muTP0 = rese0
                self.textE0.setHtml(
                    '<span style="font-size:10pt">{0:.2f}</span>'. format(e0))
                self.textE0.setPos(e0, muTP0)
                self.arrowE0.setPos(e0, muTP0)
            else:
                self.textE0.setPos(0, 0)
                self.arrowE0.setPos(0, 0)
            if self.isFYenabled:
                muFs = [np.array(mu) for mu in self.muFs]
                if not (len(muFs) == len(self.curvesMuF)):
                    raise ValueError("wrong length of raw data")
                for muF, curve in zip(muFs, self.curvesMuF):
                    curve.setData(e, muF)
        elif self.tabWidget.currentIndex() == 2:  # EXAFS
            e = np.array(self.emu)
            muT = np.array(self.muT)
            if self.isFYenabled:
                muFs = [np.array(mu) for mu in self.muFs]
                if not (len(muFs) == len(self.curvesMuF)):
                    raise ValueError("wrong length of raw data")
            else:
                muFs = []
            self.exafs_extractor.prepare([e, muT, muFs])
            self.threadEXAFS.start()

class Ge32Explorer(QtGui.QMainWindow):
    pixDataReady = QtCore.pyqtSignal(list)

    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("Acquaman EXAFS Explorer")

        self.file_menu = QtGui.QMenu('&File', self)
        self.file_menu.addAction('&Quit', self.fileQuit,
                                 QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.help_menu = QtGui.QMenu('&Help', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)
        self.exafs_monitor = MonitorWidget()
        self.exafs_monitor.show()
#        self.help_menu.addAction('&About', self.about)

        self.main_widget = QtGui.QWidget(self)
        self.DETECTOR_PIX_STR = ""
        self.DETECTOR_SUM_STR = ""

        hl = QtGui.QHBoxLayout(self.main_widget)



        canvasSplitter = QtGui.QSplitter()
        canvasSplitter.setChildrenCollapsible(False)
        canvasSplitter.setOrientation(QtCore.Qt.Horizontal)
        hl.addWidget(canvasSplitter)
        
        fileNav = QtGui.QTreeView()
        fileModel = QtGui.QDirModel()
#        fileModel = QtGui.QFileSystemModel()
#        fileModel.setRootPath(ROOTPATH)
#        fileNav.setRootIndex(fileModel.index(ROOTPATH))
        fileNav.setSortingEnabled(True)
        fileNav.header().sortIndicatorChanged.connect(fileModel.sort)
        fileNav.sortByColumn(0, QtCore.Qt.AscendingOrder)
        fileNav.setModel(fileModel)
        fileNav.doubleClicked.connect(self.load_data)

        lw = QtGui.QWidget(self.main_widget)
        l = QtGui.QVBoxLayout()
        lw.setLayout(l)


        self.slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.slider.setRange(0, DETECTOR_NCHANNELS)
        self.slider.setTickInterval(1)
        try:
            self.slider.setTickPosition(QtGui.QSlider.TicksAbove)
        except:
            print("slider: unknown settings")
        
        pixPanel = QtGui.QGroupBox(self.main_widget)
        pixPanel.setTitle('MCA Detector channels')
        self.pixLayout = QtGui.QGridLayout(self.main_widget)
        self.pixList = []
        for ipix in range(32):
            self.add_pix_cb(ipix)
        for ipix in range(31):
            self.pixList[ipix+1].setVisible(False)
        pixPanel.setLayout(self.pixLayout)
        pixPanel.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)

        self.mplFig = Figure()
        self.mplAx = self.mplFig.add_subplot(111)
        paletteWidget = FigureCanvas(self.mplFig)
        self.load_data(None)

        self.im = self.mplAx.imshow(self.sumChannel.T, cmap='jet',
                                    extent=(self.eaxis[0], self.eaxis[-1],
                                            0, 4095*FLUOBIN),
                                    aspect='auto',
                                    interpolation='none',
                                    origin='lower')

        self.slider.valueChanged.connect(self.show_frame)
        self.addToolBar(NavigationToolbar(paletteWidget, self))
        self.pixDataReady.connect(self.exafs_monitor.dataReceiver.hdf_data_listener)
        self.span = SpanSelector(self.mplAx, self.onROIselect, 'vertical',
                                 useblit=True, rectprops=dict(alpha=0.3,
                                                              facecolor='white'),
                                 button=1)
        scaleValidator = QtGui.QIntValidator()
        scaleValidator.setRange(0, DETECTOR_NCHANNELS)
        self.pixEdit = QtGui.QLineEdit('0')
        self.pixEdit.editingFinished.connect(lambda: self.slider.setValue(
                int(self.pixEdit.text())))
        self.pixEdit.setMaximumWidth(50)
        self.pixEdit.setValidator(scaleValidator)
        hlayout = QtGui.QHBoxLayout()
        hlayout.addWidget(self.slider)
        hlayout.addWidget(self.pixEdit)
        l.addLayout(hlayout)
        
        infoWidget = QtGui.QGroupBox()
        vlayout = QtGui.QVBoxLayout()

        
#        roiTabs = QtGui.QTabWidget()
#        roiPanel = QtGui.QWidget()
        roiPanel = QtGui.QGroupBox('ROI')
        panelLayout = QtGui.QVBoxLayout()
        roiScroll = QtGui.QScrollArea()
        roiWidget = QtGui.QWidget()
        self.roiButtonGroup = QtGui.QButtonGroup()
        self.roiButtonGroup.buttonClicked.connect(self.updateROIs)
        roiLayout = QtGui.QVBoxLayout()        
        self.roiPanels = []
        for chan in range(32):
            roiChan = ROIPanel(self.roiButtonGroup, chan, 0, 0, 0)
            roiChan.roimaxEdit.editingFinished.connect(lambda: self.onROIselect(None, None))
            roiChan.roiminEdit.editingFinished.connect(lambda: self.onROIselect(None, None))
            roiChan.elineEdit.editingFinished.connect(
                    lambda: self.plot_lines(None, 1))
            roiLayout.addWidget(roiChan)
            if chan > 0:
                roiChan.setVisible(False)
            self.roiPanels.append(roiChan)
        roiLayout.addStretch(1)
        self.roiButtonGroup.button(0).setChecked(True)
        self.currentROIPanel = self.roiPanels[0]
        roiWidget.setLayout(roiLayout)
        roiScroll.setWidget(roiWidget)
        roiScroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        roiScroll.setWidgetResizable(True)
        panelLayout.addWidget(roiScroll)
        roiPanel.setLayout(panelLayout)
        vlayout.addWidget(roiPanel)
#        roiTabs.addTab(roiPanel, 'ROI')

#        pfPanel = QtGui.QWidget()
        pfPanel = QtGui.QGroupBox('Peak fitting')
        pfLayout = QtGui.QVBoxLayout()
        pfCB = QtGui.QCheckBox('Enabled (really slow)')
        pfCB.setCheckState(0)
        self.PFStatus = False
        pfCB.stateChanged.connect(self.switchPF_ROI)
        pfLayout.addWidget(pfCB)
        self.PFModel = 'Gaussian'
        self.pfButtonGroup = QtGui.QButtonGroup()
        for iFS, funcStr in enumerate(['Gaussian', 'Lorentzian', 'Voigt']):
            funcRB = QtGui.QRadioButton(funcStr)
            self.pfButtonGroup.addButton(funcRB, iFS)
            if iFS == 0:
                funcRB.setChecked(True)
            pfLayout.addWidget(funcRB)

        self.pfButtonGroup.buttonClicked['int'].connect(self.updatePFs)            
        pfPanel.setLayout(pfLayout)
#        roiTabs.addTab(pfPanel,'Peak fitting')
#        vlayout.addWidget(roiTabs)
        vlayout.addWidget(pfPanel)

        linesWidget = QtGui.QGroupBox("Lines")
        hlayout = QtGui.QVBoxLayout()
        label = QtGui.QLabel('Edge position, eV')
        self.edgeEdit = QtGui.QLineEdit('0')
        self.edgeEdit.editingFinished.connect(lambda: self.plot_lines(
                float(self.edgeEdit.text()), 0))
        hlayout.addWidget(label)
        hlayout.addWidget(self.edgeEdit)
        linesWidget.setLayout(hlayout)
        vlayout.addWidget(linesWidget)

        plotPropWidget = QtGui.QGroupBox("Energy scale for 2D")
        hlayout = QtGui.QVBoxLayout()
        buttonGroup = QtGui.QButtonGroup()
        button1 = QtGui.QRadioButton("eV (real scale)")
        button1.setChecked(True)
        self.mapEax = 0
        button2 = QtGui.QRadioButton("points (faster, wrong scale)")
        buttonGroup.addButton(button1)
        buttonGroup.addButton(button2)
        button1.clicked.connect(self.updateEnergyAxis)
        button2.clicked.connect(self.updateEnergyAxis)
        hlayout.addWidget(button1)
        hlayout.addWidget(button2)
        plotPropWidget.setLayout(hlayout)
        vlayout.addWidget(plotPropWidget)

#        vlayout.addStretch(1)

        exportGroup = QtGui.QGroupBox("Export")
        agLayout = QtGui.QVBoxLayout()
        self.fileNameEdit = QtGui.QLineEdit("No data")
        self.exAllCB = QtGui.QCheckBox("Export all enabled channels")
        exportButton = QtGui.QPushButton('Export to ASCII')
        exportButton.clicked.connect(self.exportToASCII)
        agLayout.addWidget(self.fileNameEdit)
        agLayout.addWidget(self.exAllCB)
        agLayout.addWidget(exportButton)
        exportGroup.setLayout(agLayout)
        vlayout.addWidget(exportGroup)
        infoWidget.setLayout(vlayout)
        infoWidget.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Expanding)
#        self.roiCountStr = "Total ROI counts for channel {0}: {1}"
#        self.totalROIcount = QtGui.QLabel(self.roiCountStr.format(
#                'SUM', np.sum(self.totalSum[:, self.roi[0]:self.roi[-1]])))
#        vlayout.addWidget(self.totalROICount)
#        self.totalCountStr = "Total counts for channel {0}: {1}"
#        self.totalCount = QtGui.QLabel(self.totalCountStr.format(
#                'SUM', np.sum(self.totalSum)))
#        vlayout.addWidget(self.totalCount)
        hlayout = QtGui.QHBoxLayout()
        hsplit = QtGui.QSplitter()
        hsplit.setChildrenCollapsible(False)
        hsplit.setOrientation(QtCore.Qt.Horizontal)
        hsplit.setLayout(hlayout)
        hlayout.addWidget(paletteWidget)
        hlayout.addWidget(infoWidget)
        l.addWidget(pixPanel)
        l.addWidget(hsplit)
#        l.addStretch()
#        l.addWidget(dc)
        canvasSplitter.addWidget(fileNav)
        canvasSplitter.addWidget(lw)
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

    def updatePFs(self, buttonID):
        self.PFModel = buttonID
        if self.PFStatus:
            self.send_to_monitor()

    def switchPF_ROI(self, state):
        self.PFStatus = True if state else False
        self.send_to_monitor()

    def updateROIs(self, buttonId):
        try:
            sender = self.sender()
            panelId = int(sender.checkedId())
        except:
            panelId = buttonId
        self.currentROIPanel = self.roiPanels[panelId]
        ymin = np.float32(self.currentROIPanel.roiminEdit.text())
        ymax = np.float32(self.currentROIPanel.roimaxEdit.text())
        elEnergy = np.float32(self.currentROIPanel.elineEdit.text())
        self.onROIselect(ymin, ymax)
        self.plot_lines(elEnergy, 1)

    def updateEnergyAxis(self, state):
        sender = self.sender()
        if str(sender.text()).startswith("eV"):
            self.mapEax = 0
        else:
            self.mapEax = 1
        self.show_frame(self.slider.value())

    def exportToASCII(self):
        header = """XDI/1.0\nEndstation\nBioXAS Beamline - Main Endstation\n\n"""
        header += "ROI: {0} to {1} eV\n\n".format(*self.roi)
        exportName = os.path.join(os.path.dirname(
                self.fileNameStr), self.fileNameEdit.text())
        outData = []
        for dsName in ['eaxis', 'dwellTime', 'ch1', 'ch2', 'ch3', 'totalSum']:
            header += "{} ".format(dsName)
            outData.append(np.sum(getattr(self, dsName)[:, int(self.roi[0]*0.1):int(self.roi[-1]*0.1)], axis=1) if dsName == 'totalSum' else getattr(self, dsName))

        if self.exAllCB.checkState():
            for pixel in self.pixList:
                if pixel.checkState():
                    pixStr = self.DETECTOR_PIX_STR.format(pixel.text())
                    header += pixStr + " "
                    outData.append(np.sum(np.array(
                            self.scanGroup[pixStr])[:, int(self.roi[0]*0.1):int(self.roi[-1]*0.1)], axis=1))
        np.savetxt(exportName,
                   np.array(outData).T,
                   fmt='%.10e',
                   header=header)
        print("Exported to", "{}".format(exportName))        

    def onROIselect(self, ymin, ymax):

        if ymin is not None:
            self.currentROIPanel.roiminEdit.setText("{:.2f}".format(ymin))
            self.currentROIPanel.roimaxEdit.setText("{:.2f}".format(ymax))
        else:
            ymin = float(self.currentROIPanel.roiminEdit.text())
            ymax = float(self.currentROIPanel.roimaxEdit.text())

        try:
            self.roispan.remove()
        except:
            pass
        self.roispan = self.mplAx.axhspan(ymin, ymax, alpha=0.3, color='white')
        self.roi = [ymin, ymax]

        self.send_to_monitor()

    def recalculate_sum(self, state):
        sender = self.sender()
        pixNum = int(sender.text())
        modifier = 1. if state else -1.
        self.totalSum += modifier *\
            np.array(self.scanGroup[self.DETECTOR_PIX_STR.format(
                    pixNum)])
        self.send_to_monitor()
        self.sumChannel = self.totalSum / self.ch1arr
#        self.sumChannel /= np.max(self.sumChannel)
        if int(self.pixEdit.text()) == 0:
            self.show_frame(0)

    def send_to_monitor(self):
        outF = self.fit_peaks() if self.PFStatus else np.sum(
                self.totalSum[:, int(self.roi[0]*0.1):int(self.roi[-1]*0.1)],
                axis=1)

        outMsg = [self.sourceName, float(self.edgeEdit.text()), self.eaxis, self.ch1, self.ch2, self.ch3,
                  [outF]]
        self.pixDataReady.emit(outMsg)

    
    def fit_peaks(self):
#        plt.figure()
        if self.PFModel == 0:
            peak1 = models.GaussianModel(prefix='p1_')
        elif self.PFModel == 1:
            peak1 = models.LorentzianModel(prefix='p1_')
        else:
            peak1 = models.VoigtModel(prefix='p1_')

        params1 = peak1.make_params()
        totalFit = np.zeros(self.totalSum.shape[0])
#        lastfit = None
#        initfit = None
        print("Started peak fitting for", len(totalFit), "points")
        for ie in range(len(totalFit)):
            xrf = self.totalSum[ie, int(self.roi[0]*0.1):int(self.roi[-1]*0.1)]
            datax = np.linspace(0, len(xrf)-1, len(xrf))
            center = (self.roi[-1]*0.1 - self.roi[0]*0.1) * 0.5
            params1['p1_center'].set(value=center, min=center-30, max=center+30)
            params1['p1_sigma'].set(value=8, min=0.1, max=100)
            xrfmax = np.max(xrf)
            params1['p1_amplitude'].set(value=xrfmax, min=xrfmax*0.01, max=xrfmax*100)
            output = peak1.fit(xrf, params1, x=datax)
            peak_integral = output.best_values['p1_amplitude']
#            lastparams = params1
#            lastbest = output.best_values
#            print(self.eaxis[ie], peak_integral)
            totalFit[ie] = peak_integral
#            lastfit = output.best_fit
#            initfit = output.init_fit
#        print(lastparams)
#        print(lastbest)
#        plt.plot(datax, xrf, 'b+')
#        plt.plot(datax, lastfit, 'r-')
#        plt.plot(datax, initfit, 'g-')
#        plt.show()
        print("Done")
        return totalFit

    def load_data(self, index):
        self.beamline = None
        if index is None:
            self.eaxis = np.linspace(0, 100, 101)
            self.dwellTime = np.zeros_like(self.eaxis)
            self.ch1 = np.ones_like(self.eaxis)
            self.ch2 = np.zeros_like(self.eaxis)
            self.ch3 = np.zeros_like(self.eaxis)
            ch4 = np.zeros((len(self.eaxis), 1024))
            self.fileNameStr = "No Data"
            self.ch1arr = np.zeros_like(ch4)
            self.sumChannel = np.zeros_like(ch4)
            self.sumChannel[0, 0] = 1
            self.sourceName = "No Data"
            self.f = None
            self.scanGroup = []

        else:
            if not index.model().isDir(index):
                fileName = index.model().filePath(index)
            self.fileNameStr = fileName  # .split("\\")[-1]
            self.fileNameEdit.setText("{}.dat".format(
                    os.path.basename(self.fileNameStr)))
            self.sourceName = self.fileNameStr
            try:
                self.f = h5py.File(fileName, 'r')
                self.scanGroup = self.f['scan']  # BioXAS scan group name

                
                if 'BioXASMainInboardDetector' in self.scanGroup.keys() or\
                        'BioXASMainOutboardDetector' in self.scanGroup.keys():
                    self.beamline = 'BioXAS-Main'
                elif 'Ge32Element' in self.scanGroup.keys():
                    self.beamline = 'BioXAS-Side'
                elif 'KETEK' in self.scanGroup.keys():
                    self.beamline = 'IDEAS'
                elif 'Ge13El' in self.scanGroup.keys():
                    self.beamline = 'VESPERS'
                else:
                    print("Could not identify the beamline based on mca detector name")
                
                if self.beamline is None:
                    if 'beamline_configuration' in self.f.keys():
                        beamCfg = self.f['beamline_configuration']
                        if 'BioXASMainInboardGeDetectorStage' in beamCfg.keys():
                            self.beamline = 'BioXAS-Main'
                        elif 'BioXASSideGeDetectorStage' in beamCfg.keys():
                            self.beamline = 'BioXAS-Side'
                        else:
                            print("Could not identify the beamline based on beamline_configuration")                         
                            
                if self.beamline is None:
                    if 'I0Detector' in self.scanGroup.keys():
                        self.beamline = 'BioXAS-Main'  #  It can be Side as well, but the dataset names would be the same
                    elif 'I_0' in self.scanGroup.keys() or\
                            'Reference' in self.scanGroup.keys():
                        self.beamline = 'IDEAS'
                    elif 'MiniIonChamber' in self.scanGroup.keys() or\
                            'PostIonChamber' in self.scanGroup.keys() or\
                            'SplitIonChamber' in self.scanGroup.keys():
                        self.beamline = 'VESPERS'
                    else:
                        print("CAN'T IDENTIFY THE FILE STRUCTURE")  # TODO: data selection dialog
                    
                if self.beamline is None:
                    self.load_data(None)
                    return

                axisnames = re.findall(r'AxisValues.*?\'', (str(self.scanGroup.keys())))
                if len(axisnames) == 1:
                    axisname = axisnames[0].strip(('\' ,'))
                    self.eaxis = np.array(self.scanGroup[axisname])
                elif len(axisnames) == 2:
                    print("Two-dimensional scans not supported")
                else:
                    print("CAN'T IDENTIFY THE FILE STRUCTURE")  # TODO: data selection dialog

                self.dwellTime = np.array(self.scanGroup[dataFormats[self.beamline]['dwellTime']])
                self.ch1 = np.array(self.scanGroup[dataFormats[self.beamline]['I0']])
                self.ch2 = np.array(self.scanGroup[dataFormats[self.beamline]['I1']])
                self.ch3 = np.array(self.scanGroup[dataFormats[self.beamline]['I2']])

                if isinstance(dataFormats[self.beamline]['mcaTotal'], tuple):
                    self.DETECTOR_SUM_STR = dataFormats[self.beamline]['mcaTotal'][0]
                    self.DETECTOR_PIX_STR = dataFormats[self.beamline]['mcaSingle'][0]
                else:
                    self.DETECTOR_SUM_STR = dataFormats[self.beamline]['mcaTotal']
                    self.DETECTOR_PIX_STR = dataFormats[self.beamline]['mcaSingle'] if\
                        'mcaSingle' in dataFormats[self.beamline].keys() else self.DETECTOR_SUM_STR

                ch4 = np.array(self.scanGroup[self.DETECTOR_SUM_STR])
                fluolen = ch4.shape[1]
                extents = (self.eaxis[0], self.eaxis[-1], 0, (fluolen-1)*10)
#                print(extents)
                self.im.set_extent(extents)
            except:
                raise
                self.load_data(None)
                return
                
        self.roi = [0, np.max(ch4.shape[1])*FLUOBIN]
        self.totalSum = np.zeros_like(ch4)

        if index is not None:
            for ichan in range(32):
                self.pixList[ichan].setVisible(False)

            for ichan in range(32):
                pixStr = self.DETECTOR_PIX_STR.format(ichan+1)
                if pixStr in self.scanGroup.keys():
                    if self.pixList[ichan].checkState():  # and ichan not in [4,5,6,7,15]:
                        self.totalSum += np.array(self.scanGroup[pixStr])
#                    self.add_pix_cb(ichan)
                    self.pixList[ichan].setVisible(True)
                    if self.DETECTOR_PIX_STR == self.DETECTOR_SUM_STR:
                        break
                else:
                    break

            self.slider.setRange(0, ichan)
            self.pixEdit.validator().setRange(0, ichan)
            self.ch1arr = (self.ch1[:, np.newaxis]*np.ones_like(ch4))
            self.sumChannel = self.totalSum / self.ch1arr
#            self.sumChannel /= np.max(self.sumChannel)
            currentPix = int(self.pixEdit.text())
            self.show_frame(currentPix)

            try:
                if 'emission_lines' in self.f.keys():
                    lastroi = 0
                    for iroi, (roiname, roivalue) in enumerate(self.f['emission_lines'].items()):
                        roimin, roimax, roienergy =\
                            np.float32(roivalue['lowerBound']),\
                            np.float32(roivalue['upperBound']),\
                            np.float32(roivalue['energy'])
                        self.roiPanels[iroi].update_rois(roimin, roimax)
                        self.roiPanels[iroi].update_energy(roienergy)
                        self.roiPanels[iroi].roiRButton.setText(str(roiname))
                        self.roiPanels[iroi].setVisible(True)
                        lastroi = int(iroi)+1
                    for ichan in range(32-lastroi):
                        self.roiPanels[ichan+lastroi].setVisible(False)                        

                    self.currentROIPanel = self.roiPanels[0]
                    self.currentROIPanel.roiRButton.setChecked(True)

                    if lastroi > 0:
                        self.updateROIs(0)
            except:
                raise
                print("Can't load ROIs")

            edgePos = min(self.eaxis[-1], self.eaxis[0]+200)

            try:
                if 'scan_regions' in self.f.keys():
                    if len(self.f['scan_regions']) > 0:
                        for regionv in self.f['scan_regions'].values():
                            if 'edgeEnergy' in regionv:
                                edgePos = np.float32(regionv['edgeEnergy'])
                            break
            except:
                print("Can't load edge position")
            self.send_to_monitor()
            self.edgeEdit.setText("{:.2f}".format(edgePos))
            self.plot_lines(edgePos, 0)

    def add_pix_cb(self, ipix):
            row = ipix // 8
            col = ipix - row * 8
            pix = QtGui.QCheckBox(str(ipix+1))
            pix.setChecked(True)
            pix.stateChanged.connect(self.recalculate_sum)
            self.pixLayout.addWidget(pix, row, col)
            self.pixList.append(pix)

    def plot_lines(self, energy, mode):
        if mode:
            try:
                self.emissionLinePlot.remove()
            except:
                pass
            if energy is None:
                energy = np.float32(self.currentROIPanel.elineEdit.text())
            self.emissionLinePlot = self.mplAx.axhline(
                    energy, color='red', linewidth=1)
        else:
            try:
                self.absorptionLinePlot.remove()
            except:
                pass
            self.absorptionLine = energy
            self.absorptionLinePlot = self.mplAx.axvline(
                    energy, color='yellow', linewidth=1)
        self.send_to_monitor()
        self.mplFig.canvas.draw()
        self.mplFig.canvas.blit()

    def show_frame(self, ichan):
        self.pixEdit.setText(str(ichan))
        pixStr = self.DETECTOR_PIX_STR.format(ichan) if ichan > 0 else\
            self.DETECTOR_SUM_STR
        if not pixStr in self.scanGroup:
            return
        fluolen = self.totalSum.shape[1]
        fluoax = np.linspace(0, (fluolen-1)*FLUOBIN, fluolen)
        eax = np.linspace(self.eaxis[0], self.eaxis[-1],
                          int(self.eaxis[-1] - self.eaxis[0]))
        data = np.array(self.scanGroup[pixStr])/self.ch1arr if ichan > 0 else\
            self.sumChannel
#        data = self.ch4
        plotArray = data.T / np.max(data)

        #  It's a workaround for acquaman bug/feature, 
        #  when the next region starts at the same energy
        #  where the previous region ends
        xax = np.copy(self.eaxis)
        diffTest = np.diff(xax) <= 0
        if any(diffTest):
            xax[np.where(diffTest)] -= 1e-3
        #  End fix

        if self.mapEax == 0:
            f = interp2d(xax, fluoax, plotArray)
            interpData = f(eax, fluoax)
        else:
            interpData = plotArray

        self.im.set_data(interpData)
        self.im.set_cmap('jet')

        print(pixStr)
        self.mplAx.set_title(pixStr)
        self.mplFig.canvas.draw()
        self.mplFig.canvas.blit()

    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit() 


class ROIPanel(QtGui.QGroupBox):
    lineEditUpdated = QtCore.pyqtSignal(int, int)
    
    def __init__(self, bGr, chan, e0, emin, emax):
        QtGui.QGroupBox.__init__(self)
#        roiWidget = QtGui.QGroupBox()
        hlayout = QtGui.QGridLayout()
       
        roiRButton = QtGui.QRadioButton(str(chan+1))
#        roiRButton.clicked.connect(self.updateROIs)
        bGr.addButton(roiRButton, chan)
        self.roiRButton = roiRButton
        
        hlayout.addWidget(roiRButton, 0, 0, 2, 1)

        scaleValidator = QtGui.QDoubleValidator()
        scaleValidator.setRange(0, 4095*FLUOBIN, 3)
        
        roiminEdit = QtGui.QLineEdit(str(emin))
        roiminEdit.setValidator(scaleValidator)
        roiminEdit.setMaximumWidth(75)
        roimaxEdit = QtGui.QLineEdit(str(emax))
        roimaxEdit.setValidator(scaleValidator)
        roimaxEdit.setMaximumWidth(75)
        
       
        label = QtGui.QLabel('Energy')
        elineEdit = QtGui.QLineEdit(str(e0))
        elineEdit.setValidator(scaleValidator)
#        elineEdit.setMaximumWidth(75)
        hlayout.addWidget(label, 0, 1)
        hlayout.addWidget(elineEdit, 0, 2, 1, 3)
        label = QtGui.QLabel('Min')
        hlayout.addWidget(label, 1, 1)
        hlayout.addWidget(roiminEdit, 1, 2)
#        hlayout.addStretch(1)
        label = QtGui.QLabel('Max')
        hlayout.addWidget(label, 1, 3)
        hlayout.addWidget(roimaxEdit, 1, 4)

        self.roiminEdit = roiminEdit
        self.roimaxEdit = roimaxEdit
        self.elineEdit = elineEdit

        self.setLayout(hlayout)        

    def update_rois(self, roimin, roimax):
        self.roiminEdit.setText("{:.2f}".format(roimin))
        self.roimaxEdit.setText("{:.2f}".format(roimax))        

    def update_energy(self, energy):
        self.elineEdit.setText("{:.2f}".format(energy))             


if __name__ == '__main__':
    import sys
    app = QtGui.QApplication([])
#    mw = MonitorWidget()
    mw = Ge32Explorer()
#    mw.setStyleSheet("background-color: #000")
#    mw.setWindowTitle("Acquaman EXAFS Explorer")
    mw.show()
    sys.exit(app.exec_())
