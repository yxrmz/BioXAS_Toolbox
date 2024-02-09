# -*- coding: utf-8 -*-
"""
Created on Apr 26

@author: konkle
"""
import time
from functools import partial
import json
import itertools
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import LSQUnivariateSpline
from scipy.ndimage import gaussian_filter1d
from scipy.integrate import simps
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import pickle
import zmq

ZMQ_HOST = "opi1607-904"  # Main
#ZMQ_HOST = "opi1607-1004"  # Side
ZMQ_PORT = "5567"

pg.setConfigOptions(antialias=False)

isTest = True
needFY = True

FRAME_CHUNK_GENERAL = 1
FRAME_CHUNK_MU = 20
FRAME_CHUNK_EXAFS = 50

ch = 12398.419297617678  # c*h[eV*A]
eV2revA = 0.2624682843  # 2m_e(eV)/(hbar(eVs)c(Å/s))^2

# SIDE
#LABELS_TR = ['Scaler SIS3820_channels_ch16_feedback',
#             'Scaler SIS3820_channels_ch17_feedback',
#             'Scaler SIS3820_channels_ch18_feedback']
#LABELS_FY = ['Scaler SIS3820_channels_ch19_feedback']
#SCALER_ALIAS = {'Scaler SIS3820_channels_ch16_feedback': "IC 0",
#                'Scaler SIS3820_channels_ch17_feedback': "IC 1",
#                'Scaler SIS3820_channels_ch18_feedback': "IC 2",
#                'Scaler SIS3820_channels_ch19_feedback': "IC PIPS"}

# MAIN
LABELS_TR = ['Scaler_SIS3820_channels_ch16_feedback',
             'Scaler_SIS3820_channels_ch17_feedback',
             'Scaler_SIS3820_channels_ch18_feedback']
LABELS_FY = ['Scaler_SIS3820_channels_ch19_feedback']
SCALER_ALIAS = {'Scaler_SIS3820_channels_ch16_feedback': "IC 0",
                'Scaler_SIS3820_channels_ch17_feedback': "IC 1",
                'Scaler_SIS3820_channels_ch18_feedback': "IC 2",
                'Scaler_SIS3820_channels_ch19_feedback': "IC PIPS"}

SIGNAL_KIND_UNKNOWN, SIGNAL_KIND_TR, SIGNAL_KIND_FY = range(3)  # (0, 1, 2)

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

motor_dict = {'Energy': 'Energy_ev_feedback',
              'DBHR Pitch': 'DBHR Pitch'}

motor_monitors = {'Energy': 'Energy_ev_sync_feedback',
                  'DBHR Pitch': 'DBHR Pitch'}

"""
To pass configs:
door = taurus.Attribute('Balder/Door/01)
door.get_property("my_variable") after
door.put_property("my_variable")
"""


#class ZeroMQ_Listener(QtCore.QObject):
#    message = QtCore.pyqtSignal(list)
#    def __init__(self):
#        port = 5567
#        QtCore.QObject.__init__(self)
#        self.context = zmq.Context()
#        self.socket = self.context.socket(zmq.SUB)
#        self.socket.connect(f"tcp://127.0.0.1:{port}")
#        self.socket.setsockopt(zmq.SUBSCRIBE, "")
#        self.running = True
#
#    def loop(self):
#        buf = []
#        t0 = time()
#        while self.running:
#            try:
#                message = self.socket.recv()
#                prefix, name, doc = message.split(b' ', 2)
#                #document = json.loads(message)
##                print("Received Document:", prefix, name, pickle.loads(doc))
#
##                string = self.socket.recv()
##                buf.append(string)
##                if time() - t0 > 1./FPS_LIMIT:
##                    self.message.emit(buf[:])
##                    buf = []
##                    t0 = time()
#            except:
#                pass
#
#    def finish(self):
#        self.socket.close()
#        self.context.term()


class DataReceiver(QtCore.QObject):
    msgClear = QtCore.pyqtSignal()
    msgMoveable = QtCore.pyqtSignal(list)
    msgCounters = QtCore.pyqtSignal(list)
    msgRepeats = QtCore.pyqtSignal(list)
    msgHintE0 = QtCore.pyqtSignal(list)

    msgRawData = QtCore.pyqtSignal(list)
    msgMuData = QtCore.pyqtSignal(list)

    def __init__(self):
#        port = 5567
        QtCore.QObject.__init__(self)
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.SUB)
        self.socket.connect(f"tcp://{ZMQ_HOST}:{ZMQ_PORT}")
        self.socket.setsockopt_string(zmq.SUBSCRIBE, "")
        self.running = True
        self.isScanInitialized = False
#        if isTest:
#            self.testTimer = QtCore.QTimer()
#            self.testTimer.setSingleShot(True)
##            self.testTimer.timeout.connect(self.test_random_listener)
##            self.testTimer.timeout.connect(self.test_data_listener)
#            self.testTimer.timeout.connect(self.test_data_listener_FY)
#            self.testTimer.start(1000)
#            self.door = None
#        else:
#            import taurus
#            self.door = taurus.Device('Balder/Door/01')
#            dp = taurus.Attribute('Balder/Door/01/RecordData')
#            dp.addListener(self.listener)

#    def listener(self, evt_src, evt_type, evt_value):
    def listener(self, evt_src, evt_type, evt_value):
        """prefix, name, doc: beamline_name, event type (start, descriptor,
            event, stop), event data"""
#        if evt_type != 0:
#            return
#        print(evt_src, evt_type, str(evt_type))
        dataJSON = pickle.loads(evt_value)
        try:
            if evt_type == b'start':  # Preparing the plot
                self.msgClear.emit()
                print("start event", dataJSON)
                planName = str(dataJSON['plan_name'])
                self.isMultiPart = False

                if "continuous" in planName:
                    self.isMultiPart = True
                    self.tmpData = {}

                if isinstance(dataJSON['motor'], list):
                    moveable = dataJSON['motor'][0]
                else:
                    moveable = dataJSON['motor']

                movMin = dataJSON['plan_args']['start']
                movMax = dataJSON['plan_args']['stop']

                if 'energy' in moveable.lower():
                    unit = 'eV'
                elif 'pitch' in moveable.lower():
                    unit = 'deg'
#                elif moveable.endswith(('_rol', '_yaw')):  # mirrors
#                    unit = 'mrad'
#                elif moveable.endswith(('rol', 'yaw')):  # crystals
#                    unit = u'µrad'
#                elif moveable.endswith(('_x', '_y', '_z')):
#                    unit = 'mm'
                else:
                    unit = ''
                self.msgMoveable.emit([moveable, unit, movMin, movMax])
                self.moveable = moveable
                self.minmax = [movMin, movMax]

                print(moveable, unit, movMin, movMax)

                cntField = 'monitors' if self.isMultiPart else 'dark_currents'
                allCounters = list(dataJSON[cntField])

#                allCounters = list(dataJSON['counters'])
                counters, labels, units = [], [], []
                kinds, locations, colors = [], [], []
                iTR, iFY = 0, 0

#                for cd, c in zip(column_desc[3:], allCounters):
                for counterString in allCounters:
                    if not counterString.startswith("Scaler"):
                        continue
#                    cond = cd['conditioning'].lower()
                    kind = SIGNAL_KIND_UNKNOWN
                    color = COLOR_DEFAULT
                    for label in LABELS_TR:
                        if label in counterString:  # or label in cond:
                            kind = SIGNAL_KIND_TR
                            iTR += 1
                            color = colorCycleT[(iTR-1) % len(colorCycleT)]
                            break
                    for label in LABELS_FY:
#                        print(label, cond)
                        if label in counterString:  # or label in cond:
                            kind = SIGNAL_KIND_FY
                            iFY += 1
                            color = colorCycleF[(iFY-1) % len(colorCycleF)]
#                            print("found fy")
                            break
#                    if 'skip' in cond or 'hid' in cond:
#                        continue  # can't be right after 'for' for right color
#                    isRight = 'right' in cond
#                    if isRight and kind == SIGNAL_KIND_TR:
#                        loc = 2
#                    elif not isRight and kind == SIGNAL_KIND_FY:
#                        loc = 3
#                    elif isRight and kind == SIGNAL_KIND_FY:
#                        loc = 4
#                    else:
#                        loc = 1
                    if "16" in counterString:
                        loc = 1
                    elif "17" in counterString:
                        loc = 1
                    elif "18" in counterString:
                        loc = 1
                    else:
                        loc = 3

                    unit = 'counts'
                    label = SCALER_ALIAS[counterString]

                    counters.append(counterString)
                    labels.append(label)
                    units.append(unit)
                    kinds.append(kind)
                    locations.append(loc)
                    colors.append(color)

                self.counters = counters
                self.locations = locations
                self.msgCounters.emit([labels, units, locations, colors])
                print(counters, labels, units, locations, colors)

                serialNo = dataJSON['scan_id']
                repeat, repeats = None, None
#                try:
#                    if self.door is not None:
#                        repeat = self.door.get_property('repeat')['repeat']  # a list of str  # noqa
#                        if repeat != []:
#                            repeat = repeat[0]
#                        repeats = self.door.get_property('repeats')['repeats']
#                        if repeats != []:
#                            repeats = repeats[0]
#                except IndexError:
#                    pass
                self.msgRepeats.emit([serialNo, repeat, repeats])
                try:  # TODO: add E0 hint to the plan
                    hintE0 = dataJSON['e0']
                    if hintE0 != []:
                        self.msgHintE0.emit(
                            [float(hintE0)-20, float(hintE0)+20])
                except:
                    pass

                self.isScanInitialized = True
            elif evt_type == b'event':
                if not self.isScanInitialized:
                    return

                data = None
                checkFPP = False
                if self.isMultiPart:
#                    print("Multipart event")
                    allFields = [motor_monitors[self.moveable]] + self.counters
                    evtNum = str(dataJSON['seq_num'])
#                    print(evtNum, allFields)

                    if not evtNum in self.tmpData.keys():
                        self.tmpData[evtNum] = {}

                    evtData = dataJSON['data']
#                    print(evtData.items())
                    for key, value in evtData.items():
#                        print(key, self.counters)
                        if key in allFields:
                            self.tmpData[evtNum][key] = value
#                            print(self.tmpData[evtNum])
                            if len(self.tmpData[evtNum]) == len(allFields):
                                data = [self.tmpData[evtNum][d] for d in allFields]
#                                print("Data", data)
                                checkFPP = int(evtNum) == 1 and abs(data[0]-self.minmax[0]) > 1
                                if not checkFPP:
                                    self.msgRawData.emit(data)

                                try:
                                    del self.tmpData[evtNum]
                                except:
                                    pass
                                break

                else:
    #                print("event:", dataJSON)
                    allFields = [motor_dict[self.moveable]] + self.counters
                    data = [evtData[d] for d in allFields]
    #                print("data", data)
    #                for i, c in enumerate(self.counters):
    #                    if 'aem' in c:
    #                        data[i+1] *= 1e-3  # ind 0 is moveable
                    self.msgRawData.emit(data)

                if data is not None:
                    if 'energy' in self.moveable.lower():
                        i0 = data[1] # + data[2]
                        try:
                            muTData = np.log(i0 / data[2])
                        except:  # noqa
                            muTData = 0.
                        muData = [data[0], muTData]
    #                    print("mu", muData)
                        for d, loc in zip(data[1:], self.locations):
    #                        if loc in [4]:  # Fluorescence
                            if loc == 3:  # Fluorescence
                                try:
                                    muFData = d / i0
                                except:  # noqa
                                    muFData = 0.
    #                            print("muFdata", muFData)
                                muData.append(muFData)
                        if not checkFPP:
                            self.msgMuData.emit(muData)  # 1st is T
        except:  # noqa
            raise
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
        fileNames = ['test_data/S90T3_spot6_010.dat']
        for iFile, fileName in enumerate(fileNames):
            self.msgClear.emit()
            e, ch1, ch2, ch3, ch4 = np.loadtxt(
                fileName, unpack=True, usecols=(1, 3, 4, 5, 6))
            moveable = 'energy'
            movMin = e.min()
            movMax = e.max()
            self.msgMoveable.emit([moveable, 'eV', movMin, movMax])
            counters = ['ch1', 'ch2', 'ch3', 'ch4']
            locations = [1, 1, 1, 3]
            units = [u'A', u'A', u'A', u'A']
            colors = colorCycleT[0:3] + colorCycleF[0:1]
            self.msgCounters.emit([counters, units, locations, colors])
            self.msgRepeats.emit([None, iFile+1, len(fileNames)])
            self.msgHintE0.emit([7112-50, 7112+50])

            print("started")
            t0 = time.time()
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

    def loop(self):
#        buf = []
#        t0 = time()
        print("Running")
        while self.running:
            try:
                message = self.socket.recv()
                prefix, name, doc = message.split(b' ', 2)
                self.listener(prefix, name, doc)
#                self.
                #document = json.loads(message)
#                print("Received Document:", prefix, name, pickle.loads(doc))

#                string = self.socket.recv()
#                buf.append(string)
#                if time() - t0 > 1./FPS_LIMIT:
#                    self.message.emit(buf[:])
#                    buf = []
#                    t0 = time()
            except:
                raise
                pass

    def finish(self):
        self.socket.close()
        self.context.term()


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
        return c*f**3 + d

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
        try:
            mu = self.muT
            kres = self.k(self.e, mu)
            if kres is None:
                self.terminate()
                return
            ie0, k, kmin_idx = kres
            pre_edge = self.pre_edge(self.e, ie0, mu)
            if pre_edge is None:
                self.terminate()
                return
            chi = self.chi(self.e, ie0, mu, pre_edge, k, kmin_idx)
            winFT = self.hann_win(k, 2., max(k), 1.0)
            ft = self.make_ft(k, chi, winFT)
            ks = [k]
            chis = [chi]
            fts = [ft]

            for mu in self.muFs:
                kres = self.k(self.e, mu)
                if kres is None:
                    self.terminate()
                    return
                ie0, k, kmin_idx = kres
                pre_edge = self.pre_edge(self.e, ie0, mu)
                if pre_edge is None:
                    self.terminate()
                    return
                chi = self.chi(self.e, ie0, mu, pre_edge, k, kmin_idx)
                winFT = self.hann_win(k, 2., max(k), 1.0)
                ft = self.make_ft(k, chi, winFT)
                ks.append(k)
                chis.append(chi)
                fts.append(ft)
            res = [ks, chis, fts, winFT]
        except (TypeError, IndexError, ValueError):
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

        self.threadRecv = QtCore.QThread()
#        self.zeromq_listener = ZeroMQ_Listener()
#        self.zeromq_listener.moveToThread(self.thread)
#        self.thread.started.connect(self.zeromq_listener.loop)
#        self.zeromq_listener.message.connect(self.signal_received)
#        QtCore.QTimer.singleShot(0, self.thread.start)

        self.dataReceiver = DataReceiver()
        self.dataReceiver.msgClear.connect(self.setClear)
        self.dataReceiver.msgMoveable.connect(self.setMoveable)
        self.dataReceiver.msgCounters.connect(self.setCounters)
        self.dataReceiver.msgRepeats.connect(self.setRepeats)
        self.dataReceiver.msgHintE0.connect(self.setHintE0)
        self.dataReceiver.msgRawData.connect(self.setRawData)
        self.dataReceiver.msgMuData.connect(self.setMuData)

        self.dataReceiver.moveToThread(self.threadRecv)
        self.threadRecv.started.connect(self.dataReceiver.loop)
#        self.zeromq_listener.message.connect(self.signal_received)
        QtCore.QTimer.singleShot(0, self.threadRecv.start)

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
        print("hint E0", lst)
        self.exafs_extractor.hintE0 = lst

    def setRawData(self, lst):
        """lst = [moveable] + counters"""
        self.moveableValue.append(lst[0])
        if not self.counterValues:
            self.counterValues = [[] for i in lst[1:]]
        if len(lst[1:]) != len(self.counterValues):
            raise ValueError("wrong length of data")
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
        self.emu.append(lst[0])
        self.muT.append(lst[1])
        if self.isFYenabled:
            if not self.muFs:
                self.muFs = [[] for i in lst[2:]]
            if len(lst[2:]) != len(self.muFs):
                raise ValueError("wrong length of mu data")
            for data, muf in zip(lst[2:], self.muFs):
                muf.append(data)

        if self.tabWidget.currentIndex() == 1:
            if self.nPointsReceived % FRAME_CHUNK_MU == 0:
                self.updatePlots()
        elif self.tabWidget.currentIndex() == 2:
            if self.nPointsReceived % FRAME_CHUNK_EXAFS == 0:
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
                print(len(muFs), len(self.muFs), len(self.curvesMuF))
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


if __name__ == '__main__':
    import sys
    app = QtGui.QApplication([])
    mw = MonitorWidget()
#    mw.setStyleSheet("background-color: #000")
    mw.setWindowTitle("Data monitor")
    mw.show()
    sys.exit(app.exec_())
