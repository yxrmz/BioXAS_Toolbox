# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 13:31:28 2019

@author: chernir
"""
import os
#os.environ["CDF_LIB"] = '/local/dorofel/Software/cdf/lib'

from spacepy import pycdf
#import h5py
import numpy as np
#import sys
#import os
import sqlite3
import re
from collections import OrderedDict

#def applyROI(roiMinMax, )

#  Directory must contain the userdata.db file used by Acquaman
#dirname = r"X:\bioxas-m\AcquamanMainData\users\36G12634\userData"
dirname = r"X:\bioxas-m\AcquamanMainData\users\36G12815\userData"
#  Output directory. Will be created if not exists
outputDir = r"D:\BioXAS\36G12815"

#  One may convert only the scan with certain names, i.e. "MPE_Nickel%"
#  Do not specify the ordinal number of the scan here
#  Use the % wildcard to process all scans
#scanNameFilter = "%"
scanNameFilter = "W%"
#scanNameFilter = "Test datasource 4%"

#  Specific scan type can be converted
#  possible values for BioXAS are
#  'BioXASXAS%' for EXAFS/XANES scans
#  'AMXRF%' for XRF scans
#  'BioXASGeneric%' for generic scans
scanTypeFilter = "BioXASXAS%"
#scanTypeFilter = "BioXASGeneric%"
#scanTypeFilter = "AMXRF%"
#scanTypeFilter = "%"
#dirname = "/local/dorofel/localBeamlineData/BioXAS/bioxas-m/AcquamanMainData/users/30-10721/userData"
#outputDir = "/local/dorofel/tmp/HDF5/BioXAS"


dbpath = os.path.join(dirname, 'userdata.db')
if not os.path.isfile(dbpath):
    print(f"No such file {dbpath}")
    exit(1)
conn = sqlite3.connect(dbpath)

conn.row_factory = sqlite3.Row
cursor = conn.cursor()

qu = f"""SELECT name, number, filePath, id, scanInitialConditions, scanConfiguration, dateTime, endDateTime
    FROM AMScan_table WHERE name LIKE '{scanNameFilter}' AND scanConfiguration LIKE '{scanTypeFilter}'"""
cursor.execute(qu)
rows = cursor.fetchall()

# Loop filenames
#totalLen = len(rows)
for irow, row in enumerate(rows):
    print(f"Processing record {irow+1} out of {len(rows)}")
    # Extract session name and date
    AMsession = {}
    qu = "SELECT name as session_name, MAX(dateTime) as session_start FROM AMRun_table WHERE dateTime<=?"
    cursor.execute(qu, (row['dateTime'],))
    sess = cursor.fetchone()
    if sess:
        AMsession.update({k: sess[k] for k in sess.keys()})

    DCdict = {}
    ROIdict = {}
    BCdict = {}
    SRdict = {}
    ELdict = {}
    for k in row.keys():
        print(f"{k}: {row[k]}")
    cdfPath = os.path.join(dirname, row['filePath'])

    scanid = (row['id'],)

    # Extract dark currents
#    qu_dc = """SELECT dataName, darkCurrent, timeUnitMultiplier FROM AM1DDarkCurrentCorrectionAB_table
#     JOIN AMScan_table_analyzedDataSources ON AMScan_table_analyzedDataSources.id2=AM1DDarkCurrentCorrectionAB_table.id
#     WHERE table2='AM1DDarkCurrentCorrectionAB_table' AND id1 LIKE ?"""

    qu_dc = """SELECT dataName, darkCurrent, timeUnitMultiplier FROM AM1DDarkCurrentCorrectionAB_table,
     AMScan_table_analyzedDataSources WHERE table2='AM1DDarkCurrentCorrectionAB_table' AND
     AMScan_table_analyzedDataSources.id2=AM1DDarkCurrentCorrectionAB_table.id AND id1 LIKE ?"""

    for dcrow in cursor.execute(qu_dc, scanid):
#        print(dcrow['dataName'])
        DCdict[dcrow['dataName']] = (dcrow['darkCurrent'], dcrow['timeUnitMultiplier'])
#    sys.exit()
    # Extract ROI info
#    qu_roi = """SELECT AMRegionOfInterestAB_table.name as roi_name,
#     binningRangeLowerBound, binningRangeUpperBound FROM AMRegionOfInterestAB_table
#     JOIN AMScan_table_analyzedDataSources ON AMScan_table_analyzedDataSources.id2=AMRegionOfInterestAB_table.id
#     WHERE table2='AMRegionOfInterestAB_table' AND id1=?"""

#    qu_roi = """SELECT DISTINCT AMRegionOfInterestAB_table.name as roi_name,
#     binningRangeLowerBound, binningRangeUpperBound, AMRegionOfInterest_table.energy FROM AMRegionOfInterestAB_table,
#     AMRegionOfInterest_table, AMScan_table_analyzedDataSources WHERE AMRegionOfInterestAB_table.id=AMScan_table_analyzedDataSources.id2 AND
#     AMRegionOfInterest_table.id=AMRegionOfInterestAB_table.id1 AND AMScan_table_analyzedDataSources.id1=?"""

    qu_roi = """SELECT DISTINCT AMRegionOfInterestAB_table.name as roi_name,
     binningRangeLowerBound, binningRangeUpperBound FROM AMRegionOfInterestAB_table,
     AMScan_table_analyzedDataSources WHERE AMRegionOfInterestAB_table.id=AMScan_table_analyzedDataSources.id2 AND
     AMScan_table_analyzedDataSources.id1=?"""

    for roirow in cursor.execute(qu_roi, scanid):
        ROIdict[roirow['roi_name']] = (roirow['binningRangeLowerBound'], roirow['binningRangeUpperBound'])

    exclude_keys = ['id', 'name', 'AMDbObjectType', 'thumbnailCount', 'thumbnailFirstId']
    exclude_values = [None, '']
    # Extract beamline configuration
    table1, idx1 = row['scanInitialConditions'].split(';')
    query_params = (idx1, table1)
    qu_bc = """SELECT AMControlInfo_table.* FROM AMControlInfo_table JOIN AMControlInfoList_table_controlInfos
     ON AMControlInfo_table.id=AMControlInfoList_table_controlInfos.id2
     WHERE AMControlInfoList_table_controlInfos.id1=? AND AMControlInfoList_table_controlInfos.table1=?"""

    for row1 in cursor.execute(qu_bc, query_params):
        ctrl_name = row1['name']
        BCdict[ctrl_name] = {}
        for k1 in row1.keys():
            if k1 not in exclude_keys:
                if row1[k1] not in exclude_values:
                    BCdict[ctrl_name].update({k1: row1[k1]})

    # Extract scan regions
    table2, idx2 = row['scanConfiguration'].split(';')
    query_params = (idx2, table2)
#    print(query_params)
    
    qu_sr = "SELECT AMScanAxis_table_axisRegions.id2, AMScanAxis_table_axisRegions.table2"
    qu_sr += f" FROM AMScanAxis_table_axisRegions JOIN {table2}_scanAxes ON"
    qu_sr += f" {table2}_scanAxes.id2=AMScanAxis_table_axisRegions.id1"
    qu_sr += f" WHERE {table2}_scanAxes.id1=? AND {table2}_scanAxes.table1=?"
    try:
        cursor.execute(qu_sr, query_params)
        rows1 = cursor.fetchall()
        for row1 in rows1:
            params2 = (row1['id2'],)
            qu = f"SELECT * FROM {row1['table2']} WHERE id=?"
            for row2 in cursor.execute(qu, params2):
                region_name = row2['name']
                SRdict[region_name] = {}
                for k2 in row2.keys():
                    if k2 not in exclude_keys:
                        if row2[k2] not in exclude_values:
                            SRdict[region_name].update({k2: row2[k2]})
    except sqlite3.OperationalError:
        # there might be no such table for certain scan types
        pass


    qu_sr = f"SELECT configurationDbObject, header FROM {table2} WHERE id={idx2}"
    xasheader, scan_energy, scan_edge = None, None, None
    try:
        cursor.execute(qu_sr)
        rows1 = cursor.fetchall()
        for row1 in rows1:
            xasheader = row1['header']
            dboTable, dbot_row = row1['configurationDbObject'].split(';')
            for row2 in cursor.execute(
                    f"SELECT energy, edge FROM {dboTable} WHERE id={dbot_row}"):
                scan_energy, scan_edge = row2['energy'], row2['edge']
#        print(scan_energy, scan_edge)
    except sqlite3.OperationalError:
        # there might be no such table for certain scan types
        pass

    query_params = (idx2,)
    qu_sr = f"""SELECT name, energy, lowerBound, upperBound FROM AMRegionOfInterest_table, {table2}_regionsOfInterest
    WHERE AMRegionOfInterest_table.id={table2}_regionsOfInterest.id2 AND
    {table2}_regionsOfInterest.id1=?"""

    try:
        cursor.execute(qu_sr, query_params)
        rows1 = cursor.fetchall()
        for row1 in rows1:
            lineName = row1['name']
            ELdict[lineName] = {'energy': row1['energy'],
                                'lowerBound': row1['lowerBound'],
                                'upperBound': row1['upperBound']}
    except sqlite3.OperationalError:
        # there might be no such table for certain scan types
        pass

    h5FileName = f"{row['name']}_{row['number']}.dat"

    try:
        os.makedirs(outputDir)
    except FileExistsError:
        pass

    header = """XDI/1.0\n"""
#    header = """XDI/1.0\nEndstation\nBioXAS Beamline - Main Endstation\n\n"""
    if scan_edge:
        header += "Element symbol: {0}\nElement edge: {1}\n".format(*(scan_edge.split(" ")))
#    with h5py.File(os.path.join(outputDir, h5FileName), 'w') as h5file:
    edgeEnergy = None
    print(f"Converting {cdfPath} to {h5FileName}")
    if SRdict:
#            scanregions = h5file.create_group(u'scan_regions')
        for sr_key, sr_value in SRdict.items():
#                srentry = scanregions.create_group(f'{sr_key}')
            for rk, rv in sr_value.items():
#                    if isinstance(rv, str):
#                        rv = rv.encode('utf-8')
                if str(rk) == "edgeEnergy":
                    edgeEnergy = float(rv)
                    break
            break

        if ELdict:
#            emissionLines = h5file.create_group(u'emission_lines')
            for el_key, el_value in ELdict.items():
                print(el_key)
#                eLines
#                elgrp = emissionLines.create_group(el_key)
#                for elval_key, elval_val in el_value.items():
#                    elgrp.create_dataset(elval_key, data=np.array(elval_val))


    dataOut = OrderedDict()

    if True:
#            roientry.create_dataset(roi_key, data=np.array(roi_value))
        
        with pycdf.CDF(cdfPath) as cdf:
            sliceStartGlo = 0
            for key in cdf.keys():
                if len(cdf['ScalerDwellTimeFeedback']) > len(cdf[key]):
                    sliceStartGlo = 1
                    break
            for col in ['EnergyFeedback',
                        'EnergySetpoint',
                        'ScalerDwellTimeFeedback']:
                dataOut[col] = np.array(cdf[col])[sliceStartGlo:]
            for ikey, key in enumerate(cdf.keys()):
                if len(re.findall('ICR|Raw', str(key))) > 0 or key in dataOut.keys():
                    continue
        #        if ikey == 0:
        #            print("Total counts in channel 1:", np.sum(np.array(cdf[key])))
                if len(cdf[key]) > 0:
#                    print(key)
                    if str(key).endswith("boardDetector"):
                        data2D = np.array(cdf[key])
#                        print(data2D.shape)
                        for roi_key, roi_value in ROIdict.items():
                            if str(roi_key).endswith("InB") or str(roi_key).endswith("OutB"):
#                                if float(roi_value[-1]) < edgeEnergy:
                                if str(roi_key).startswith(scan_edge.split(" ")[0]):
                                    colname = roi_key
                                    data1D = np.sum(
                                            data2D[:, int(roi_value[0]/10):int(roi_value[-1]/10)], axis=1)
                                    dataOut[colname] = data1D
#                                    print(colname, data1D.shape, edgeEnergy, roi_value)
                


#                    nxentry.create_dataset(key, data=np.array(cdf[key]))
                    if key in DCdict.keys():
                        print("Found dark currents for ", key)
                        sliceStart = 1 if len(cdf['ScalerDwellTimeFeedback']) > len(cdf[key]) else 0
                        c_name = "{}_darkCorrected".format(key)
                        data=np.array(cdf[key]) -\
                             (DCdict[key][0]*DCdict[key][1] *
                             np.array(cdf['ScalerDwellTimeFeedback']))[sliceStart:]
                        dataOut[c_name] = data
#                        print(c_name)
#                        nxentry.create_dataset(
#                             "{}_darkCorrected".format(key),
#                             data=np.array(cdf[key]) -
#                             (DCdict[key][0]*DCdict[key][1] *
#                             np.array(cdf['ScalerDwellTimeFeedback']))[sliceStart:])
    #                    nxentry.create_dataset(
    #                        "{}_darkCorrection".format(key),
    #                        data=DCdict[key][0] * DCdict[key][1] *
    #                             np.array(cdf['ScalerDwellTimeFeedback']))

#        nxentry.attrs.create("start_time", np.string_(f"{row['dateTime']}"))
#        nxentry.attrs.create("end_time", np.string_(f"{row['endDateTime']}"))
#        if AMsession:
#            for sk in AMsession:
#                nxentry.attrs.create(sk.encode('utf-8'), np.string_(f"{AMsession[sk]}"))
#
#        dcentry = h5file.create_group(u'dark_currents')
#        for dc_key, dc_value in DCdict.items():
#            dcentry.create_dataset(dc_key, data=np.array(dc_value))
#        roientry = h5file.create_group(u'ROIs')
#        for roi_key, roi_value in ROIdict.items():
#            roientry.create_dataset(roi_key, data=np.array(roi_value))
#
#        bcall = h5file.create_group(u'beamline_configuration')
#        for bc_key, bc_value in BCdict.items():
#            bcentry = bcall.create_group(f'{bc_key}')
#            for ck, cv in bc_value.items():
#                if isinstance(cv, str):
#                    cv = cv.encode('utf-8')
#                bcentry.create_dataset(ck.encode('utf-8'), data=np.array(cv,))
#        if SRdict:
#            scanregions = h5file.create_group(u'scan_regions')
#            for sr_key, sr_value in SRdict.items():
#                srentry = scanregions.create_group(f'{sr_key}')
#                for rk, rv in sr_value.items():
#                    if isinstance(rv, str):
#                        rv = rv.encode('utf-8')
#                    srentry.create_dataset(rk.encode('utf-8'), data=np.array(rv, ))
#        if ELdict:
#            emissionLines = h5file.create_group(u'emission_lines')
#            for el_key, el_value in ELdict.items():
#                elgrp = emissionLines.create_group(el_key)
#                for elval_key, elval_val in el_value.items():
#                    elgrp.create_dataset(elval_key, data=np.array(elval_val))

#    print(cdf.keys())
#    print(DCdict)
#    print(ROIdict)
    for icol, col in enumerate(dataOut.keys()):
        header += f"Column {icol}: {col}\n"
    
    header += "/////////////////////////////////\n"
    header += f"Scan: {row['name']} #{row['number']}\n"
    header += f"Date: {row['dateTime']}\n\n"
    
    dataArr = None
    datacol = "\n----------------------------------\n"
    for colname, col in dataOut.items():
        datacol += "{}\t".format(colname)
        print(colname)
        dataArr = np.vstack((dataArr, col)) if dataArr is not None else np.array(col)
    
    if xasheader:    
        header += xasheader
    
    header += datacol
    
    np.savetxt(os.path.join(outputDir, h5FileName),
               dataArr.T,
               header=header,
               delimiter="\t")

