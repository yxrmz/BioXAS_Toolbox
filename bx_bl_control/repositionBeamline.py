# -*- coding: utf-8 -*-
"""
Created on Sun Jul 23 17:43:44 2017

@author: bioxasmain
"""

from blElements import BeamLine, Motor, MotorK


#M1_DS = MotorK(pvStr='SMTR1607-5-I21-03', fieldStr="mm", fbk="fbk",
#               kStr="Bl1607-I21:M1:PCS:sensors:C", kOffset=8.108145)
#print(M1_DS.position, M1_DS.positionK)


mainBL = BeamLine(configuration='Main')
mainBL.repositionToEnergy(15000, silent=True)


