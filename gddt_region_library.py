''' gddt_region_library
Region library for growing degree days and treeline analysis.
Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
import sys

from icecream import ic
import numpy as np

########################################################################
####  ECOREGIONS 2017 DEFINITIONS
########################################################################
def ArcticCoastalTundra():
    reg_dict = {
        "reg_str": 'Arctic Coastal Tundra',
        "reg_abv": 'ArctCoastTundra',
        "reg_lats": np.array(
            [69.43, 70.32, 70.89, 71.36, 70.85, 70.19, 69.89]),
        "reg_lons": np.array(
            [-163.081, -161.851, -159.346, -156.577, -152.314, -145.986, 
            -142.207]) % 360,
    }
    
    return reg_dict

def ArcticFoothillsTundra():
    reg_dict = {
        "reg_str": 'Arctic Foothills Tundra',
        "reg_abv": 'ArctFootTun',
        "reg_lats": np.array([68., 68., 70., 70., 68.87, 68.]),
        "reg_lons": np.array(
            [-163., -167., -167., -147.74, -147.74, -163.]) % 360,
    }
    
    return reg_dict

def BrooksBritishRange():
    ''' Brooks-British Range Tundra '''
    reg_dict = {
        "reg_str": 'Brooks Range',
        "reg_abv": 'BrooksRange',
        "reg_lats": np.array(
            [67.36, 68.69, 69.57, 68.55, 65.95, 65.88, 68.20]),
        "reg_lons": np.array(
            [-162.729, -162.729, -140., -135.835, -135.835, -136.67, 
             -139.28]) % 360,
    }
    
    return reg_dict

def CanadianLowArcticTundra():
    ''' Canadian Low Arctic Tundra '''
    reg_dict = {
        "reg_str": 'Canadian Low Arctic Tundra',
        "reg_abv": 'CanLowArctTun',
        "reg_lats": np.array([69., 59.72, 63.36, 65.34, 68.34]),
        "reg_lons": np.array(
            [-135.75, -94.88, -90.308, -96.548, -93.47]) % 360,
        }
    
    return reg_dict

def atlas_ecoregions2017():
    ''' Return all Ecoregions 2017 (Dinerstein et al. 2017) '''
    atlas = (ArcticCoastalTundra(), ArcticFoothillsTundra(), 
             BrooksBritishRange(), CanadianLowArcticTundra())
    return atlas

########################################################################
####  CUSTOM REGION DEFINITIONS
########################################################################
def AntarcticCoast():
    reg_dict = {
        "reg_str": 'Antarctic Coast',
        "reg_abv": 'AntarcticCoast',
        "reg_lats": np.array([-72, -60]),
        "reg_lons": np.array([0, 359])
    }

    return reg_dict

def Arctic50N():
    reg_dict = {
        "reg_str": 'Arctic 50+N',
        "reg_abv": 'Arctic50N',
        "reg_lats": np.array([50, 90]),
        "reg_lons": np.array([0, 359])
    }
    
    return reg_dict

def AyiyakRiver():
    reg_dict = {
        "reg_str": 'Ayiyak River',
        "reg_abv": 'Ayiyak',
        "reg_lats": np.array([68.53]),
        "reg_lons": np.array([208.31])
    }
    
    return reg_dict    
    
def BelangerIslandCell():
    reg_dict = {
        "reg_str": 'Belanger Island Cell',
        "reg_abv": 'BelangerIslCell',
        "reg_lats": np.array([56.07]),
        "reg_lons": np.array([283.54])
    }

    return reg_dict

def BelangerIslandNorth():
    reg_dict = {
        "reg_str": 'Belanger Island North',
        "reg_abv": 'BelangerIslN',
        "reg_lats": np.array([55, 91]),
        "reg_lons": np.array([192, 348])
    }

    return reg_dict

def BrooksRange():
    reg_dict = {
        "reg_str": 'Brooks Range',
        "reg_abv": 'BrooksRange',
        "reg_lats": np.array([65.5, 69.5]),
        "reg_lons": np.array([195, 209.5])
    }

    return reg_dict

def BrooksRange_point():
    reg_dict = {
        "reg_str": 'Brooks Range point',
        "reg_abv": 'BrooksRangePoint',
        "reg_lats": np.array([68.5]),
        "reg_lons": np.array([196.])
    }

    return reg_dict

def BrooksRange_colonist():
    reg_dict = {
        "reg_str": 'Brooks Range (67.38N, 202.5E)',
        "reg_abv": 'BrooksRangeCol',
        "reg_lats": np.array([67.5]),
        "reg_lons": np.array([202.])
    }

    return reg_dict
    
def BrooksCoast():
    reg_dict = {
        "reg_str": 'Brooks Coast',
        "reg_abv": 'BrooksToCoast',
        "reg_lats": np.array([70, 73]),
        "reg_lons": np.array([195, 209.5])
    }

    return reg_dict
    
def CentralSiberia():
    #  Estimated from Dinerstein et al. 2017
    reg_dict = {
        "reg_str": 'Central Siberia',
        "reg_abv": 'CentralSiberia',
        "reg_lats": np.array([61.47, 73.30, 77.88, 75.54, 70.67, 59.17, 55.05]),
        "reg_lons": np.array([67.06, 74.37, 106.9, 149.71, 158.48, 156.58, 134.89])
    }
    
    return reg_dict
    
def ChukchiPeninsula():
    #  Estimated from Dinerstein et al. 2017
    reg_dict = {
        "reg_str": 'Chukchi Peninsula',
        "reg_abv": 'Chukchi',
        "reg_lats": np.array([71.81, 66.01, 59.94, 50.86, 52.75, 60.74, 70.78]),
        "reg_lons": np.array([181.81, 190.71, 170.35, 158.07, 153.14, 160.98, 158.13])
    }
    
    return reg_dict
    
def EdmontonIsh():
    reg_dict = {
        "reg_str": 'Edmunton (roughly)',
        "reg_abv": 'Edmunton',
        "reg_lats": np.array([53.55]),
        "reg_lons": np.array([246.45])
    }
    
    return reg_dict
    
def EurasianTundra():
    #  Estimated from Dinerstein et al. 2017
    reg_dict = {
        "reg_str": 'Eurasian Tundra',
        "reg_abv": 'EurasianTundra',
        "reg_lats": np.array(
            [58.27, 62.23, 79.59, 82.58, 80.53, 66.09, 50.00, 55.15, 64.55,
             69.47, 69.21, 64.64, 69.51]),
        "reg_lons": np.array(
            [6.90, 5.03, 10.57, 64.23, 145.29, 190.91, 160.00, 154.69, 171.25,
            133.56, 100.08, 71.69, 25.08]),
    }
    
    return reg_dict
    
def GreenlandIsh():
    reg_dict = {
        "reg_str": 'GreenlandIsh',
        "reg_abv": 'GreenlandIsh',
        "reg_lats": np.array([60, 90]),
        "reg_lons": np.array([288, 350.5])
    }

    return reg_dict
    
def HighArctic():
    reg_dict = {
        "reg_str": 'HighArctic',
        "reg_abv": 'HighArctic',
        "reg_lats": np.array([66, 90]),
        "reg_lons": np.array([0, 359])
    }

    return reg_dict
    
def Iceland():
    reg_dict = {
        "reg_str": 'Iceland',
        "reg_abv": 'Iceland',
        "reg_lats": [63.3, 66.5],
        "reg_lons": [335.39, 346.5]
    }
    
    return reg_dict
    
def Ilulissat():
    reg_dict = {
        "reg_str": 'Ilulissat',
        "reg_abv": 'Ilulissat',
        "reg_lats": np.array([69.22]),
        "reg_lons": np.array([308.95])
    }

    return reg_dict
    
def IPCC_NorthEurope():
    reg_dict = {
        "reg_str": 'IPCC North Europe',
        "reg_abv": 'IPCCNorthEurope',
        "reg_lats": (np.array([48, 75, 75, 51]), np.array([51, 75, 75, 61.3])),
        "reg_lons": (np.array([350, 350, 360, 360]), np.array([-1, -1, 40, 40]))
    }

    return reg_dict
    
def KotelnyIsland():
    reg_dict = {
        "reg_str": 'Kotelny Island',
        "reg_abv": 'KotelnyIsland',
        "reg_lats": np.array([75.41]),
        "reg_lons": np.array([142.0])
    }

    return reg_dict

def NAmArctic():
    reg_dict = {
        "reg_str": 'North American Arctic',
        "reg_abv": 'NAmArctic',
        "reg_lats": np.array([59, 91]),
        "reg_lons": np.array([192, 348])
    }

    return reg_dict
    
def ScandinaviaAndWesternRussia():
    #  Estimated from Dinerstein et al. 2017
    reg_dict = {
        "reg_str": 'Scandinavia and Western Russia',
        "reg_abv": 'ScAndWestRus',
        "reg_lats": np.array([58.34, 62.14, 70.67, 73.44, 61.73, 61.88]),
        "reg_lons": np.array([5.86, 4.57, 18.99, 75.33, 66.87, 48.72])
    }
    
    return reg_dict
    
def SiberiaForestBox():
    #  Definition from Tchebakova et al. 2009
    reg_dict = {
        "reg_str": 'Siberian Forests',
        "reg_abv": 'SiberiaForestBox',
        "reg_lats": np.array([50, 75]),
        "reg_lons": np.array([60, 140])
    }
    
    return reg_dict

def SiberiaGuesstimate():
    #  All "Siberia" refers to this one before August ~5 2024
    reg_dict = {
        "reg_str": 'SiberiaGuesstimate',
        "reg_abv": 'SiberiaGuesstimate',
        "reg_lats": np.array([70, 80]),
        "reg_lons": np.array([80, 150])
    }
    
    return reg_dict
        
def SiberiaNorthBorealForest():
    reg_dict = {
        "reg_str": 'Siberia North of Boreal Forest',
        "reg_abv": 'SibNorBorFor',
        "reg_lats": np.array([60, 67.5]),
        "reg_lons": np.array([51.4, 82.4])
    }
    
    return reg_dict

def SiberiaPoint():
    reg_dict = {
        "reg_str": 'SiberiaPoint',
        "reg_abv": 'SiberiaPoint',
        "reg_lats": np.array([66.5]),
        "reg_lons": np.array([70.4])
    }
    
    return reg_dict
    
def Svalbard():
    reg_dict = {
        "reg_str": 'Svalbard',
        "reg_abv": 'Svalbard',
        "reg_lats": np.array([78.14]),
        "reg_lons": np.array([22.12])
    }
    
    return reg_dict

def Tasiilaq():
    reg_dict = {
        "reg_str": 'Tasiilaq',
        "reg_abv": 'Tasiilaq',
        "reg_lats": np.array([65.6]),
        "reg_lons": np.array([322.37])
    }

    return reg_dict
    
def Umiujaq():
    reg_dict = {
        "reg_str": 'Umiujaq',
        "reg_abv": 'Umiujaq',
        "reg_lats": np.array([56.3]),
        "reg_lons": np.array([283.7])
    }
    
    return reg_dict
    

def UpperAlaska():
    reg_dict = {
        "reg_str": 'Upper Alaska',
        "reg_abv": 'UpperAlaska',
        "reg_lats": np.array([68, 72]),
        "reg_lons": np.array([190, 220])
    }

    return reg_dict
    
def UpperQuebec():
    reg_dict = {
        "reg_str": 'Upper Quebec',
        "reg_abv": 'UpperQuebec',
        "reg_lats": np.array([50.157, 54.829, 62.740, 60.642, 52.213, 50.179]),
        "reg_lons": np.array(
            [293.541, 280.555, 281.859, 295.434, 304.402, 300.139])
    }

    return reg_dict
    
def WholeArctic():
    reg_dict = {
        "reg_str": 'Arctic',
        "reg_abv": 'Arctic',
        "reg_lats": np.array([50, 90]),
        "reg_lons": np.array([0, 360])
    }
    
    return reg_dict




########################################################################
####  HELPER FUNCTIONS
########################################################################

def west180_to_360(west180):
    ''' Convert from deg 180 to deg 360 '''
    east360 = west180 % 360 #wrap to 360 degrees

    return east360