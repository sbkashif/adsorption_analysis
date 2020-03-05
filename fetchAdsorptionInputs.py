import numpy as np

import statistics
import sys
import argparse
import os



def readData(input_file):
    global x
    global y
    global z
#check if input file exists
    count=0
    c=0
    
    numerics=[]
    
    if os.path.isfile(input_file):
        #print("Input file exists")
        f=open(input_file,'r')
        lines=f.readlines()
        f.close()
    else:
        print("Input file does not exist.Exiting..")
        sys.exit(0)
    for i in lines:
        if (not i.startswith('#') and not i.startswith('@')):
            numerics.append(i)
    #print("Read data from file:",lines)
    return numerics
    
def fetchSurfacePoints(lines):
    surface_point1=(float)(lines[0])
    surface_point2=(float)(lines[1])
    
    return surface_point1, surface_point2

def fetchInflectionPoints(lines,surface_point1,surface_point2):
    inflection_point1=(float)(lines[-1])
    difference1=surface_point1-inflection_point1
    inflection_point2=(float)(lines[2])
    mid_point=(inflection_point1+inflection_point2)/2
    difference2=inflection_point2-surface_point2
    out=open("Adsorption_points.xvg","w")
    out.write("%f\n%f\n%f\n%f\n%f"%(inflection_point1,inflection_point2,difference1,difference2,mid_point))
    return inflection_point1,inflection_point2,difference1,difference2,mid_point

surface_file='surface.xvg'
inflection_file='inflection.xvg'
fetchInflectionPoints(readData(inflection_file),fetchSurfacePoints(readData(surface_file))[0],fetchSurfacePoints(readData(surface_file))[1])
