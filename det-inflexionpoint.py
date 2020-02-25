import numpy as np

import statistics
import sys
import argparse
import os


def create_parser():

    parser = argparse.ArgumentParser(prog = 'stat_operations', usage = '%(prog)s [-h for help]', \
                                      description = 'Calculate distance between two atoms')
    parser.add_argument('-f', "--f", help = 'Input .xvg file (Required).')
    parser.add_argument('-o','--o',help="Output .xvg file", default='stat_output.xvg')
    parser.add_argument('-start','--start',default=int(float(2)))


    return parser


def main(args):
    filename=args.f
    outfilename=args.o
    start=int(float(args.start))
    readData(filename, outfilename)
					

def readData(input_file,output_file):
    global x
    global y
    global z
    #check if input file exists
    count=0
    
    if os.path.isfile(input_file):
        print("Input file exists")
        f=open(input_file,'r')
        lines=f.readlines()
        f.close()
    else:
        print("Input file does not exist")

    x=[]
    y=[]
    z=[]
    angles=[]
    slopes=[]
    inflexion=[]
    
    for i in lines:
    #if linestart != @ or #
        if (not i.startswith('#') and not i.startswith('@')):
            line=i.split()
            x.append(float(line[0].strip()))
            y.append(float(line[1].strip()))
        else:
            count+=1

    start=count
    end=count+len(x)

    print("Size",len(x))
    for i in range(0,len(x)-1):
        yd=y[i+1]-y[i]
        xd=x[i+1]-x[i]
        slope=yd/xd
        slopes.append(slope)
        angle = np.rad2deg(np.arctan(slope))
        angles.append(angle)
    for i in range(0,len(angles)):
        if (angles[i]==0 and angles[i-1]<0):
            surface=x[i]
        if ((slopes[i]<slopes[i-1])):
            inflexion.append(x[i])


    print("Angles:",angles)
    out=open(output_file,"w")
    out.write("#Surface point:\n%s"%(surface))
    out.write('\n') 
    out.write("#Decreasing slope points:\n")
    for i in inflexion:
        out.write("%f\n"%(i))
    out.close()

    print("Output file written")
	




def readData(input_file,output_file):
    global x
    global y
    global z
    #check if input file exists
    count=0
    
    if os.path.isfile(input_file):
        print("Input file exists")
        f=open(input_file,'r')
        lines=f.readlines()
        f.close()
    else:
        print("Input file does not exist")

    x=[]
    y=[]
    z=[]
    angles=[]
    slopes=[]
    inflexion=[]
    
    for i in lines:
    #if linestart != @ or #
        if (not i.startswith('#') and not i.startswith('@')):
            line=i.split()
            x.append(float(line[0].strip()))
            y.append(float(line[1].strip()))
        else:
            count+=1

    start=count
    end=count+len(x)

    print("Size",len(x))
    for i in range(0,len(x)-1):
        yd=y[i+1]-y[i]
        xd=x[i+1]-x[i]
        slope=yd/xd
        slopes.append(slope)
        angle = np.rad2deg(np.arctan(slope))
        angles.append(angle)
    for i in range(0,len(angles)):
        if (angles[i]==0 and angles[i-1]<0):
            surface=x[i]
        if ((slopes[i]<slopes[i-1])):
            inflexion.append(x[i])


    print("Angles:",angles)
    out=open(output_file,"w")
    out.write("#Surface point:\n%s"%(surface))
    out.write('\n') 
    out.write("#Decreasing slope points:\n")
    for i in inflexion:
        out.write("%f\n"%(i))
    out.close()

    print("Output file written")
	
if __name__=='__main__':

    args = create_parser().parse_args()
    main(args)



