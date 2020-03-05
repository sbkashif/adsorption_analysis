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
    parser.add_argument('-step','--step',default=int(float(1)))

    return parser


def main(args):
    filename=args.f
    outfilename=args.o
    step=int(float(args.step))
    readData(filename,outfilename, step)
					

def readData(input_file,output_file,step):
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
    inflection=[]
    second_derivative=[]
    coordinate_x=[]
    coordinate_y=[]
    surface2=[]

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
    for i in range(0,len(x)-step,step):
        yd=y[i+step]-y[i]
        xd=x[i+step]-x[i]
        slope=yd/xd
        slopes.append(slope)
        coordinate_x.append(x[i])
        coordinate_y.append(y[i])
        angle = np.rad2deg(np.arctan(slope))
        angles.append(angle)

    print("Trimmed size", len(angles))
    for i in range(1,len(angles)):
        print(coordinate_x[i],angles[i])
        if (angles[i]>(-15) and angles[i-1]<-15):
            surface1=coordinate_x[i]
        if (angles[i]>(15) and angles[i-1]<15):
            surface2.append(coordinate_x[i])
    for i in range (0, len(coordinate_x)-1):
        sd=(coordinate_y[i+1]-2*coordinate_y[i]+coordinate_y[i-1])/((coordinate_x[i+1]-coordinate_x[i])*(coordinate_x[i+1]-coordinate_x[i]))
        #sd=(y[i+step]-2*y[i]+y[i-step])/((x[i+step]-x[i])*(x[i+step]-x[i]))
        print(coordinate_x[i],sd)
        second_derivative.append(sd)
    print(len(second_derivative))
    for i in range(0, len(second_derivative)-1):
        if(((np.sign(second_derivative[i])==-1 and np.sign(second_derivative[i+1])==1) or (np.sign(second_derivative[i])==1 and np.sign(second_derivative[i+1])==-1)) or np.isclose(second_derivative[i],0.0)==True and np.isclose(coordinate_y[i],0.0)==False):
            inflection.append(coordinate_x[i])



    out=open(output_file,"w")
    out.write("#Surface point1:\n%s"%(surface1))
    out.write('\n')
    out.write("#Surface point2:\n%s"%(surface2[0]))
    out.write('\n') 
    out.write("#Sign change second derivative points:\n")
    for i in inflection:
        out.write("%f\n"%(i))
    out.close()

    print("Output file written")
	




	
if __name__=='__main__':

    args = create_parser().parse_args()
    main(args)



