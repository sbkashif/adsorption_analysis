# Salman Bin Kashif
# Sarupria Group
# Clemson University

#Date created: January 18, 2020
#Last modified: February 11, 2020

import MDAnalysis as mda
import numpy as np
import math
import statistics
import sys
import argparse
import os

def create_parser():
        parser = argparse.ArgumentParser(prog = 'HeatMapUsingMDAnalysis', usage = '%(prog)s [-h for help]', \
                                      description = 'Generate the heat map for a trajectory')
        parser.add_argument('-f', "--f", help = 'Input xtc file (Required).')
        parser.add_argument('-s',"--s",help='Input pdb file (required).')
        parser.add_argument('-bt',"--bt",help="Begin Frame (Frame ID)",default=int(0))
        parser.add_argument('-et',"--et",help="End Frame (Frame ID)")
        parser.add_argument('-length',"--length")
        parser.add_argument('-coverage',"--coverage")
        parser.add_argument('-o',"--o")        	
 
        return parser


def main(args):
        XTCFile=args.f
        GROFile=args.s
        OUTFile=args.o
        bt=int(args.bt)
        et=int(args.et)
        length=float(args.length)
        coverage=int(args.coverage)
        data=np.loadtxt("Adsorption_points.xvg")
        getTrajectoryDistribution(bt,et,GROFile,XTCFile,data[0],data[1],data[2],data[3],length,data[4],coverage,OUTFile)


def obtainTrajectoryData(bt,et,grofile,xtcfile):
        global coordinates_pa
        global coordinates_mbl
        global coordinates_mbl_sep
        global size
        global box
        global x_n_bins
        global y_n_bins
        global Lx
        global Ly
        global Lz
        global Lz_max
        global logfile
        global n_molecules
        global bf
        global ef

        logfile=open("output.log","w")	
        #Read the full trajectory from the xtc file
        u=mda.Universe(grofile,xtcfile)
        print(u)

        #Assign the polyamide atom coordinates 'pa' variable
        pa=u.select_atoms("resname MPD1 MPD2 TMC1 TMC2 TMC3")
        logfile.write("\nRead input trajectory Polyamdie coordinates:\n%s"%(pa))
        
        #Assign the methylene blue atom coordinates to 'mbl' variable
        mbl=u.select_atoms("resname MBL")
        logfile.write("\nRead input trajectory Methylene Blue coordinates:\n%s"%(mbl))
        coordinates_mbl=[]
        coordinates_pa=[]
        box=[]

        freq=[]
        for ts in u.trajectory[0:2]:
            freq.append(u.trajectory.time)
        
        time_diff=freq[1]-freq[0];

        print(time_diff)

        bf=(int)((bt-freq[0])/time_diff)
        ef=(int)((et-freq[0])/time_diff)
        print("bf:",bf,"ef:",ef)

        print("Reading trajectory...")
        for ts in u.trajectory[bf:ef+1]:
            logfile.write("\nReading timestep:%s\n"%(ts))
            coordinates_mbl.append(mbl.positions)
            coordinates_pa.append(pa.positions)
            box.append(ts.dimensions[:3])
            #print(ts.dimensions[:3])
        
        coordinates_pa=np.array(coordinates_pa)
        coordinates_mbl=np.array(coordinates_mbl)
        

        box=np.array(box)
    	
        #Since MDAnalysis give output in Angstorm, the units are converted to nm for consistency
        coordinates_pa=np.divide(coordinates_pa,10.0)
        
        coordinates_mbl=np.divide(coordinates_mbl,10.0)

        coordinates_mbl_sep=[[]]*len(coordinates_mbl)


        n_molecules=(int)(len(coordinates_mbl[0])/38)

        print(n_molecules)        

        for i in range(0,len(coordinates_mbl)):
            placeholder=np.split(coordinates_mbl[i],n_molecules)
            coordinates_mbl_sep[i]=placeholder
                   
        print(np.shape(coordinates_mbl_sep))
        print(np.shape(coordinates_pa))
 
        box=np.divide(box,10.0)
       
        #print("COG from MD analysis:\n")
        #print(cog)         

        #Extracting the box size thorughout the trajectory
        Lx=box[:,0]
        Ly=box[:,1]
        Lz=box[:,2]
        Lx_max=max(Lx)
        Ly_max=max(Ly)
        Lz_max=max(Lz)

        print("LX:",Lx)
        size=len(coordinates_mbl)
        print("\nNumber of frames read:%s\n"%(size))
        logfile.write("\nNumber of frames read:%s\n"%(size))
        #The use of this part is in heatmap and it will be revoked later

        #print (coordinates_mbl_sep[2][0][:,2]);        
 
        return coordinates_pa, coordinates_mbl_sep, box

def calculateFrameCOG(frame, input_group_coordinates):
	
	#Obtaining the x,y and z coordinates for the given frame
        #The coordiates are extracted to treat raw or with PBC
        group_coordinates=input_group_coordinates.copy()
        x=input_group_coordinates[frame][:,0]
        y=input_group_coordinates[frame][:,1]
        z=input_group_coordinates[frame][:,2]
        p_x=group_coordinates[frame][:,0]
        p_y=group_coordinates[frame][:,1]
        p_z=group_coordinates[frame][:,2]
 
        #print("readFrame X coordinates read:\n",x)
        #print("readFrame Y coordinates read:\n",y)
        #Lx=box[:,0]
        #Ly=box[:,1]
        #Lz=box[:,2]
        
        #print("Box dimension of frame: %s\t is %s,%s,%s"%(frame,Lx[frame],Ly[frame],Lz[frame]))
        

        #Applying the pbc to p_* variables 
        for i in range(0,len(x)):
                #Uncomment next line to print the x-coordinate value before applying the PBC
                #print(x[i])
                data=p_x[i]/Lx[frame]
                fd=math.floor(data)
                p_x[i]=p_x[i]-(fd*Lx[frame])
        for i in range(0,len(y)):
                #Uncomment next line to print the y-coordinate value before the applyinf the PBC
                #print(y[i])
                data=p_y[i]/Ly[frame]
                fd=math.floor(data)
                p_y[i]=p_y[i]-(fd*Ly[frame])
        for i in range(0,len(z)):
                data=p_z[i]/Lz[frame]
                fd=math.floor(data)
                p_z[i]=p_z[i]-(fd*Lz[frame])
       
        if(frame==0 or (frame+1)%10==0 or (frame+1)==size):
            logfile.write("\ncalculateFrameCOG function X coordinates read:\n%s"%(x))
            logfile.write("\ncalculateFrameCOG function Y coordinates read:\n%s"%(y))
        
        
        #COG calculated for the raw trajectory
        cog_x=sum(x)/len(x)
        cog_y=sum(y)/len(y)
        cog_z=sum(z)/len(z)
        
        #COG converted to PBC
        data_x=cog_x/Lx[frame]
        fd=math.floor(data_x)
        cog_x=cog_x-(fd*Lx[frame])

        data_y=cog_y/Ly[frame]
        fd=math.floor(data_y)
        cog_y=cog_y-(fd*Ly[frame])

        data_z=cog_z/Lz[frame]
        fd=math.floor(data_z)
        cog_z=cog_z-(fd*Lz[frame])

        return cog_x, cog_y,cog_z

def calculate_frame_atoms(frame,molecule,input_group_coordinates):
        
        group_coordinates=input_group_coordinates.copy()
        p_z=group_coordinates[frame][molecule][:,2]

        for i in range(0,len(p_z)):
            data=p_z[i]/Lz[frame]
            fd=math.floor(data)
            p_z[i]=p_z[i]-(fd*Lz[frame])
       
        if(frame==0 or (frame+1)%10==0 or (frame+1)==size):
            logfile.write("\ncalculateFrameMolecule function Z coordinates read:\n%s"%(p_z))


        return p_z




def getTrajectoryDistribution(bt,et,grofile,xtcfile,inflexion1,inflexion2,distance1,distance2,length,mid_point,coverage,outfile):
        
        #global coglist_pa
        #global coglist_mbl	
        obtainTrajectoryData(bt,et,grofile,xtcfile)	
        print(ef)
        print(bf)
        coglist_pa=np.full((ef,1),-100.0)
        count_mbl_adsorbed=np.full((ef,1),0)
        #freq_cog_pa=np.full((x_n_bins,y_n_bins),0)
        #freq_cog_mbl=np.full((x_n_bins,y_n_bins),0)
        
        #bf=bf-1
 

        mid_point=mid_point+Lz_max/2
        inflexion1=inflexion1+Lz_max/2
        inflexion2=inflexion2+Lz_max/2
        
        for i in range(bf,ef):
            i=i-bf
            print ("Looping in frame corresponding to extracted trajectory timestep:%d\n"%(i))
            #Calculating the MBL cog for frame i
	   
        #Calculating the PA cog for frame i
            logfile.write("\nCalculating COG for PA for frame corresponding to extracted trajectory timestep:%s\n"%(i))
            #coglist_pa[i+bf]=cog_z

            #Calculating number of adsorbed atoms in frame i
            logfile.write("Writing the number of adsorbed atoms in each molecule\n")
            for j in range(0,n_molecules):
                a=calculate_frame_atoms(i,j,coordinates_mbl_sep)
                #count = sum(map(lambda x : abs(x-inflexion1) < (length+distance1) and x>mid_point, a))
                count = sum(map(lambda x : ((abs(x-inflexion1) < (length+distance1) and x>mid_point) or (abs(inflexion2-x) <(length+distance2) and x<mid_point)), a))
                logfile.write("Molecule%d \t Adsorpbed atoms%d\n"%(j,count))
                if count>coverage:
                    count_mbl_adsorbed[i+bf]+=1 
            
	   
 

        coglist_w=open(outfile,"w")
        coglist_w.write('#Frame\t MBL count\n')
        
        for i in range(bf,ef):
            coglist_w.write('%d\t%8.3f\n'%((int)(i+1),count_mbl_adsorbed[i]))
	   
if __name__=='__main__':
	args = create_parser().parse_args()
	main(args)
    
