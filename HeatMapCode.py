# Salman Bin Kashif
# Sarupria Group
# Clemson University

#Date created:- 

import MDAnalysis as mda
import numpy as np
import math
import statistics
import sys
import argparse
import os

def create_parser():
	"""

	This is the function to pass the input values
	
	:type binx: int
	:param binx: No. of bins in X-direction
	
	"""

	parser = argparse.ArgumentParser(prog = 'HeatMapUsingMDAnalysis', usage = '%(prog)s [-h for help]', \
                                      description = 'Generate the heat map for a trajectory')
	parser.add_argument('-f', "--f", help = 'Input xtc file (Required).')
	parser.add_argument('-s',"--s",help='Input pdb file (required).')
	parser.add_argument('-binx', "--binx", help = "Enter x bin size", default=int(0.2))
	parser.add_argument('-biny', "--biny", help = "Enter y bin size", default=int(0.2))
	parser.add_argument('-bf',"--bf",help="Begin Frame (Frame ID)",default=int(1))
	parser.add_argument('-ef',"--ef",help="End Frame (Frame ID)")
	return parser

def main(args):
    XTCFile=args.f  
    GROFile=args.s
    binx=float(args.binx)
    biny=float(args.biny)
    bf=int(args.bf)
    ef=int(args.ef)
    #obtainTrajectoryData(binx,biny,bf,ef)
    getTrajectoryDistribution(binx,biny,bf,ef)


def obtainTrajectoryData(x_bin_size,y_bin_size,bf,ef):
    global coordinates
    global size
    global box
    global x_n_bins_max
    global y_n_bins_max
    global logfile
    global Lx_avg
    global Ly_avg

    logfile=open("output.log","w")	
    #Read the full trajectory from the xtc file
    u=mda.Universe("./npt-equil-whole-pacxb.gro","./npt-equil-whole-pacxb.xtc")
    print(u) 

    #Assign the polyamide atom coordinates 'pa' variable
    pa=u.select_atoms("resname MPD1 MPD2 TMC1 TMC2 TMC3")
    logfile.write("\nRead input trajectory Polyamdie coordinates:\n%s"%(pa))
        
    coordinates=[]
    box=[]

    print(u.trajectory)
        

    for ts in u.trajectory[bf:ef]:
        print(ts)
        logfile.write("\nReading timestep:%s\n"%(ts))
        coordinates.append(pa.positions)
        box.append(pa.dimensions)
        
    coordinates=np.array(coordinates)
    box=np.array(box)

    print("Total number of frames",len(coordinates))
 
    	
    #Since MDAnalysis give output in Angstorm, the units are converted to nm for consistency
    coordinates_pa=np.divide(coordinates,10.0)
    box=np.divide(box,10.0)

    Lx=box[:,0]
    Ly=box[:,1]
    Lz=box[:,2]
    Lx_max=max(Lx)
    Ly_max=max(Ly)

    Lx_avg=sum(Lx)/len(Lx)
    Ly_avg=sum(Ly)/len(Lx)

    x_n_bins_max=math.ceil(Lx_max/x_bin_size)
    y_n_bins_max=math.ceil(Ly_max/y_bin_size)

    print("Maximum number of x bins: %s"%(x_n_bins_max))
    print("Maximum number of y bins: %s"%(y_n_bins_max))

def readFrame(frame):
    global x
    global y
    global z
    global Lx
    global Ly
    global Lz
    global Lx_max
    global Ly_max
    	
	
	#Obtaining the x,y and z coordinated from the frame
    x=coordinates[frame][:,0]
    y=coordinates[frame][:,1]
    z=coordinates[frame][:,2]
    Lx=box[:,0]
    Ly=box[:,1]
    Lz=box[:,2]

    
    print("Box dimension:%s%s%s\n"%(Lx[frame],Ly[frame],Lz[frame]))
    for i in range(0,len(x)):
        data=x[i]/Lx[frame]
        fd=math.floor(data)
        x[i]=x[i]-(fd*Lx[frame])
    for i in range(0,len(y)):
        data=y[i]/Ly[frame]
        fd=math.floor(data)
        y[i]=y[i]-(fd*Ly[frame])
    for i in range(0,len(z)):
        data=z[i]/Lz[frame]
        fd=math.floor(data)
        z[i]=z[i]-(fd*Lz[frame])

def getDistbn(frame,x_bin_size,y_bin_size):
    global max_z
    global min_z
    readFrame(frame)

    x_n_bins=math.ceil(Lx[frame]/x_bin_size)
    y_n_bins=math.ceil(Ly[frame]/y_bin_size)

    print("No. of x bins:%s, No. of y bins:%s"%(x_n_bins,y_n_bins))

    max_z=np.full((x_n_bins,y_n_bins),-100.0)   
    min_z=np.full((x_n_bins,y_n_bins),100.0)
    for i in range(0,len(x)):
        bin_x_ID=math.floor(x[i]/x_bin_size)
        bin_y_ID=math.floor(y[i]/y_bin_size)
        max_z[bin_x_ID][bin_y_ID]=max(max_z[bin_x_ID][bin_y_ID],z[i])
        min_z[bin_x_ID][bin_y_ID]=min(min_z[bin_x_ID][bin_y_ID],z[i])
    
def getTrajectoryDistribution(x_bin_size,y_bin_size,bf,ef):
    global z_maxvalues
    global z_maxvalues_avg
    global z_minvalues  
    global z_minvalues_avg

    ef=int(ef)
    bf=int(bf-1)
    obtainTrajectoryData(x_bin_size,y_bin_size,bf,ef)	
   
    z_maxvalues=np.full((ef-bf,x_n_bins_max,y_n_bins_max),-100.0)
    z_minvalues=np.full((ef-bf,x_n_bins_max,y_n_bins_max),100.0)
    
    print(z_maxvalues.shape)
   
    n_count_max=np.zeros((x_n_bins_max,y_n_bins_max))
    n_count_min=np.zeros((x_n_bins_max,y_n_bins_max))
    z_maxvalues_avg=np.zeros((x_n_bins_max,y_n_bins_max))
    z_minvalues_avg=np.zeros((x_n_bins_max,y_n_bins_max))
    
    for i in range(bf,ef):
        i=i-bf
        print ("Frame%d\n"%(i))
        getDistbn(i,x_bin_size,y_bin_size)
        z_maxvalues[i]=max_z
        z_minvalues[i]=min_z

    for i in range(bf,ef):
        i=i-bf
        for j in range(0,x_n_bins_max):
            for k in range(0,y_n_bins_max):
                if(j*x_bin_size <=Lx[i] and k*y_bin_size <= Ly[i]):
                    n_count_max[j][k]+=1
                    n_count_min[j][k]+=1
    print(n_count_max)
    for i in range(bf,ef):
        i=i-bf
        for j in range(0,x_n_bins_max):
            for k in range(0,y_n_bins_max):
                if z_maxvalues[i][j][k]!=-100:
                    #n_count_max[j][k]+=1
                    z_maxvalues_avg[j][k]+=z_maxvalues[i][j][k]
                if z_minvalues[i][j][k]!=100:
                    #n_count_min[j][k]+=1
                    z_minvalues_avg[j][k]+=z_minvalues[i][j][k]
    print(z_maxvalues_avg)
    for i in range(0,len(z_maxvalues_avg)):
        for j in range(0,len(z_maxvalues_avg[i])):
            #if (i==(len(z_maxvalues_avg)-1)):
            #    z_maxvalues_avg[i][j]=z_maxvalues_avg[i][j]+(x_n_bins_max-Lx_avg/x_bin_size)*z_maxvalues_avg[i][j]
            #if(i!=(len(z_maxvalues_avg)-1) and j==(len(z_maxvalues_avg)-1)):
            #    z_maxvalues[i][j]=z_maxvalues_avg[i][j]+(y_n_bins-Ly_avg/y_bin_size)*z_maxvalues_avg[i][j]
            z_maxvalues_avg[i][j]/=n_count_max[i][j]
            z_minvalues_avg[i][j]/=n_count_min[i][j]
    z_maxvalues_stat_w=open('./MaxZValuesStats.xvg',"w")
    z_minvalues_stat_w=open('./MinZValuesStats.xvg',"w")
    for i in range(0,len(z_maxvalues_avg)-1):
        for j in range(0,len(z_maxvalues_avg[i])-1):
            z_maxvalues_stat_w.write('%8.3f\t%8.3f\t%8.3f\n'%(i*x_bin_size+x_bin_size/2,j*y_bin_size+y_bin_size/2,z_maxvalues_avg[i][j]))
            z_minvalues_stat_w.write('%8.3f\t%8.3f\t%8.3f\n'%(i+x_bin_size/2,j+y_bin_size/2,z_minvalues_avg[i][j]))
        if(i!=len(z_maxvalues_avg)-1):
            z_maxvalues_stat_w.write("\n")
            z_minvalues_stat_w.write("\n")

if __name__=='__main__':
    args = create_parser().parse_args()
    main(args)
    
#obtainTrajectoryData(0.2,0.2,19000,20000)
