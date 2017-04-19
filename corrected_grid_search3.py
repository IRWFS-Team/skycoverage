#Author: Eric Shore
#Purpose to determine the fration of sky observable with an IR WFS (Fast Method)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import os

plt.close("all")
###Define System Parameters###
M_min = 15 #Limiting Magnitude
band = "J" #Wavelength Band
theta_0 = 1./(60.)/(3.) #isoplanetic angle
breakdown = 15 #size of individual patches of sky in degrees RA by degrees DEC
xbinsize = 1.#size of RA bin in degrees
ybinsize = 1.#size of DEC bin in degrees

#create arrays of RA and DEC values at the edge of each patch
DECs = np.arange(-90,90,breakdown)
RAs = np.arange(0,360,breakdown)


#determine number of bins in each patch of sky
N_xbins = breakdown/xbinsize
N_ybins = breakdown/ybinsize

#determine the folder the data will be stored in
outfolder = band+"-"+str(M_min)+"_"+str(int(theta_0*3600))+" fast"
if not os.path.exists(outfolder): #if folder does not exist, create a new one
        print outfolder
        os.mkdir(outfolder)
        os.mkdir(outfolder+"/figures")
        os.mkdir(outfolder+"/datafiles")
else:  #otherwise send a warning that this program may overwrite data
        print "WARNING."
        print "Path already exists. May overwrite data."
        plt.pause(10) #gives time to stop program

#determine side length of each isoplanatic area (in degrees)
theta_0 = theta_0/np.sqrt(2.)
#loop through every patch of sky
for d in DECs:
	for a in RAs:
                #determine file path of data
		path = "Dec"+str(d)+"-"+str(d+breakdown)+"_RA"+str(a)+"-"+str(a+breakdown)+"/" 

                ###Create Grid###
               
                filled = []
                y_bins = np.arange(d,d+breakdown+ybinsize,ybinsize)#location of the edges of the bins
                L_cells = np.zeros(len(y_bins)-1) #empty array to hold the length (RA extent) of each isoplanatic area in degrees 
                bins = []
                #loop through every DEC bin to determine the RA bins at those DECs
                for i in range(1,len(y_bins)):
                        bins.append([])
                        length  = max(np.cos(y_bins[i]*np.pi/180.),np.cos((y_bins[i-1])*np.pi/180.)) #length of 1 degree of RA at DEC y_bins[i] (relative to DEC=0).  chosen to be the maximum of the length at both edges
                        N = int(xbinsize/(theta_0/length))*N_xbins +N_xbins #total Number of theta_0 cells at said DEC (# of cells per bin * # of bins)
                        #determine the length (RA extent) of each isoplanatic area in degrees  
                        L_cells[i-1] = breakdown/(int(N/N_xbins)*N_xbins)
                        #loop through every RA bin at that DEC
                        for j in range(int(N_xbins)):
                                #add empty grid of isoplantic regions
                                bins[i-1].append(np.zeros([int(ybinsize/theta_0),N/N_xbins],dtype=bool))

                #loop through every magnitude 
                for i in range(-4,M_min):
                         #open the file containing every star in that patch with a magnitude between i and i+1
                        if i==-4:

                                infile = "../2MASS/sorted3/"+path+band+"_<-3.txt"
                        elif M_min>17 and i>=17:

                                infile = "../2MASS/sorted3/"+path+band+"_>17.txt"
                        else:

                                infile = "../2MASS/sorted3/"+path+band+"_"+str(i+1)+"-"+str(i)+".txt"
                        #open the file and store as an array
                        data = np.loadtxt(infile,unpack=True)
                        if len(data)==0: # if the file is empty move to next one
                                continue
                        #recover the RA and DEC of every star
                        r = np.array([data[0],data[1]])
                        
                        ###determine which bin every star is in###       
                        y = np.floor((r[1]-d)/ybinsize)
                        x = np.floor((r[0]-a)/xbinsize)
                        y=y.astype(int)
                        x=x.astype(int)
                        
                        ###determine which cell in its bin every star is located at###
                        dec_bin_length = ybinsize/int(ybinsize/theta_0) #the height (DEC extent) in degrees of each isoplanatic area 
                        dec2 = np.floor((r[1]-d)/dec_bin_length)
                        ra2 = np.floor((r[0]-a)/(L_cells[y]))
                        dec = dec2 - y*int(ybinsize/theta_0)
                        ra = ra2 - x*(xbinsize/L_cells[y])
                        
                        
                        ###fill in grid cells###
                        if np.size(y)==1: #if there is only 1 star
                                bins[int(y)][int(x)][int(dec)][int(ra)]=True
                        else: #if there is more than 1 star
                                for j in range(np.size(y)):
                                        bins[int(y[j])][int(x[j])][int(dec[j])][int(ra[j])]=True
			

                ###bin the data###
                binned = np.zeros([N_xbins,N_ybins]) #empty array to store binned data
                #loop through every bin
                for i in range(int(N_ybins)):
                        for j in range(int(N_xbins)):
                                #determine fraction of that bin with stars
                                binned[j,i] = float(np.sum(bins[i][j]))/(len(bins[i][j][0,:])*len(bins[i][j][:,0]))
                                
                ###plot the binned data for that patch of sky###
                plt.figure()
                plt.imshow(np.transpose(binned),origin="lower",interpolation="nearest")
                cb = plt.colorbar()

                plt.xticks(np.linspace(-0.5,breakdown/xbinsize-0.5,breakdown/3+1),np.linspace(a,a+breakdown,breakdown/3+1))
                plt.yticks(np.linspace(-0.5,breakdown/ybinsize-0.5,breakdown/3+1),np.linspace(d,d+breakdown,breakdown/3+1))
                #save data
                plt.savefig(outfolder+"/figures/"+path[:-1]+"_seeing.png")
                np.savetxt(outfolder+"/datafiles/"+path[:-1]+"_seeing.txt",np.transpose(binned))

plt.close("all")


