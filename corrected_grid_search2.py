#Author: Eric Shore
#Purpose to determine the fration of sky observable with an IR WFS (Slow Method)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import os


plt.close("all")
###Define System Parameters###
M_min = 15 #Limiting Magnitude
band = "J" #Wavelength Band
theta_0 = 1./(60.)/(3.) #isoplanetic angle
breakdown = 15  #size of individual patches of sky in degrees RA by degrees DEC
xbinsize = 1.#size of RA bin in degrees
ybinsize = 1.#size of DEC bin in degrees

#create arrays of RA and DEC values at the edge of each patch
DECs = np.arange(-90,90,breakdown)
RAs = np.arange(0,360,breakdown)

#determine number of bins in each patch of sky
N_xbins = breakdown/xbinsize
N_ybins = breakdown/ybinsize

#determine the folder the data will be stored in
outfolder = band+"-"+str(M_min)+"_"+str(int(theta_0*3600))+" slow"
if not os.path.exists(outfolder): #if folder does not exist, create a new one
        print outfolder
        os.mkdir(outfolder)
	os.mkdir(outfolder+"/figures")
        os.mkdir(outfolder+"/datafiles")
else:  #otherwise send a warning that this program may overwrite data
        print "WARNING."
        print "Path already exists. May overwrite data."
        plt.pause(10)

#determine side length of each isoplanatic area (in degrees)
theta_0 = theta_0/np.sqrt(2.)
#loop through every patch of sky
for d in DECs:
	for a in RAs:
		#determine file path of data
		path = "Dec"+str(d)+"-"+str(d+breakdown)+"_RA"+str(a)+"-"+str(a+breakdown)+"/" 

                ###Create Grid###
                
                DEC = np.linspace(d,d+breakdown,int(ybinsize/theta_0)*N_ybins+1)#DEC gridlines in degrees
                RA = []
                filled = []
                
                y_bins = np.arange(d,d+breakdown+ybinsize,ybinsize)#location of the edges of the bins
                bins = []
                #loop through every DEC bin to determine the RA bins at those DECs
                for i in range(1,len(y_bins)):
                        bins.append([])
                        length  = max(np.cos(y_bins[i]*np.pi/180.),np.cos((y_bins[i-1])*np.pi/180.)) #length of 1 degree of RA at DEC y_bins[i] (relative to DEC=0). chosen to be the maximum of the length at both edges
                        N = int(xbinsize/(theta_0/length))*N_xbins +N_xbins #total Number of theta_0 cells at said DEC (# of cells per bin * # of bins)
                        RA.append(np.linspace(a,a+breakdown,N+1)) #add a list of all (RA) cell boundaries in that bin
                        #loop through every RA bin at that DEC
                        for j in range(int(N_xbins)):
                                #add empty grid of isoplantic regions
                                bins[i-1].append(np.zeros([int(ybinsize/theta_0),N/N_xbins],dtype=bool))


                #loop through every magnitude 
                for i in range(-4,M_min):
	
                        if i==-4:
                                dfile = open("../2MASS/sorted2/"+path+band+"_<-3.txt","r")
                        elif M_min>17 and i>=17:
                                dfile = open("../2MASS/sorted2/"+path+band+"_>17.txt","r")
                        else:
                                dfile = open("../2MASS/sorted2/"+path+band+"_"+str(i+1)+"-"+str(i)+".txt","r")
                        #read the file into a list
                        lines = dfile.readlines()
                        #loop through every line
                        for j in range(len(lines)):
                                #recover the information about RA and DEC of that star
                                data = lines[j].split("|")
                                r = np.array([float(data[0]),float(data[1])])
                               
                                ###determine which bin the star is in###
                                y = int((r[1]-d)//ybinsize)# which DEC bin is the star in?
                                x = int((r[0]-a)//xbinsize)# which RA bin is the star in?
                                if y== breakdown/ybinsize: #deal with edge case
                                        y-=1
                                if x== breakdown/xbinsize:#deal with edge case
                                        x-=1

                                ###determine which cell in that bin the star is located at###
                                dec2 = max(np.where(DEC<=r[1])[0]) #largest DEC value less than star's DEC
                                dec = dec2 - (y)*len(bins[y][x][:,0])#index of star's DEC in its bin (above - #of values in all previous bins )
                                ra2 = max(np.where(RA[y]<=r[0])[0])#largest RA value less than star's RA
                                ra = ra2 - (x)*len(bins[y][x][0,:])#index pf star's RA in its bin (above - #of values in all previous bins )
                                
                                #fill in that grid cell
                                if not bins[y][x][dec][ra]:
                                        bins[y][x][dec][ra]=True
			

                ###bin the data###
                binned = np.zeros([N_xbins,N_ybins])#empty array to store binned data
                #loop through every bin
                for i in range(int(N_ybins)):
                        for j in range(int(N_xbins)):
                                #determine fraction of that bin with stars
                                binned[j,i] = float(np.sum(bins[i][j]))/(len(bins[i][j][0,:])*len(bins[i][j][:,0]))

                ###plot the binned data for that patch of sky###
                plt.figure()
                plt.imshow(np.transpose(binned),origin="lower",)
                cb = plt.colorbar(ticks=np.linspace(0,1,11))
                cb.ax.set_xlim(0,1)
                cb.set_clim(vmin=0,vmax=1)
                
                plt.xticks(np.linspace(-0.5,breakdown/xbinsize-0.5,breakdown/3+1),np.linspace(a,a+breakdown,breakdown/3+1))
                plt.yticks(np.linspace(-0.5,breakdown/ybinsize-0.5,breakdown/3+1),np.linspace(d,d+breakdown,breakdown/3+1))
                #save data
                plt.savefig(outfolder+"/figures/"+path[:-1]+"_seeing.png")
                np.savetxt(outfolder+"/datafiles/"+path[:-1]+"_seeing.txt",np.transpose(binned))
               
plt.close("all")

