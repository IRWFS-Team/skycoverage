#Author: Eric Shroe
#Purpose: To produce a zoomed in plot of a 1 degree RA by 1 degree DEC region of the sky (uses slow method outlined in "Method Writeup.txt")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import os


plt.close("all")
###Define System Parameters###
M_min = 15 #Limiting Mangitude
band = "J" #Wavelength Band
theta_0 = 1./(60.)/(3.) #isoplanetic angle in degrees
breakdown = 15 #size of individual patches of sky in degrees RA by degrees DEC
xbinsize = 1.#size of RA bin in degrees
ybinsize = 1.#size of DEC bin in degrees

#create arrays of RA and DEC values at the edge of each patch
DECs = np.arange(-90,90,breakdown)
RAs = np.arange(0,360,breakdown)

#determine number of bins in each patch of sky
N_xbins = breakdown/xbinsize
N_ybins = breakdown/ybinsize

zoom = [-28,52] #coordinates of lower left point of the zoomed in region (DEC,RA) 
coords = np.array([250.423458 ,  36.461306]) #coordinates of certain objects in zoomed in region (optional)

plot_coords =False #flag to add the location of certain objects in the zoomed in region to plot

#determine the folder the data will be stored in
outfolder = band+"-"+str(M_min)+"_"+str(int(theta_0*3600))+" slow"
if not os.path.exists(outfolder): #if folder does not exist, create a new one
        print outfolder
        os.mkdir(outfolder)
	os.mkdir(outfolder+"/figures")
        os.mkdir(outfolder+"/datafiles")
else: #otherwise send a warning that this program may overwrite data
        print "WARNING."
        print "Path already exists. May overwrite data."
        #plt.pause(10)

#determine side length of each isoplanatic area (in degrees)
theta_0 = theta_0/np.sqrt(2.) 

#loop through every patch which contains the zoomed in coordinates
for d in [DECs[max(np.where(DECs<=zoom[0])[0])]]:
	for a in [RAs[max(np.where(RAs<=zoom[1])[0])]]:
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
                        N = int(xbinsize/(theta_0/length))*N_xbins +N_xbins #total Number of isoplanatic cells at said DEC (# of cells per bin * # of bins)
                        RA.append(np.linspace(a,a+breakdown,N+1)) #add a list of all (RA) cell boundaries in that bin
                        #loop through every RA bin at that DEC
                        for j in range(int(N_xbins)):
                                #add empty grid of isoplantic regions
                                bins[i-1].append(np.zeros([int(ybinsize/theta_0),N/N_xbins],dtype=bool))


                #loop through every magnitude 
                for i in range(-4,M_min):
                        #open the file containing every star in that patch with a magnitude between i and i+1
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
                                r = np.array([float(data[0]),float(data[1])]) #(RA,DEC)
                                if r[0]<zoom[1] or r[0]>=zoom[1]+1 or  r[1]<zoom[0] or r[1]>=zoom[0]+1: #if the star isn't in the zoomed in region, continue to the next one
                                    continue
                                
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
                                    

                ###Create a plot of that zoomed in region###
                plt.figure()
                plt.imshow(bins[y][x],origin="lower",cmap=plt.get_cmap("gray"))
                if plot_coords:
                    ###plot the location of certain objects in that bin
                    plt.plot((coords[0]-zoom[1])*len(bins[y][x][0,:]),(coords[1]-zoom[0])*len(bins[y][x][:,0]),"r.")
                    plt.xlim(-0.5,len(bins[y][x][0,:])-0.5)
                    plt.ylim(-0.5,len(bins[y][x][:,0])-0.5)
                plt.xticks(np.linspace(-0.5,len(bins[y][x][0,:])-0.5,6),np.linspace(zoom[1],zoom[1]+1,6))
                plt.yticks(np.linspace(-0.5,len(bins[y][x][:,0])-0.5,6),np.linspace(zoom[0],zoom[0]+1,6)
                plt.xlabel("RA [degrees]")
                plt.ylabel("DEC [degrees]")
                #save files
                plt.savefig(outfolder+"/figures/zoom_in_r"+str(zoom[1])+"_d"+str(zoom[0])+".png")
                np.savetxt(outfolder+"/datafiles/zoom_in_r"+str(zoom[1])+"_d"+str(zoom[0])+".txt",np.transpose(bins[y][x]))
#Show plot
plt.show()
plt.close("all")

