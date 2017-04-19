#Author: Eric Shore
#Purpose: To plot a zoomed in region of every Globular Cluster in the file "Globular_Clusters.txt"
import numpy as np
import matplotlib.pyplot as plt
import os

#Function for converting RA and DEC from hms and dms to degrees
def string_to_degrees(ra_string,dec_string):
    ###Convert RA to degrees###
    #Split the ra coordinate into  h, m, and s
    ra_split1 = ra_string.split("h")
    ra_split2 = ra_split1[1].split("m")
    h = float(ra_split1[0])
    m = float(ra_split2[0])
    s = float(ra_split2[1].split("s")[0])
    #convert to degrees
    ra_degrees = h*15.+ m/4. + s/240.
    
    ###Convert DEC to degrees###
    #Split the dec coordinate into  d, m, and s
    dec_split1 = dec_string.split("d")
    dec_split2 = dec_split1[1].split("m")
    d = float(dec_split1[0])
    m = float(dec_split2[0])
    s = float(dec_split2[1].split("s")[0])
    #convert to degrees
    dec_degrees = -100. #meaningless value just in case
    if dec_string[0]=="-": #if DEC is negative
        dec_degrees = d-m/60.-s/3600.
    else: #if DEC is positive
        dec_degrees = d+ m/60. + s/3600.
    return [ra_degrees,dec_degrees]
 
#main function   
def main():
    ###Define System Parameters###
    M_min = 15 #limitng magnitude
    band = "J" #wavelength band
    theta_0 = 1./(60.)/(3.) #isoplanetic angle
    breakdown = 15 #size of individual patches of sky in degrees RA by degrees DEC
    xbinsize = 1.#size of RA bin in degrees
    ybinsize = 1.#size of DEC bin in degrees

    #create arrays of RA and DEC values at the edge of each patch
    DECs = np.arange(-90,90,breakdown)
    RAs = np.arange(0,360,breakdown)

    #create arrays of RA and DEC values at the edge of each patch
    N_xbins = breakdown/xbinsize
    N_ybins = breakdown/ybinsize
    
    #determine the folder the data is stored in
    outfolder = band+"-"+str(M_min)+"_"+str(int(theta_0*3600))+" slow"
    #determine side length of each isoplanatic area (in degrees)
    theta_0 = theta_0/np.sqrt(2.)

    plot_coords=True  #flag to add the location of certain objects in the zoomed in region to plot
    #open file containing location of Globular clusters
    fout = open("Globular_Clusters.txt","r")
    #read 1st 2 lines (they don't contain any data)
    fout.readline()
    fout.readline()
    #loop through every Globular cluster
    while True:
        #read in line
        data = fout.readline()
        if data=="": #if the line is empty, this means that it is at the end of the file, so exit loop
            break
        #recover Cluster name, RA, DEC, and whether or not the Cluster is deemed "visible" from the data
        name,ra_string,dec_string,visible = data.split("\t")
        remove_spaces =  name.split(" ")
        if remove_spaces[1]==0:
            name = remove_spaces[0]
        else:
            name = remove_spaces[0]+" "+remove_spaces[1]

        #recover RA and DEC in degrees
        ra,dec = string_to_degrees(ra_string,dec_string)

        #determine zoomed in region
        zoom = [np.floor(dec),np.floor(ra)]
        coords = [ra,dec]

    
        #loop through every patch which contains the zoomed in coordinates
        for d in [DECs[max(np.where(DECs<=zoom[0])[0])]]:
            for a in [RAs[max(np.where(RAs<=zoom[1])[0])]]:
                #determine file path of data
		path = "Dec"+str(d)+"-"+str(d+breakdown)+"_RA"+str(a)+"-"+str(a+breakdown)+"/" 

                ###Create Grid###
                
                DEC = np.linspace(d,d+breakdown,int(ybinsize/theta_0)*N_ybins+1)#DEC gridlines in degrees
                RA = []
                
                y_bins = np.arange(d,d+breakdown+ybinsize,ybinsize)#location of the edges of the bins
                bins = []
                #loop through every DEC bin to determine the RA bins at those DECs
                for i in range(1,len(y_bins)):
                        bins.append([])
                        length  = max(np.cos(y_bins[i]*np.pi/180.),np.cos((y_bins[i-1])*np.pi/180.)) #length of 1 degree of RA at DEC y_bins[i] (relative to DEC=0). chosen to be the maximum of the length at both edges
                        N = int(xbinsize/(theta_0/length))*N_xbins +N_xbins #total Number of theta_0 cells at said DEC (# of cells per bin * # of bins)
                        RA.append(np.linspace(a,a+breakdown,N+1)) #add a list of all (RA) cell boundaries in that bin
                        #loop through every RA bin at that DE
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
                                r = np.array([float(data[0]),float(data[1])])
                                if r[0]<zoom[1] or r[0]>=zoom[1]+1 or  r[1]<zoom[0] or r[1]>=zoom[0]+1: #if the star isn't in the zoomed in region, continue to the next one
                                    continue

                                ###determine which bin the star is in###
                                y = int((r[1]-d)//ybinsize)# which DEC bin is the star in?
                                x = int((r[0]-a)//xbinsize)# which RA bin is the star in?
                                if y== breakdown/ybinsize: #deal with edge case
                                    y-=1
                                if x== breakdown/xbinsize:#deal with edge case
                                    x-=1

                                ###determine which cell in that bin the star is star is located at###
                                dec2 = max(np.where(DEC<=r[1])[0]) #largest DEC value less than star's DEC
                                dec = dec2 - (y)*len(bins[y][x][:,0])#index of star's DEC in its bin (above - #of values in all previous bins )
                                ra2 = max(np.where(RA[y]<=r[0])[0])#largest RA value less than star's RA
                                ra = ra2 - (x)*len(bins[y][x][0,:])#index pf star's RA in its bin (above - #of values in all previous bins )
                                #fill in that grid cell
                                if not bins[y][x][dec][ra]:
                                    bins[y][x][dec][ra]=True
                ###Create a plot of that zoomed in region###
                plt.figure()
                plt.imshow(bins[y][x],origin="lower",cmap=plt.get_cmap("gray"),interpolation="nearest")
                if plot_coords:
                    ###plot the location of the globular cluster###
                    plt.plot((coords[0]-zoom[1])*len(bins[y][x][0,:]),(coords[1]-zoom[0])*len(bins[y][x][:,0]),"r.")
                    plt.xlim(-0.5,len(bins[y][x][0,:])-0.5)
                    plt.ylim(-0.5,len(bins[y][x][:,0])-0.5)
                    plt.title(name)
                plt.xticks(np.linspace(-0.5,len(bins[y][x][0,:])-0.5,6),np.linspace(zoom[1],zoom[1]+1,6))
                plt.yticks(np.linspace(-0.5,len(bins[y][x][:,0])-0.5,6),np.linspace(zoom[0],zoom[0]+1,6))
                plt.xlabel("RA [degrees]")
                plt.ylabel("DEC [degrees]")
                #save files
                plt.savefig(outfolder+"/Globular_Clusters/"+str(name)+".png")
                np.savetxt(outfolder+"/Globular_Clusters/"+str(name)+".txt",np.transpose(bins[y][x]))


                

        plt.close("all")


main()       
