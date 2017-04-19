#Author: Eric Shore
#Purpose: to plot a combined image for sky visibility using both methods outlined in Method_Writeup.txt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from astropy.io import fits
from Globular_Clusters import string_to_degrees
from matplotlib import rc
#Figure text formating
rc('text',usetex=True)
rc('font',family='serif')
rc('font',serif='cm')
rc('font',size=28)


###converts hms (RA) and dms (DEC) into decimal degrees###
def string_to_degrees(ra_string,dec_string):
    ###first handle RA###
    #Split up the string into 3 numbers delimited by "h","m","s" characters 
    ra_split1 = ra_string.split("h") 
    ra_split2 = ra_split1[1].split("m")
    h = float(ra_split1[0])
    m = float(ra_split2[0])
    s = float(ra_split2[1].split("s")[0])
    #calculate value of angle in degrees
    ra_degrees = h*15.+ m/4. + s/240.


    ###Now handle DEC###
    #Split up the string into 3 numbers delimited by "d","m","s" characters 
    dec_split1 = dec_string.split("d")
    dec_split2 = dec_split1[1].split("m")
    d = float(dec_split1[0])
    m = float(dec_split2[0])
    s = float(dec_split2[1].split("s")[0])
    #Calculate value of angle in degrees
    dec_degrees = -100. #meaningless placeholder value
    if dec_string[0]=="-": #subtract minutes and seconds from degrees if d is negative
        dec_degrees = d-m/60.-s/3600.
    else: #add minutes and seconds from degrees if d is negative
        dec_degrees = d+ m/60. + s/3600.
    return [ra_degrees,dec_degrees]


#Close all figures in case some are open
plt.close("all")

###Set Parameters###
xbinsize = 1. #length of bins (in degrees RA)
ybinsize = 1. #height of bins (in degrees DEC)
breakdown = 15 #size of regions that were considered individually
band = "J" #Wavelength band of calculations
theta_0 = "20" #isoplanatic angle
M_min = 15  #limiting magnitude


#initiate matrices to hold the sky visibility calculated using both methods (for the entire sky, not each region)
matrix = np.zeros([int(360/xbinsize),int(180/ybinsize)]) #slow method
matrix2 = np.zeros([int(360/xbinsize),int(180/ybinsize)]) #fast method
#determine folder the files are stored in  
folder = band +"-"+str(M_min)+"_"+str(theta_0)+" slow"
folder2 =  band +"-"+str(M_min)+"_"+str(theta_0)+" fast"

### Fill in the matrices ###
#loop for every region
for i,d in enumerate(np.arange(-90,90,breakdown,dtype=int)): 
    for j,a in enumerate(np.arange(0,360,breakdown,dtype=int)):
	#determine filename of sky visibility data
        path = "Dec"+str(d)+"-"+str(d+breakdown)+"_RA"+str(a)+"-"+str(a+breakdown)
	#determine file paths for both methods
        infile = folder+"/datafiles/"+path+"_seeing.txt"
        infile2 = folder2+"/datafiles/"+path+"_seeing.txt"
	#load data (for the region) for both methods into dummy arrays
        data = np.loadtxt(infile,unpack=True)
        data2 = np.loadtxt(infile2,unpack=True)
	#insert data into the matrices for the entire sky
        matrix[j*breakdown:(j+1)*breakdown,i*breakdown:(i+1)*breakdown] = data
        matrix2[j*breakdown:(j+1)*breakdown,i*breakdown:(i+1)*breakdown] = data2
	#plot the difference between the values determined by each method for a specific region
        if d == -75 and a == 180:
            plt.figure()
            plt.imshow(np.transpose(data-data2),origin="lower",cmap="RdBu")
            cb = plt.colorbar()

### determine positions of all observable and unobservable Globular Clusters in the Galaxy###        
fout = open("Globular_Clusters.txt","r")
fout.readline() #first 2 lines contain no data 
fout.readline()	#first 2 lines contain no data
#create empty lists for positions of clusters
ras1,ras2 = [],[]
decs1,decs2 = [],[]
#Loop through file
while True:
    #read out a line
    data = fout.readline()
    #if line is empty, then we are at the end of the file, so break out of loop
    if data=="":
        break
    #split the line into the requisite columns
    name,ra_string,dec_string,can_see = data.split("\t")
    #format the name into 1 string with no spaces
    remove_spaces =  name.split(" ")
    if remove_spaces[1]==0:
        name = remove_spaces[0]
    else:
        name = remove_spaces[0]+" "+remove_spaces[1]
    #determine the positions of the globular cluster
    ra,dec = string_to_degrees(ra_string,dec_string)
    #sort into correct list depending if the globular cluster is observable.
    if can_see[0]=="y": #observable
        ras1.append(ra)
        decs1.append(dec)
    elif can_see[0]== "n":  #not observable
        ras2.append(ra)
        decs2.append(dec)
#turn the lists into arrays
ras1,decs1 = np.array(ras1),np.array(decs1)
ras2,decs2 = np.array(ras2),np.array(decs2)
#close file
fout.close()


###plot the sky visibility using slow method###
plt.figure()
plt.imshow(np.transpose(np.log10(matrix)),origin="lower",interpolation="nearest") #plots log(visibility)

##follwing two lines overplot the positions of the clusters (commented out)
#plt.plot((ras1-0.5)/360.*len(matrix[:,0]),(decs1+89.5)/180.*len(matrix[0]),"k.")
#plt.plot((ras2-0.5)/360.*len(matrix[:,0]),(decs2+89.5)/180.*len(matrix[0]),".",color=(0.5,0.5,0.5))
##

plt.xlim(-0.5,len(matrix[:,0])-0.5)
plt.ylim(-0.5,len(matrix[0])-0.5)

#create a logorithmic colorbar
cb = plt.colorbar(ticks =np.log10(np.array([0.003,0.01,0.1,0.9])))
cb.set_clim(vmin=-3,vmax=0)
cb.ax.set_yticklabels([0.003,0.01,0.1,0.9])
cb.set_label("Fraction of sky visible [per square degree]")

plt.xticks(np.linspace(-0.5,360/xbinsize-0.5,13),np.arange(0,390,30))
plt.yticks(np.linspace(-0.5,180/ybinsize-0.5,13),np.arange(-90,105,15))

plt.xlabel("RA [deg]")
plt.ylabel("DEC [deg]")

##save figure differently if overplotting globular clusters (commented out)
#plt.savefig(folder+"/Globular_Clusters/composite.png")
#plt.savefig(folder+"/figures/composite.png")

###Save data as a fits file### (commented out)

#hdr = fits.Header()
#hdr["band"] = (band, "the wavelength band ")
#hdr["Mmin"] = (str(M_min), "the dimmest magnitude stars capable of being used for wavefront corrections")
#hdr["theta0"] = (str(theta_0)+" arcseconds","the isoplanatic angle of the Adaptive Optics system")
#hdr["RAbin"] = (str(xbinsize)+" degrees","The length of one bin in RA")
#hdr["DECbin"] = (str(ybinsize)+" degrees","the length of one bin in DEC")
#hdr["patch"] = (str(breakdown**2)+" degrees","the size of each individual patch of sky used to create the composite image")
#hdr["method"] = ("slow", "The method used to dertermine sky coverage")
#fits.writeto("Sky_Coverage.fits",matrix,hdr)



###plot the sky visibility using fast method ###
plt.figure()
plt.imshow(np.transpose(matrix2),origin="lower",interpolation="nearest")
#set up colobar
cb = plt.colorbar(ticks=np.linspace(0,1,11))
cb.ax.set_xlim(0,1)
cb.set_clim(vmin=0,vmax=1)
cb.set_label("Fraction of sky visible [per square degree]")

plt.xticks(np.linspace(-0.5,360/xbinsize-0.5,13),np.arange(0,390,30))
plt.yticks(np.linspace(-0.5,180/ybinsize-0.5,13),np.arange(-90,105,15))

plt.xlabel("RA [deg]")
plt.ylabel("DEC [deg]")


plt.savefig(folder2+"/figures/composite.png")

###plot the difference between the two methods###
plt.figure()
plt.imshow(np.transpose((matrix-matrix2)/matrix),origin="lower",cmap="RdBu")

cb = plt.colorbar()
cb.set_label("Difference in Fraction of Sky Visible [per square degree]")

plt.xticks(np.linspace(-0.5,360/xbinsize-0.5,13),np.arange(0,390,30))
plt.yticks(np.linspace(-0.5,180/ybinsize-0.5,13),np.arange(-90,105,15))

plt.xlabel("RA [deg]")
plt.ylabel("DEC [deg]")
plt.title("(slow method - Fast method) / slow method")
plt.savefig("difference.png")


#Display plots
plt.show()
