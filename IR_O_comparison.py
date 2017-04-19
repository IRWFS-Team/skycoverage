#Author: Eric Shore
#Purpose: To compare the visibility of of an IR and Optical WFS for several 1 degree RA by 1 degree DEC regions.
#######NOTE: THE PROGRAMS zoom_in.py and optical_comparison.py MUST BE RUN PRIOR TO THIS######
import numpy as np
import matplotlib.pyplot as plt

#Define regions of comparison. format (lower RA, lower DEC). extent is 1 degree RA by 1 degree DEC
zoom = [[60,80],[44,200],[36,250]]
#define titles for plot (unused)
titles = ["Just Off Plane","North Galactic Pole","Site of Globular Cluster"]


#empty lists for infrared and optical data
infrared = []
optical = []


#define number of comparison regions
N=2

#loop through every comparison region
for i in range(N):
	#determine folder of the files. Assmumes limiting magnitude of 15 in J-band (IR), and 14 in R-band (Optical)
	outfolder_i = "J-15_"+str(int(1./(60.)/(3.)*3600))+" slow"
	outfolder_o = "R-14_"+str(int(1./(60.)/(3.)*3600))+" slow"
	#load files into lists
	infrared.append(np.loadtxt(outfolder_i+"/datafiles/zoom_in_r"+str(zoom[i][1])+"_d"+str(zoom[i][0])+".txt",unpack=True))
	optical.append(np.loadtxt(outfolder_o+"/datafiles/zoom_in_r"+str(zoom[i][1])+"_d"+str(zoom[i][0])+".txt",unpack=True))
	
	
	#calculate total number of visible isoplanatic areas in both optical and infrared
	s11 = np.sum(infrared[i])
	s12 = np.sum(optical[i])
	#print the fraction of them. This is the improvement of the IR WFS over the Optical one
	print "improvement: ",s11,s12,float(s11)/float(s12)
	#plot the infrared in one subplot...
	plt.title(titles[i])
	plt.subplot(2,2,2*i+1)
	plt.imshow(infrared[-1],cmap="gray",origin="Lower",interpolation="nearest")
	plt.xticks(np.linspace(-0.5,len(infrared[-1][0,:])-0.5,3),np.linspace(zoom[i][1],zoom[i][1]+1,3))
        plt.yticks(np.linspace(-0.5,len(infrared[-1][:,0])-0.5,3),np.linspace(zoom[i][0],zoom[i][0]+1,3))
	plt.title("Infrared")
	plt.xlabel("RA")
	plt.ylabel("DEC")
	

	#... and the optical in the other
	plt.subplot(2,2,2*(i+1))
	plt.imshow(optical[-1],cmap="gray",origin="Lower",interpolation="nearest")
	plt.xticks(np.linspace(-0.5,len(optical[-1][0,:])-0.5,3),np.linspace(zoom[i][1],zoom[i][1]+1,3))
	plt.yticks([10000],[10000])#junk entry to remove tick labels 
       
	plt.title("Optical")
	plt.xlabel("RA")

	#save figure (currently commented out)
	#plt.savefig("comparison/"+titles[i]+".png")	
#Display plots.
plt.show()

		
