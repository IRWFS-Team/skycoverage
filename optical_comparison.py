import numpy as np
import matplotlib.pyplot as plt
import os

filename = "comparison/pole.csv"

optical = np.loadtxt(filename,delimiter=",",unpack=True, skiprows = 2)


###Define Minimum Brightness###
M_min = 14
band = "R"
theta_0 = 1./(60.)/(3.) ###isoplanetic angle
breakdown = 15


DECs = np.arange(-90,90,breakdown)
RAs = np.arange(0,360,breakdown)

xbinsize = 1.
ybinsize = 1.
N_xbins = breakdown/xbinsize
N_ybins = breakdown/ybinsize


zoom = [44,200] # (lower DEC,lower RA)


outfolder = band+"-"+str(M_min)+"_"+str(int(theta_0*3600))+" slow"
if not os.path.exists(outfolder):
        print outfolder
        os.mkdir(outfolder)
	os.mkdir(outfolder+"/figures")
        os.mkdir(outfolder+"/datafiles")
else: 
        print "WARNING."
        print "Path already exists. May overwrite data."
        #plt.pause(10)

theta_0 = theta_0/np.sqrt(2.)

for d in [DECs[max(np.where(DECs<=zoom[0])[0])]]:
	for a in [RAs[max(np.where(RAs<=zoom[1])[0])]]:
		print d," \t ",a
		path = "Dec"+str(d)+"-"+str(d+breakdown)+"_RA"+str(a)+"-"+str(a+breakdown)+"/" 

                ###Create Grid###
                
                DEC = np.linspace(d,d+breakdown,int(ybinsize/theta_0)*N_ybins+1)#IN DEGREES
                RA = []
                filled = []
                
                y_bins = np.arange(d,d+breakdown+ybinsize,ybinsize)
                bins = []
                for i in range(1,len(y_bins)):
                        bins.append([])
                        length  = max(np.cos(y_bins[i]*np.pi/180.),np.cos((y_bins[i-1])*np.pi/180.)) #length of 1 degree of RA at DEC y_bins[i] (relative to DEC=0)
                        N = int(xbinsize/(theta_0/length))*N_xbins +N_xbins #total Number of theta_0 cells at said DEC (# of cells per bin * # of bins)
                        #for j in range(int(ybinsize//theta_0)): #for each DEC bin
                        RA.append(np.linspace(a,a+breakdown,N+1)) #add a list of all (RA) cell boundaries in that bin
                       # print i
                        for j in range(int(N_xbins)):
                                bins[i-1].append(np.zeros([int(ybinsize/theta_0),N/N_xbins],dtype=bool))
                #plt.cla()
                #fig=plt.figure()
                #ax=fig.add_subplot(111)
                #for i in range(len(DEC)-1):
                #	for j in range(len(RA[i])):
                #		ax.add_patch(patches.Rectangle((RA[i][j],DEC[i]), 360./(len(RA[i])-1), theta_0, color='red', alpha=0.4))
                #plt.xlim(0,360)
                #plt.ylim(-90,90)
                #plt.show(block=False)
                ###Load data###
                #fig=plt.figure()
                #ax=fig.add_subplot(111)


                counter=0
		for lazy in range(1):


			    
                            for j in range(len(optical[0])):
                                r = np.array([optical[1,j],optical[2,j]])
                                if r[0]<zoom[1] or r[0]>=zoom[1]+1 or  r[1]<zoom[0] or r[1]>=zoom[0]+1:
                                    continue

                                counter+=1
                                y = int((r[1]-d)//ybinsize)# which DEC bin is the star in?
                                x = int((r[0]-a)//xbinsize)# which RA bin is the star in?

                                if y== breakdown/ybinsize: #deal with edge case
                                    y-=1
                                if x== breakdown/xbinsize:#deal with edge case
                                    x-=1

                                dec2 = max(np.where(DEC<=r[1])[0]) #largest DEC value less than star's DEC
                                dec = dec2 - (y)*len(bins[y][x][:,0])#index of star's DEC in its bin (above - #of values in all previous bins )
                                ra2 = max(np.where(RA[y]<=r[0])[0])#largest RA value less than star's RA
                                ra = ra2 - (x)*len(bins[y][x][0,:])#index pf star's RA in its bin (above - #of values in all previous bins )
	
                                if not bins[y][x][dec][ra]:
                                    bins[y][x][dec][ra]=True
                                    


                plt.figure()
                print counter
                together = np.zeros([len(bins[y][x][:,0]),2*len(bins[y][x][0,:])])
                together[:,:len(bins[y][x][0,:])]=bins[y][x-1]
                together[:,len(bins[y][x][0,:]):]=bins[y][x]
                #plt.imshow(together,origin="lower",cmap=plt.get_cmap("winter"))
                plt.imshow(bins[y][x],origin="lower",cmap=plt.get_cmap("winter"))
                cb=plt.colorbar()
               
                plt.xticks(np.linspace(-0.5,len(bins[y][x][0,:])-0.5,6),np.linspace(zoom[1],zoom[1]+1,6))
                plt.yticks(np.linspace(-0.5,len(bins[y][x][:,0])-0.5,6),np.linspace(zoom[0],zoom[0]+1,6))
                #plt.xticks(np.linspace(-0.5,2*len(bins[y][x][0,:])-0.5,11),np.linspace(zoom[1],zoom[1]+2,11))
                #plt.yticks(np.linspace(-0.5,len(bins[y][x][:,0])-0.5,6),np.linspace(zoom[0],zoom[0]+1,6))
                plt.xlabel("RA [degrees]")
                plt.ylabel("DEC [degrees]")
                plt.savefig(outfolder+"/figures/zoom_in_r"+str(zoom[1])+"_d"+str(zoom[0])+".png")
                np.savetxt(outfolder+"/datafiles/zoom_in_r"+str(zoom[1])+"_d"+str(zoom[0])+".txt",np.transpose(bins[y][x]))
                #plt.savefig(outfolder+"/figures/HUDF-FP.png")
                #np.savetxt(outfolder+"/datafiles/HUDF-FP_seeing",np.transpose(together))

                       

plt.show()
plt.close("all")

