#Author: Eric Shore
#Purpose: To sort the 2MASS catalogue into different files based on position on the sky, and magnitude (for "fast" method outlined in "Method Writeup.txt")
import numpy as np
from math import ceil
import os

###Define parameters###
breakdown = 15 #size of individual patches of sky in degrees RA by degrees DEC

#create arrays of RA and DEC values at the edge of each patch
DECs = np.arange(-90,90,breakdown,dtype=float)
RAs = np.arange(0,360,breakdown,dtype=float) 

#loop through each patch
for d in DECs:
	for r in RAs:
		#determine file path of data
		path = "Dec"+str(int(d))+"-"+str(int(d+breakdown))+"_RA"+str(int(r))+"-"+str(int(r+breakdown))+"/"
		print path
		if not os.path.exists(path):#if folder does not exist create one
			os.mkdir(path)
		os.chdir(path)
		

                ### create all the files ###
                #fout = open("K_>17.txt","w")
                #fout.close()
                #fout = open("H_>17.txt","w")
                #fout.close()
                fout = open("J_>17.txt","w")
                fout.close()
                for i in range(20):
                    #fout = open("K_"+str(17-i)+"-"+str(16-i)+".txt",mode="w")
                    #fout.close()
                    #fout = open("H_"+str(17-i)+"-"+str(16-i)+".txt",mode="w")
                    #fout.close()
                    fout = open("J_"+str(17-i)+"-"+str(16-i)+".txt",mode="w")
                    fout.close()
                #fout = open("K_<-3.txt",mode="w")
                #fout.close()
                #fout = open("H_<-3.txt",mode="w")
                #fout.close()
                fout = open("J_<-3.txt",mode="w")
                fout.close()
                os.chdir("..")

#loop through every 2MASS file
for i in range(92):
        #determine the 3 letter code at the beginning of the filename
	if i<26:
		string = "aa"+chr(97+i)
	elif i<52:
		string = "ab"+chr(97+i-26)
	elif i<57:
		string = "ac"+chr(97+i-52)
	elif i<83:
		string = "ba"+chr(97+i-57)
	else:
		string = "bb"+chr(97+i-83)
	print string
        #open file
	cat_part = open("../psc_"+string,mode="r") 
        #loop through file
	while(True): #attempt to read from file, if this fails, then we are at the end of the file
    		try:
                        #read a line
			text = cat_part.readline()
                        #recover RA and DEC of the star
			data = text.split("|")
			RA = data[0]
			DEC = data[1]
                        #determine the patch of sky the star is in
                        ra = RAs[max(np.where(RAs<=float(RA))[0])]
                        dec = DECs[max(np.where(DECs<=float(DEC))[0])]

                        #determine the folder where the star should be stored
                        folder = "Dec"+str(int(dec))+"-"+str(int(dec+breakdown))+"_RA"+str(int(ra))+"-"+str(int(ra+breakdown))+"/"
                        #determine the brightness of the star
#			if data[14]=='\N':
#				K_mag = np.nan
#			else:
#				K_mag = float(data[14])
#			if data[10]=='\N':
#				H_mag = np.nan
#			else:
#				H_mag = float(data[10])
			if data[6]=='\N':
				J_mag = np.nan
			else:
				J_mag = float(data[6])
                        
                        #determine what to ouput into the sorted file
			output = "%f %f %f \n" %(float(RA),float(DEC),J_mag)
			###########################################
			###########  H-band sort  #################
			###########################################
#			if np.isnan(H_mag):
#				pass
#			elif H_mag>17:
#				fout = open(folder+"H_>17.txt","a")
#				fout.write(output)
#				fout.close()
			
#			elif H_mag<-3:
#				fout = open(folder+"H_<-3.txt",mode="a")
#				fout.write(output)
#				fout.close()
				
#			else:
#				H = int(ceil(H_mag))
#				fout = open(folder+"H_"+str(H)+"-"+str(H-1)+".txt",mode="a")
#				fout.write(output)
#   				fout.close()
		
			###########################################
			###########  J-band sort  #################cieling in python
			###########################################
			if np.isnan(J_mag):
				pass
			elif J_mag>17:
				fout = open(folder+"J_>17.txt","a")
				fout.write(output)
				fout.close()
				
			elif J_mag<-3:
				fout = open(folder+"J_<-3.txt",mode="a")
				fout.write(output)
				fout.close()
				
			else:
				J = int(ceil(J_mag))
				fout = open(folder+"J_"+str(J)+"-"+str(J-1)+".txt",mode="a")
				fout.write(output)
   				fout.close()
		
			###########################################
			###########  K-band sort  #################
			###########################################
#			if np.isnan(K_mag):
#				pass
#			elif K_mag>17:
#				fout = open("K_>17.txt","a")              
#				fout.write(text)
#				fout.close()
#				
#			elif K_mag<-3:
#				fout = open("K_<-3.txt",mode="a")
#				fout.write(text)
#				fout.close()
#			
#			else:
#				K = int(ceil(K_mag))
#				fout = open("K_"+str(K)+"-"+str(K-1)+".txt",mode="a")
#				fout.write(text)
 #   				fout.close()
		except(IndexError):
			if data==[""]:
				print "End of file"
				break
			else: 
				print "Error"
				continue
	cat_part.close()
	
	
		
