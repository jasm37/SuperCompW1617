#Runs with python 2.7 
#This script should be run in the folder where all the matrices files are.
#As output it will create binary versions of the files .vec and .mat
from array import array
def toBinary(n, fileformat):
    dim = str(n)
    #Load files
    in_filename = 'size'+dim+'x'+dim+'.'+fileformat #input ASCII
    out_filename = 'size'+dim+'x'+dim+'_bin.'+fileformat #output binary
    in_file = open(in_filename,'r')
    out_file = open(out_filename, 'wb')
    l = [float(i) for i in in_file.read().split()]
    #print l
    double_l = array('d',l)
    #to binary
    double_l.tofile(out_file)
    
    ##Close files
    in_file.close()
    out_file.close()


dimensions = [8, 16, 32, 64, 512, 1024, 2048, 4096, 8192]
formats = ['vec','mat']
#### there might memory issues with big dimensions(bigger than 8192)!
for d in dimensions:
    for fmt in formats:
        toBinary(d,fmt)
#toBinary(8192,'mat')
