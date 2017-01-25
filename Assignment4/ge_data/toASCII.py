
###
###This script can only be run in a machine with linux OS, capable of running hexdump
###
import os
def toASCII(n):
    dim = str(n)
    #Load files#
    # sample input : size64x64_bin.sol
    # sample output: size64x64_ascii.sol
    in_filename = 'size'+dim+'x'+dim+'_bin.sol' #output binary
    out_filename = 'size'+dim+'x'+dim+'_ascii.sol' #input ASCII
    command =  r"""hexdump -v -e '1/8 "%9f"' -e '"\n"' """+ in_filename +">" + out_filename
    os.system(command)

#If you want to print all matrices to ASCII at once, run the following:
#dimensions = [8, 16, 32, 64, 512, 1024, 2048, 4096, 8192]
#for d in dimensions:
#    toASCII(d)

#Otherwise run the following:
d = 64
toASCII(d)
