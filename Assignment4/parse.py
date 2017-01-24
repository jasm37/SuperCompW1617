
#scipt to parse output file into csv

#args = getopt.getopt(argv,"hi:o:",["ifile="])

inputfile = open('io_results.out');
outputfile = open('io_all.csv','w');

problemsizes = [64,512,1024,2048,4096,8192];
procsizes = [8,16,32,64];
i=0;
j=0;
newblock = 0;
rankcount = 0;


for line in inputfile:
    if line.startswith("[R"):
        if newblock == 1:
            j=j+1;
            if j>3:
                j=0;
                i=i+1;
        #if line[0:5] == "[R00]":
        if rankcount > procsizes[j]-2:
            newblock = 1
            rankcount = 0;
        else:
            newblock = 0
            rankcount=rankcount+1;
        problems=problemsizes[i]
        proc=procsizes[j]
        line=line.replace(";","")
        line=line.replace("[R","")
        line=line.replace("]","")
        times=line.split(' ')
        outputfile.write('{}, {}, {}, {}, {}, {}, {}, {}\n'.format(problems,proc,times[0],times[3],times[5],times[7],times[9],times[11]))
#outputfile.close()
