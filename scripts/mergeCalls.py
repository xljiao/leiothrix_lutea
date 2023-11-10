### ****************************************************************************
###
### mergeCalls
###
### Simon Martin shm45@cam.ac.uk
###
### ****************************** USAGE ***************************************
###
### python mergeCalls.py -f file1,file2,etc. -m <method> -I <fasta_index.fai> -o <output_file> [--outSep <space/tab/comma> --minCalls]
###
### input file format must be:
### scaffold position ind1 ind2 etc...
###
### header line must be present
###
### scaffold order must be the same as in the .fai file - this will be the case if GATK was used
###
### method (-m) can be 'union' or 'intersect' or 'all'
###
### 'intersect' will return only sites PRESENT IN BOTH FILES
###
### 'union' will return sites PRESTENT IN EITHER FILE - and the other will be padded with Ns
###
### 'all' will return ALL SITES, padding Ns where necessary
###
### --outSep can be used to define the seaprator in the output file ("space", "comma" or "tab" (default))
###
### if input or output files ens in ".gz", they will be assumed to be gzipped.
###
### ****************************************************************************

import sys
import time
import gzip

def getOptionValue(option):
  optionPos = [i for i,j in enumerate(sys.argv) if j == option][0]
  optionValue = sys.argv[optionPos + 1]
  return optionValue

if "-f" in sys.argv:
  fileNames = getOptionValue("-f").split(",")
else:
  print "\nplease specify file(s) using -f file1,file2,etc. \n"
  sys.exit()

if "-I" in sys.argv:
  faiName = getOptionValue("-I")
else:
  print "\nplease specify fasta index file name using -I <file_name> \n"
  sys.exit()

if "-o" in sys.argv:
  outName = getOptionValue("-o")
else:
  print "\nplease specify output file name using -o <file_name> \n"
  sys.exit()
  
if "-m" in sys.argv:
  if getOptionValue("-m") == "union":
    method = "union"
  elif getOptionValue("-m") == "intersect":
    method = "intersect"
  elif getOptionValue("-m") == "all":
    method = "all"
else:
  print "\nplease specify method as -m union / intersect / all \n"
  sys.exit()

if "--minCalls" in sys.argv:
  minCalls = int(getOptionValue("--minCalls"))
else:
  minCalls = None

if "--outSep" in sys.argv:
  if getOptionValue("--outSep") == "comma":
    outSep = ","
  elif getOptionValue("--outSep") == "tab":
    outSep = "\t"
  elif getOptionValue("--outSep") == "space":
    outSep = " "
  else:
    print "\nThe only options for --outSep are 'comma', 'tab' or 'space'\n"
    sys.exit()
else:
  outSep = "\t"

#create a list of file handles
files = []
for fileName in fileNames:
  if fileName[-3:] == ".gz":
    files.append(gzip.open(fileName, "r"))
  else:
    files.append(open(fileName, "rU"))

numFiles = len(files)

if outName[-3:] == ".gz":
  out = gzip.open(outName, "w")
else:
  out = open(outName, "w")

linesWritten = 0

#build dictionay of scaffolds
print "building scaffold length dictionary..."
fai = open(faiName, "rU")
scafLines = fai.readlines()
scafs = []
scafLens = {}
for line in scafLines:
  scaf,scafLen = line.split()[0:2]
  scafs.append(scaf)
  scafLens[scaf] = int(scafLen)
  
fai.close()

headers = [file.readline().split() for file in files]

#write headers
out.write(outSep.join([outSep.join(headers[0][0:2]), outSep.join([outSep.join(header[2:]) for header in headers])]) + "\n")

objectsList = [file.readline().split() for file in files]

print "Merging files..."
for scaf in scafs:
  for site in range(1,scafLens[scaf]+1):
    site = str(site)
    filesRepresented = 0
    outObjects = [scaf, site]
    for x in range(numFiles):
      if [scaf,site] == objectsList[x][0:2]:
        #if its a match, add the objects to the output and read in the next line for that file
        outObjects = outObjects + objectsList[x][2:]
        objectsList[x] = files[x].readline().split()
        filesRepresented += 1
      else:
        #if not a match, add Ns for this file, and dont read next line
        #but can't break because we need to move to next for those that do have the current line 
        outObjects = outObjects + ["A"]*(len(headers[x])-2)
    #so now we've created the output, but need to decide if we can write it
    if method == "all" or (method == "union" and filesRepresented >= 1) or (method == "intersect" and filesRepresented == numFiles):
      if not minCalls or len([call for call in outObjects[2:] if call != "A"]) >= minCalls:
        out.write(outSep.join(outObjects) + "\n")
        linesWritten += 1
        if linesWritten % 1000000 == 0:
          print linesWritten, "Lines written to output..."
    #and thats it. Move on to the next site in the genome

for f in files:
  f.close()

out.close()

