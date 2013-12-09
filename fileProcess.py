import sys,math
import linecache

def countLine(inFile):
#return line number in a file
        nlines = sum(1 for line in open(inFile,'r'))
        return nlines


def extractRow(inFile, header=True, matchList = [''], 
search='', lineInitial='', start = 0, until='', step=1):
#return col to be extracted in list, from 0 to m
#matchList can be a list
        result = [] 
	if until == '': until = countLine(inFile)
	rows = range(start, until)      #range to look at, from 1
	if header: result.append(0)
        for n in rows:
                line = linecache.getline(inFile,n+1)
                if not line.startswith(lineInitial): continue
                if line.find(search)==-1:continue
                for item in matchList:
                        if line.find(item)!=-1 and (n-start)%step ==0: result.append(n)
        return result

##	for line in open(inFile,'r'): #the line now is i+1, last line is i
## 		i+=1
##		if header and i==0: continue
##		if i<start or i>until: continue
##		if not line.startswith(lineInitial): continue
##		if line.find(search)==-1: continue
##		for item in matchList: 
##				if line.find(item)!=-1 and (i-start)%step==0 : row.append(i)
##	return row


def extractCol(inFile,bufferLine=0,split="\t"):
#return col to be extracted, from 0 to n
	i = -1
        for line in open(inFile,'r'):
		i+=1
		if i==0: headers = line.strip('\n').split('\t')
		if i > bufferLine : continue
		print line
	for h,j in enumerate(headers): print h,"\t", j
	print "Total Columns: ", len(headers),'\n'
	col = input("Column Wanted (as variable in py): ")
	print "You Picked: ", col, "\n"
	try: col =  list(col)
	except TypeError: col=[col]
	return(col)


def getList(inFile,col,row=[],harshIt=False,mutual=False, split="\t", unique=True):
	mylist = [];
	myharsh={}
	i = -1
	for line in open(inFile,'r'):
		i+=1
		if row!=[] and row.count(i)==0: continue
		chopLine = line.strip('\n').split('\t')
		if harshIt:
			myharsh[chopLine[col[0]]] = chopLine[col[1]]
			if mutual: myharsh[chopLine(col[1])] = chopLine(col[0])
			continue
                if unique:
                        if mylist.count(chopLine[col[0]])!=0 : continue
                mylist.append(chopLine[col[0]])
	if harshIt: 
		return myharsh
	return mylist


def extractTable(inFile, outFile, row=[], col=[]):
	if row==[]: 
		nline = countLine(inFile)
		row = range(0,nline)
	outf = open(outFile, 'w')
	i=-1
	for line in open(inFile,'r'):
		i+=1
		if row.count(i)==0: continue
		tmp = row.pop(row.index(i))	#reduce list to make program faster
		if col==[]: 
			outf.write(line)
			continue
		chopLine=line.strip("\n").split("\t")
		for i in col:	#some day I will make a class named table
			outf.write("%s\t"%chopLine[i])
		outf.write('\n')
	outf.close()
	return outFile


