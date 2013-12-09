import sys, math, linecache, numpy
import fileProcess as fp
import basicFunction as bf
import blastSeq as bs

#log=bf.savelog("AnalysisTrees")
#log.go()

##STEP1 TRANSLATE PROTEIN ID INTO GENE KEYS
tranx2gene=fp.getList("data/annotation/full.annotation",col=[2,1],harshIt=True)

oufile = open("data/align/treeCor.translate",'w')
infile = "data/align/treeCor"

for i in range(3,fp.countLine(infile)+1):
	l = linecache.getline(infile,i).split("\t")
	try: 
		l[0] = tranx2gene[l[0]]
		l[1] = tranx2gene[l[1]]
	except KeyError:
		print "\tNo keys found for: ",l[0],l[1] 
		continue
	line = "\t".join(l)
	oufile.write(line)

oufile.close()

##STEP2 NORMALIZE COEFFICIENTS
infile = "data/align/treeCor.translate"
oufile = open("%s.normalize"%infile,'w')

p1 = fp.getList(infile,[0],unique=False)
p2 = fp.getList(infile,[1],unique=False)
pair = []
for i in range(len(p1)):
	pair.append([p1[i],p2[i]])

pairOrder=[i[0] for i in sorted(enumerate(pair), key=lambda x:x[1])]


while len(pairOrder)!=0:
	i=pairOrder.pop()
	group,values,p = [],[],pair[i]	#group of pair with same name, values of group, name of the group
	for j in pairOrder:
		if pair[j]==p: group.append(j)	
	tmp = map(lambda x: pairOrder.remove(x),group)
	group.append(i)
	for j in group:
		values.append(float(linecache.getline(infile,j+1).strip("\n").split("\t")[2]))
	valueM = sum(values)/len(values)
	P = "\t".join(p)
	oufile.write("%s\t%.5f\n"%(P,valueM))

oufile.close()

infile = "%s.normalize"%infile
#infile = "data/align/treeCor.translate"
ppifile = "data/annotation/Pilot_PPI_dataset_unfiltered.txt"
p1 = fp.getList(ppifile,[0],unique=False)
p2 = fp.getList(ppifile,[2],unique=False)
p2 = map(lambda x: x.strip("\r"),p2)

oufile = open("./data/relation/mirrorTreeCoef",'w')

bg={}
for i in range(len(p1)):
	if p1[i] ==p2[i]:continue
        for j in range(1,fp.countLine(infile)+1):
                l = linecache.getline(infile,j).strip("\n").split("\t")
		if l.count(p1[i])!=0 and l.count(p2[i])!=0:
			line = "\t".join(l)
			print p1[i],p2[i],line
			oufile.write("%s\n"%line)
			continue

oufile.close()

#infile = "./data/relation/mirrorTreeCoef"
#for i in range(1,fp.countLine(infile)+1):
#	l = linecache.getline(infile,j).strip("\n").split("\t")
	

#log.done()


