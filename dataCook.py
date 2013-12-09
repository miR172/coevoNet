#this script lead you all the way to done with reciprocal blast

import sys, math, linecache
import fileProcess as fp
import basicFunction as bf
import blastSeq as bs

jimmyPPIFile = "data/annotation/Pilot_PPI_dataset_unfiltered.txt"
annotationFile = "data/annotation/full.annotation"

genelist = fp.getList(jimmyPPIFile,[0])
genelist.extend(fp.getList(jimmyPPIFile,[2]))
genelist = bf.uniqueArray(genelist)
genelist = map(lambda x: x.strip("\r"),genelist)

if genelist.count("Bait") !=0:genelist.remove("Bait")
if genelist.count("Protein_symbol") !=0:genelist.remove("Protein_symbol")

annotationlist=fp.extractRow(annotationFile,header=True,matchList=genelist) 
len(annotationlist)
sequenceIDList = []
for i in annotationlist:
	 sequenceIDList.append(linecache.getline(annotationFile,i+1).split("\t")[2])

sequenceIDList = bf.uniqueArray(sequenceIDList)

#got sequence ID list

sequenceFile = "Caenorhabditis_elegans.WBcel235.73.pep.all.fa"
sampleFasta = "data/sequence/CE.sample.fa"

genome = bs.genomeSequenceFasta(sequenceFile)
genome.setQuery(sequenceIDList)
genome.getSeqBlockFasta(sampleFasta)


#### GO TO HPC RUN BLAST ####
#### AND CONTINUE AFTER  ####


blastResultList = open("data/blast/list.txt").readlines()

gfList =[]
for i in blastResultList:
	i=i.strip("\n").strip(" ")
	f="".join(["data/blast/",i])
	query = fp.getList(f,[1])
	query = bf.uniqueArray(query)
	gf = bs.genomeSequenceFasta("%s.fa"%f.split(".fa.")[1])
	gf.setQuery(query)
	gf.getSeqBlockFasta("data/recip/%s.sample.fa"%gf.species) #prepare reciprocal fasta files from target species genome
	gfList.append(gf)


#### GO TO HPC RUN BLAST|CHANGE NAME TO BE CLEAR IF NEED ####
#### AND CONTINUE AFTER GET RECIPROCAL ####


forwardlist = []	#forward direction of blast
backwardlist = []	#backward direction of blast

for i in range(0,len(blastResultList)):
	fileName = blastResultList[i].strip("\n").strip(" ") 
	br = bs.blastResult(fileName)
	print "constructed blastResult object %s"%br.fileName ##check
	br.setQueryList()
	print "setted queryList"	##check 
	br.setTargetList()
	print "setted targetList"	##check
	if br.querySpecies == "CE.sample.fa":
		forwardlist.append(br)
	else:
		backwardlist.append(br)


##check blastresult construction##
for i in backwardlist:
      print i.querySpecies, "\t", i.targetSpecies, "\t", len(i.getQueryList()),"\t", len(i.getTargetList()), i.fileName

for i in forwardlist:
      print i.querySpecies, "\t", i.targetSpecies, "\t", len(i.getQueryList()),"\t", len(i.getTargetList()), i.fileName

##	query		target		queries	targets
##	CBRE.sample.fa 	CEsample.fa 	11488 	16717
##	CBRI.sample.fa 	CEsample.fa 	8463 	16744
##	CJ.sample.fa 	CEsample.fa 	10275 	15857
##	CR.sample.fa 	CEsample.fa 	12123 	16973

##	CEsample.fa 	CBRE.sample.fa 	4479 	11499
##	CEsample.fa 	CJ.sample.fa 	4366 	10298
##	CEsample.fa 	CBRI.sample.fa 	4511 	8482
##	CEsample.fa 	CR.sample.fa 	4546 	12137


for query in forwardlist:
	queryFasta = bs.genomeSequenceFasta(query.querySpecies)
	for target in backwardlist:
		targetFasta = bs.genomeSequenceFasta(target.querySpecies)
		if target.targetSpecies==query.querySpecies and target.querySpecies==query.targetSpecies:
			recip = bs.reciprocalBlast(query,target)
			for q in query.getQueryList():
				hits = recip.findReciprocal(q,top=1)
				if hits==[]:continue
				recip.getMultipleAlignFasta(queryFasta,targetFasta)
			#done each target for this query
		#if not recip blast result next
	#done each target

#use qHit.setEvalueCut(value)... and .setBest() and getBestTargets() to access top hits, for this time, no need to do 
					


