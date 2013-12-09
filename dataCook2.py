import sys, math, linecache
import fileProcess as fp
import basicFunction as bf
import blastSeq as bs

log=bf.savelog("makeMultipleAlign")
log.go()

blastResultList = open("data/blast/list.txt").readlines()

forwardlist = []        #forward direction of blast
backwardlist = []       #backward direction of blast

for i in range(0,len(blastResultList)):
        fileName = blastResultList[i].strip("\n").strip(" ")
        br = bs.blastResult(fileName)
        print "constructed blastResult object %s"%br.fileName ##check
        br.setQueryList()
        print "setted queryList"        ##check 
        br.setTargetList()
        print "setted targetList"       ##check
        if br.querySpecies == "CE.sample.fa":
                forwardlist.append(br)
        else:
                backwardlist.append(br)


##check blastresult construction##
for i in backwardlist:
      print i.querySpecies, "\t", i.targetSpecies, "\t", len(i.getQueryList()),"\t", len(i.getTargetList()), i.fileName

for i in forwardlist:
      print i.querySpecies, "\t", i.targetSpecies, "\t", len(i.getQueryList()),"\t", len(i.getTargetList()), i.fileName


for query in forwardlist:
	print query.querySpecies,":"
        queryFasta = bs.genomeSequenceFasta(query.querySpecies)
        for target in backwardlist:
                targetFasta = bs.genomeSequenceFasta(target.querySpecies)
                if target.targetSpecies==query.querySpecies and target.querySpecies==query.targetSpecies:
                        recip = bs.reciprocalBlast(query,target)
                        for q in query.getQueryList():	
                                recip.findReciprocal(q)
                                if len(recip.getHits())==0: continue
                                recip.getMultipleAlignFasta(queryFasta,targetFasta)
                        #done each target for this query
                #if not recip blast result next


log.done()
