# -*- coding: utf-8 -*-
#this is the procedure that prepares and manage data for blast and sequence alignment

import sys, math, linecache
import fileProcess as fp
import basicFunction as bf
import numpy

class genomeSequenceFasta:	#deal with fasta files
	
	def __init__(self,speciesFile,d="./data/sequence/"):
		self.species = speciesFile
		self.fileName = "".join([d,speciesFile])
		self.queryTitles = []
		self.queryEnds = []
		self.queryList = []
		self.seqBlocks = []

        def reset(self):
                self.queryList=[]
                self.queryTitles=[]
                self.queryEnds=[]
                self.seqBlocks=[]
	
	def setQuery(self,queryList):   #queryList is array!
                self.queryList = []
		self.queryList = queryList
		self.queryTitles = []
		self.setQueryTitles()
                self.setQueryEnds()
                self.setSeqBlocks()

	def setQueryTitles(self):
		if self.queryList==[]:
			print "You need to set query list first. Use genomeSequenceFasta.setQuery(list)"
			return
		self.queryTitles = fp.extractRow(self.fileName,header=False,lineInitial='>', matchList=self.queryList)

	def setQueryEnds(self,appxBlockSize=100):
                self.queryEnds = []
                maxline = fp.countLine(self.fileName)-1
		if self.queryTitles==[]:
			print "query not in file %s: "%self.fileName, self.queryList[0]," ..."
			return
		for title in self.queryTitles:
                        start = title+2
                        end = title
                        while start < maxline: 
                                tmp = linecache.getline(self.fileName,start)
                                start += 1
                                if tmp.startswith(">"):
                                        end = start-2
                                        break
			if start >=maxline: end = maxline
			self.queryEnds.append(end)

	def setSeqBlocks(self):
                self.queryBlocks = []
		if self.queryTitles==[] or self.queryEnds==[] :
			print "No match of queryList: ",self.queryList[0]," ..."
			return
		for i in range(0,len(self.queryTitles)):
        		block = range(self.queryTitles[i],self.queryEnds[i])
			self.seqBlocks.extend(block)
		self.seqBlocks = bf.uniqueArray(self.seqBlocks)
		self.seqBlocks.sort()        

        def getQueryList(self):
                return self.queryList

       	def getQueryTitles(self):
                return self.queryTitles
      
	def getQueryEnds(self):
                return self.queryEnds
            
	def getSeqBlockFasta(self,fileName,append=True):
                if append: outf = open(fileName,'a')
                else: outf = open(fileName,'w')
		for b in self.seqBlocks:
        		outf.write(linecache.getline(self.fileName, b+1))
		outf.close()


	
class blastResult:

	def __init__(self,resultFile,d="./data/blast/"):
		self.fileName = "".join([d,resultFile.strip("\n").strip(" ")]) #make right directory
		self.targets = []
		self.allQueries = []
		self.allTargets = []
		try: self.querySpecies = "%s.fa"%resultFile.split(".fa.")[0]
                except IndexError:
                    self.querySpecies = ""
                    print "unkown querySpecies, for file: %s"%self.fileName
		try: self.targetSpecies = "%s.fa"%resultFile.split(".fa.")[1]
                except IndexError:
                    self.targetSpecies = ""
                    print "unkown targetSpecies, for file: %s"%self.fileName
                try: self.cutoff = self.fileName.split(".fa.")[2]
                except IndexError:
                    self.cutoff = ""
                    print "unkown cutoff"
#               self.setQueryList()
#               self.setTargetList()
                
	def setQueryList(self):
                if self.allQueries != []:
                    d = input("Take Long: Seriously Reload All QueryList? [y/n] :\nFor file: %s"%self.resultFile)
                    if d!="y": return
		self.allQueries = fp.getList(self.fileName,[0]) #unique in this method default true

	def setTargetList(self):
		self.allTargets = fp.getList(self.fileName,[1]) #unique in this method default true
	
        def setQuerySpecies(self,qspecies):
                self.querySpecies = qspecies

        def setTargetSpecies(self,tspecies):
                self.targetSpecies = tspecies

        def setCutoff(self,cutoff):
                self.cutoff = cutoff
                
	def getQueryList(self):
		return(self.allQueries)
	
	def getTargetList(self):
		return(self.allTargets)

        def getFileName(self):
                return(self.fileName)

            
class hit:

	def __init__(self,query="", aBlastResult=""):
		self.query = query
		self.mom = aBlastResult #aBlastResult is a pointer to class blastResult
		self.hits = []
		self.targets = []
                self.best = []  #list of best hits
                self.identityCutoff = 0
                self.alignLengthCutoff = 0
                self.evalueCutoff = 0
                self.identities=[]
                self.alignLengths=[]
                self.evalues=[]
                             
	def setMom(self,result):
		self.mom = result
                if self.mom.getQueryList()==[]: self.mom.setQueryList()

	def setQuery(self,query):
                self.query=""
                self.hits=[]
                self.targets=[]
                self.identities=[]
                self.alignLengths=[]
                self.evlues=[]
		if self.mom=="":
			print "You need to set a blastResults to look at, Use hit.setMom=blastResults(resultfile)"
			return
		if self.mom.getQueryList().count(query)==0:
			print "Query: %s has no hits in Species: %s above cutoff"%(query, self.mom.targetSpecies)
			self.targets = []
                        return
		self.query = query
		self.setHits(self.query)
		
        def setHits(self,query=""):
		if self.query=="":
			print "You need to set a query. Use blastResults.setQuery(query)"
			return
		self.hits = fp.extractRow(self.mom.fileName, header=False, lineInitial="%s\t"%self.query)
		self.targets = bf.uniqueArray(fp.getList(self.mom.fileName,row=self.hits,col=[1]))      #if improve, these 1,2,3...numbers of cols should be masked in class blastResult
                self.identities = bf.uniqueArray(fp.getList(self.mom.fileName,row=self.hits,col=[2]))
                self.alignLengths = bf.uniqueArray(fp.getList(self.mom.fileName,row=self.hits,col=[3]))
                self.evalues = bf.uniqueArray(fp.getList(self.mom.fileName,row=self.hits,col=[11]))

##        def setIdentityCut(self,cutoff=0):
##                self.identityCutoff = cutoff
##
##        def setLengthCut(self,cutoff=0):
##                self.alignLengthCutoff = cutoff
##
##        def setEvalueCut(self,cutoff=0):
##                self.evalueCutoff = cutoff
##   
##        def setBest(self,inclusive=True):
##                if (inclusive == False):
##                        self.best = filter(lambda x: identity[x] > self.identityCutoff
##                                and alignLength[x] > self.alignLengthCutoff
##                                and evalue[x] < self.evalueCutoff, range(0,len(self.targets)) )
##                else:   self.best = filter(lambda x: identity[x] >= self.identityCutoff
##                           and alignLength[x] >= self.alignLengthCutoff
##                           and evalue[x] <= self.evalueCutoff, range(0,len(self.targets)) ) 
##                

        def getMom(self):
                return(self.mom)
            
        def getQuery(self):
                return(self.query)
            
	def getTargetList(self):
                return(self.targets) #self.targets are list of target gene names(col 2 in blast result files)

	def getHits(self):
		return(self.hits)

        def getBest(self):
                return(self.best)

        def getBestTargets(self):
                if self.best==[]: return []
                result = bf.copyArray(self.targets,self.best)
                return self.bestHits
                          
class reciprocalBlast:
        def __init__(self, query, target):      #take two blastResult object as input
                self.queryBlast = query
                self.targetBlast = target
                self.querySpecies = self.queryBlast.querySpecies
                self.targetSpecies = self.targetBlast.querySpecies
                self.queryHit = hit(aBlastResult=self.queryBlast)
                self.targetHit = hit(aBlastResult=self.targetBlast)
                self.finalHits = []
                self.currentQuery = ""

        def setCutoff(self):
        #to be update if necessary
        #need to call cutoff setters in class hit according to option?!
                return
        
        def findReciprocal(self,queryName,top=1): #queryname is the query for QueryBlast(forward!)
                self.queryHit.setQuery(queryName)
                print "checking: ",queryName
                origin = self.queryHit.getTargetList() #target in target, query in query
                final = []
                for t in origin:
                        self.targetHit.setQuery(t)
                        reflect = self.targetHit.getTargetList() #target in query, query in target
                        if reflect.count(queryName) != 0 :
                                print "\t got ", t
                                final.append(t)
                                if len(final) ==top:
                                        break #only look at top seq
                final = bf.uniqueArray(final)
                self.finalHits = final
                self.currentQuery = queryName

        def getHits(self):
                return self.finalHits

        def getMultipleAlignFasta(self,faQ,faT):        #take two genomeSequenceFasta object as input
                if self.finalHits ==[]:
                        print "No reciprocal hits to align."
                        return
                faQ.reset();faT.reset()
                faQ.setQuery([self.currentQuery])
                faQ.getSeqBlockFasta("./data/align/%s.self"%self.currentQuery, False)
                faT.setQuery(self.finalHits)
                faT.getSeqBlockFasta("./data/align/%s.targets"%self.currentQuery,True)
                
class alignment:
        def __init__(self,gene,d="./data/align/"):
                self.gene = gene
                self.fileName = "%s%s.fasta_aln"%(d,self.gene)
                self.starts = fp.extractRow(self.fileName,header=False,lineInitial=">")
                self.titles = []
                self.species=[]
                self.ends= []
                self.seq=[]
                self.auxMatrix = []
                self.setSeq()
                self.distMatrix =[[0 for x in xrange(len(self.titles))] for x in xrange(len(self.titles)+1)]
                self.distMatrix[0] = self.titles
                self.mean,self.varSum,self.var = 0,0,0
                self.setDistMatrix()

        def setSeq(self): #called in __init__
                for i in range(0,len(self.starts)):
                        t = linecache.getline(self.fileName,self.starts[i]+1).strip(">").split("  ")[0]
                        self.titles.append(t)
                self.species=bf.copyArray(self.titles)
                for i in range(0,len(self.starts)-1):
                        self.ends.append(self.starts[i+1]-1)
                self.ends.append(fp.countLine(self.fileName)-1)
                for i in range(0,len(self.titles)):
                        self.seq.append("")
                        for j in range(self.starts[i]+1,self.ends[i]+1):
                                self.seq[i]="".join([self.seq[i],linecache.getline(self.fileName,j+1).strip("\n")])
                for i in range(0,len(self.species)):
                        if self.species[i].startswith("CBG"):
                                self.species[i]="CBRI";continue
                        if self.species[i].startswith("CBN"):
                                self.species[i]="CBRE";continue
                        if self.species[i].startswith("CJA"):
                                self.species[i]="CJ"; continue
                        if self.species[i].startswith("CRE"):
                                self.species[i]="CR"; continue
                        self.species[i]="CE"
                #self.species in the same order with self.titles, self.starts, self.ends, self.sequence
                self.auxMatrix.append([i[0] for i in sorted(enumerate(self.species), key=lambda x:x[1])])
                #things in () will return index of origin array with array untouched!        

        def computeDist(self,s1,s2): #called in setDistMatrix()
                n = len(s1)
                m = 0
                for i in range(0,len(s1)): 
                        if s1[i]==s2[i]:
                                if s1[i]=="-":
                                        n -= 1
                                        continue
                                m += 1
                try: return 1-1.0*m/n #"严钦骅你是不是有病，明明直接算mismatch-if语句!=就好了非要算match然后用1-过去……2货！"
                except ZeroDivisionError:
                        print "\tEmptyAlignment:",s1[::10],"...\t", s2[::10],"...\n"
                        return 1
                
        def setDistMatrix(self): #called in __init__
                total,n = 0,0
                for i in range(0,len(self.titles)):
                        for j in range(0,len(self.titles)):
                                if i==j:
                                        self.distMatrix[i+1][j] = 1
                                else:
                                        if self.distMatrix[j+1][i]!=0:
                                                self.distMatrix[i+1][j] = self.distMatrix[j+1][i]
                                        else:
                                                self.distMatrix[i+1][j] = self.computeDist(self.seq[i],self.seq[j])
                                                total+=self.distMatrix[i+1][j]
                                                n += 1
                try: self.mean = total/n
                except ZeroDivisionError:
                        print "0 records in ",self.gene, self.distMatrix
                        self.distMatrix=[0]
                        self.auxMatrix=[0]
                        return
                for i in self.auxMatrix[0]:
                        tmp = map(lambda x: self.distMatrix[i+1][x]-self.mean,self.auxMatrix[0])
                        self.auxMatrix.append(tmp)
                self.auxMatrix[0]=map(lambda x: self.species[x],self.auxMatrix[0])
                for r in range(1,len(self.distMatrix)):
                        for c in range(r,len(self.distMatrix[r])):
                                self.varSum += (self.distMatrix[r][c]-self.mean)**2
                self.var = self.varSum/n
        
        def getAuxMatrix(self): 
                return self.auxMatrix

        def printDistMatrix(self,d="./data/distanceMatrix/"):
                outf = open("%s%s.dm"%(d,self.gene),'w')
                for i in self.distMatrix:
                        for j in range(0,len(i)):
                                outf.write("%s\t"%i[j])
                        outf.write("\n")
                outf.close()

            
