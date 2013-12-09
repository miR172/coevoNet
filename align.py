import blastSeq as bs
import math,time
import basicFunction as bf
###
#test code, fasta_aln file: ./data/alignTest/test.fasta_aln 
#a = bs.alignment("test",d="./data/alignTest/")
#a.getAuxMatrix()


genelist = open("data/align/genelist",'r').readlines()
genelist = map(lambda x: x.strip("\n"),genelist)

alnlist = [] #algnment list in the order of genelist
for g in genelist:
	try: aln = bs.alignment(g,d="./data/align/trees/")
	except IOError: continue
#	aln.printDistMatrix(d="./data/align/distanceMatrix/") 
	alnlist.append(aln)


outf = open("./data/align/treeCor",'w')
outf.write(time.ctime())
outf.write("\nNonTolCorrectionRelation\n") #####change notes here!

for i in range(0,len(alnlist)-1): 
	A = alnlist[i].getAuxMatrix() #capital for Matrix
	s1 = A[0]
	v1 = math.sqrt(alnlist[i].varSum*10000)	#protect small number
	g1 = alnlist[i].gene
	print g1
	for j in range(i+1,len(alnlist)):
		B = alnlist[j].getAuxMatrix()
		s2 = B[0]
		v2 = math.sqrt(alnlist[j].varSum*10000) #protect small number
		g2 = alnlist[j].gene
		R = 0
		myS = bf.intersection(s1,s2)
		if len(myS) <= 1:continue #not enough species set	
		print "\t", g2
		As = map(lambda x: s1.index(x),myS)
		Bs = map(lambda x: s2.index(x),myS)
		Ai = [myS]
		Bi = [myS]
		for s in myS:
			Ai.append(map(lambda x: A[A[0].index(s)+1][x],As))
			Bi.append(map(lambda x: B[B[0].index(s)+1][x],Bs))		
		for n in range(1,len(Ai)):
			for m in range(n,len(Ai[0])):
				R += Ai[n][m]*Bi[n][m]
		try: R = R*10000/(v1*v2)
		except ZeroDivisionError: 
			print "\tsmall variance X100:",g1,v1,g2,v2
		outf.write("%s\t%s\t%.5f\n"%(g1,g2,R))

outf.close()
		
		
	
	
