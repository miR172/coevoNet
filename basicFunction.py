import sys,math
import time
import threading

def uniqueArray(a):
        b=[]
        for i in a:
                if b.count(i)==0:b.append(i)
        return(b)


def copyArray(a,b=[]):
	new = []
	if b==[] : b=range(0,len(a))
	for i in b: new.append(a[i])
	return(new)

def intersection(a,b,unique=True): #a, b are lists, return intersection, if none return []
        i = filter(lambda x: b.count(x)!=0,a)
        i = uniqueArray(i)
        return i

class savelog:	
	def __init__(self,b='./log/net_function_tmp'):
		self.fsock = open('%s.log'%b,'a')
		self.log=b
		self.save_out=sys.stdout
		
	def go(self): 
		sys.stdout=self.fsock
		print 'Session START: ', time.ctime(),'by YQH'

	def done(self):
		print 'Session END: ', time.ctime(), 'by YQH'
		self.fsock.close()
		sys.stdout=self.save_out


#Make class heap:
def heapSort(a):
	buildMaxHeap(a)
        for i in range(len(a),-1,0):
                size -=1
		maxHeapify(a,1)

def buildMaxHeap(a,size):
        for i in range(len(a)/2-1,-1,-1):
                maxHeapify(a,i)

def maxHeapify(a,size,i):
	l=2*(i+1)-1
	r=2*(i+1)
	if l<=size and a[l] > a[i]: largest = l
	else: largest =i
	if r<=size and a[r] > a[largest]: largest =r
	if largest != i:
		a[i],a[largest] = a[largest],a[i]
		maxHeapify(a,size,largest)


