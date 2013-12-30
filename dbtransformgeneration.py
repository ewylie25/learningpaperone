'''
@author: Andrea Cadeddu, andrea.cadeddu@northwestern.edu
'''
import pymongo



class CaseListGen(object):
    '''
    general db crawler class
    '''


    def __init__(self):
        '''
        Constructor
        '''
        #connect to localhost
        self.dbr = pymongo.MongoClient().datatest
        #self.dbr = pymongo.MongoClient().data
    
    def loadReactionsIDs(self, filename):
	f=open(filename, "r")
	self.ids=[]
	for lines in f:
	    line =lines.split()
	    self.ids.append(int(line[0]))  
        
    def genRetroEntries(self):
        '''
        this generate a list which contains ALL the db entry;
        '''
        cursor = self.dbr.retro.find()
        self.retroentries = [a for a in cursor]
        cursor.close()
    
    def getRetroEntries(self):
        if not self.retroentries:
            self.genRetroEntries()
        return self.retroentries
    
    def genReactionList(self):
        cursor = self.dbr.retro.find()
        self.reactionlist=[a['rxn'] for a in cursor if a["rxn"] is not None]
        cursor.close()
        
    def getReactionList(self):
        if not self.reactionlist:
            self.genReactionlist()
        return self.reactionlist

    def rmRetroEntriesMatchingID(self):
	self.retroentries=[entry for entry in self.retroentries if entry["_id"] not in self.ids]
	
            
    def rmRetroEntriesNotMatchingID(self):
	self.retroentries=[entry for entry in self.retroentries if entry["_id"] in self.ids]


    def output(self, fname):
        '''
        '''
        f=open(fname,"w")
        for each in self.retroentries:
		eid=each["_id"]
	    	ename= each["name"]
		try:
	    	   etrans=each["reaction_smarts"]
		except:
		   etrans=""
	    	printline=str(eid)+"\t"+ename + "\t" + etrans + "\n"
	        try:
            	   f.write(printline)
		except:
		   print(""+str(eid)+"failed, probably weird chars present")

        
if __name__ == "__main__":
    # Testing
    cgl=CaseListGen()
    cgl.loadReactionsIDs("cycliz.txt")
    cgl.genRetroEntries()
    cgl.rmRetroEntriesNotMatchingID()
    cgl.output("reactions.txt")
    
    
