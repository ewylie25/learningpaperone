'''
@author: Andrea Cadeddu, andrea.cadeddu@northwestern.edu
'''
import pymongo
from stereofix import run_reactants



class CaseListGen(object):
    '''
    general class to deal with the db:
    read all the reaction entries in db and for each apply the transform on top of three matching molecules. It generates then a txt file with the results.
    '''


    def __init__(self, params):
        '''
        Constructor
        '''
        
        self.dbr = pymongo.MongoClient().datatest
        self.dbc = pymongo.MongoClient().data
        
        
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
        
    def getFlags(self,itemindex):
        r_check, e_check = False, False
        if self.retroentries.index(itemindex).haskey("r_check"):
            r_check=True
        if self.retroentries.index(itemindex).haskey("r_check"):
            e_check=True
        return r_check, e_check
            
    def genReactantList(self):
        '''
        returns a huge list of RXN:TEST:TEST:TEST (STRINGS!)
        '''
        totest=[]
        for each in self.getReactionList():
            rxn=rdkit.Chem.AllChem.rxnFromSmarts(each)
            totest.append(each, self.getMatchingChemicals(rxn))
        self.totest=totest

    def genResults(self):
        res=[]
        for each in self.totest:
            eachres=["","",""]
            for i in each[1:]: 
                try:
                    eachres[i]=run_reactants(each[0],each[1],false,false)
                except:
                    eachres[i]="RRfailed"
            res.append(eachres)
        self.res=res
        
    def getResults(self):
        if not self.res:
            self.genResults()
        return self.res
        
    def reactionHasMultipleReactants(self,reaction):
        if reaction.GetNumReactantTemplate() >= 2:
            return True
        else:
            return False
        
    def getMatchingChemicals(self, reaction):
        '''
        return 3 matching chemicals (in rkdit_mol format) for the reaction
        this should be run only once per reaction...
        '''
        cursor = self.dbc.chemical.find()#three need only three
        try:
            reactants=[]
            for chemical in cursor:
                mol_obj=rdkit.Chem.MolFromSmiles(chemical["smiles"])
                motif=reaction.GetReactionTemplate(0)
                if mol_obj.HasSubstructMatch(motif): # I'm not sure how to specify the motif
                    reactants.append(rdkit.MolToSmiles(mol_obj))
                    #if multiple reactants are needed for this reaction, we need all of them!
                    if self.reactionHasMultipleReactants(reaction):
                        other_reactants=[rdkit.MolToSmiles(reaction.GetReactionTemplate(i)) for i in range(reaction.GetNumReactantTemplates())]
                        rlist=[mol_obj, other_reactants[1:]]
                        reactants[-1]=".".join(rlist)
                if len(reactants)>2:
                    break
        finally:
            cursor.close()
        return reactants


            
    def output(self, fname):
        '''
        '''
        f=open(fname,"w")
        for each in self.res:
            f.write("\t".join(each)+"\n")
        

        
if __name__ == "__main__":
    # Testing
    cgl=CaseListGen()
    cgl.genRetroEntries()
    
    
    
