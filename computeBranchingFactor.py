import pymongo
import numpy as np
db=pymongo.MongoClient().new_data
#cursor=db.chemical.find({"cost":{"$exists":True}})
cursor=db.chemical.find({"cost":{"$ne":None}})
totinconnections=0
totoutconnections=0
totentries=0
entries=np.array([[entry['in_degree'],entry['out_degree']] for entry in cursor])
#print entries[:,0]
#for each in cursor:
#    totalinconnections=totalinconnections+each['in_degree']
#    totaloutconnections=totaloutconnections+each['out_degree']
#    totentries=totentries+1
inconnections=entries[:,0]
outconnections=entries[:,1] #np.array([entry['in_degree'] for entry in entries])
#outconnections=np.array([entry['out_degree'] for entry in entries])
print "tot in connections: ", np.sum(inconnections), "over ",len(inconnections),"compounds, avg: ", np.average(inconnections), " std: ", np.std(inconnections)
print "tot out connections: ", np.sum(outconnections), "over ",len(inconnections),"compounds, avg: ", np.average(outconnections), " std: ", np.std(outconnections)
#print "tot_in: ", totalinconnections, "tot_out: ", totaloutconnections, "tot_entries: ", totentries
