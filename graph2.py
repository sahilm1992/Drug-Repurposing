#OM

from math import*
import pickle
import networkx as nx

def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)



def jaccard_similarity(x,y):
	if(len(x) ==0 or len(y)==0):
		return 0
	intersection_cardinality = len(set.intersection(*[set(x), set(y)]))
	union_cardinality = len(set.union(*[set(x), set(y)]))
	return intersection_cardinality/float(union_cardinality)

drug_Features = load_obj('drug_Features')

drugKeys = drug_Features.keys()




G = nx.Graph()


drugKeys = list(drugKeys)
for drug1 in drugKeys[0:600]:
	#print (drug1)
	gene1 = drug_Features[drug1]['GENES']
	SE1 = drug_Features[drug1]['SE']

	for drug2 in drugKeys[0:2000]:
		if(drug1 is not drug2):

			gene2 = drug_Features[drug2]['GENES']
			SE2 = drug_Features[drug2]['SE']
			
			jaccard_similarity_GENE = jaccard_similarity(gene1,gene2)
			jaccard_similarity_SE = jaccard_similarity(SE1,SE2)
			#print (jaccard_similarity_GENE)
			if(jaccard_similarity_GENE >0.5):
				G.add_edge(drug1,drug2,weight=jaccard_similarity_GENE)
				print (drug1,drug2,jaccard_similarity_GENE)

print (G.number_of_nodes())

#x=1
#y=2
#G.add_edge(x,y,weight=len(ans))



print ( nx.connected_components(G))


nx.write_adjlist(G,"graphs/graph.adjlist")

nx.write_graphml(G, "graphs/cyto.gml")
nx.write_graphml(G, "graphs/test.graphml")
