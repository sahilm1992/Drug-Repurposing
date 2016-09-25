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



kegg_dict_drugs_features=load_obj('kegg_dict_drug_features+SE')
kegg_dict_disease_features=load_obj('kegg_dict_disease_features+SE')




G = nx.Graph()


drugKeys = list(kegg_dict_drugs_features)
diseaseKeys = list(kegg_dict_disease_features)
for drug1 in drugKeys:
	#print (drug1)
	gene1 = kegg_dict_drugs_features[drug1]['GENE']
	SE1 = kegg_dict_drugs_features[drug1]['SE']
	PATHWAYS1=kegg_dict_drugs_features[drug1]['PATHWAY']
	
	for drug2 in drugKeys:
		if(drug1 is not drug2):
#			if(len(gene1) !=1):
#				continue
		#	print (len (gene1))
			
			gene2 = kegg_dict_drugs_features[drug2]['GENE']
			SE2 = kegg_dict_drugs_features[drug2]['SE']
			#if(len(SE2)>1):
			#	print (SE2)
			PATHWAYS2=kegg_dict_drugs_features[drug2]['PATHWAY']
			jaccard_similarity_GENE = jaccard_similarity(gene1,gene2)
			jaccard_similarity_SE = jaccard_similarity(SE1,SE2)
			jaccard_similarity_PATHWAYS	 = jaccard_similarity(PATHWAYS1,PATHWAYS2)
		
			wt1=0
			wt2=0
			wt3=0
			
			wt1=1 if jaccard_similarity_GENE >0.5 else 0

			wt2=1 if jaccard_similarity_SE >0.5 else 0

			wt3=1 if jaccard_similarity_PATHWAYS >0.5 else 0
		#	if(wt2==1):
		#		print ("ues")
	
			edgeWt = wt1+wt2+wt3
			#print (jaccard_similarity_GENE)
			if(edgeWt >1):
				G.add_edge(drug1,drug2,weight=edgeWt)
				#print ('+',kegg_dict_drugs_features[drug1] ,'-',kegg_dict_drugs_features[drug2])
				print ('+',drug1,'-',drug2,edgeWt)




for disease1 in diseaseKeys:
	#print (drug1)
	gene1 = kegg_dict_disease_features[disease1]['GENE']
	#SE1 = kegg_dict_disease_features[drug1]['SE']
	PATHWAYS1=kegg_dict_disease_features[disease1]['PATHWAY']

	for disease2 in diseaseKeys:
		if(disease1 is not disease2):

			gene2 = kegg_dict_disease_features[disease2]['GENE']
#			SE2 = kegg_dict_disease_features[disease2]['SE']
			PATHWAYS2=kegg_dict_disease_features[disease2]['PATHWAY']
			jaccard_similarity_GENE = jaccard_similarity(gene1,gene2)
#			jaccard_similarity_SE = jaccard_similarity(SE1,SE2)
			jaccard_similarity_PATHWAYS	 = jaccard_similarity(PATHWAYS1,PATHWAYS2)
			#print (jaccard_similarity_GENE)


			wt1=0
			wt2=0
			wt3=0
			
			wt1=1 if jaccard_similarity_GENE >0.5 else 0

			#wt2=1 if jaccard_similarity_SE >0.5 else 0

			wt3=1 if jaccard_similarity_PATHWAYS >0.5 else 0

			edgeWt = wt1+wt2+wt3
			#print (jaccard_similarity_GENE)
			if(edgeWt >0.2):
				G.add_edge(disease1,disease2,weight=edgeWt)
				print (disease1,disease2,edgeWt)




for drug in drugKeys:
	#print (drug1)
	gene1 = kegg_dict_drugs_features[drug]['GENE']
	PATHWAYS1=kegg_dict_drugs_features[drug]['PATHWAY']

	for disease in diseaseKeys:
		#print (drug1)
		gene2 = kegg_dict_disease_features[disease]['GENE']
		#SE1 = kegg_dict_disease_features[drug1]['SE']
		PATHWAYS2=kegg_dict_disease_features[disease]['PATHWAY']

		jaccard_similarity_GENE = jaccard_similarity(gene1,gene2)
#		jaccard_similarity_SE = jaccard_similarity(SE1,SE2)
		jaccard_similarity_PATHWAYS	 = jaccard_similarity(PATHWAYS1,PATHWAYS2)
		#print (jaccard_similarity_GENE)


		wt1=0
		wt2=0
		wt3=0
		
		wt1=1 if jaccard_similarity_GENE >0.3 else 0

		#wt2=1 if jaccard_similarity_SE >0.5 else 0

		wt3=1 if jaccard_similarity_PATHWAYS >0.3 else 0

		edgeWt = wt1+wt2+wt3

		#print (jaccard_similarity_GENE)
		if(edgeWt >1):
			G.add_edge(drug,disease,weight=edgeWt)
			print (drug,disease,edgeWt)
#			print ('+',kegg_dict_drugs_features[drug] ,'-',kegg_dict_disease_features[disease])


print (G.number_of_nodes())

#x=1
#y=2
#G.add_edge(x,y,weight=len(ans))



print ( nx.connected_components(G))


nx.write_adjlist(G,"graphs/graph.adjlist")

nx.write_graphml(G, "graphs/cyto.gml")
nx.write_graphml(G, "graphs/test.graphml")
