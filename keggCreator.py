#OM
import pandas
import pickle
def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


tsv_kegg_dr_genes = pandas.read_csv('download/Kegg/kegg-dr-gene.csv')

#tsv_pubChem_Drug_Bank= tsv_pubChem_Drug_Bank

keggDrugCols=tsv_kegg_dr_genes[tsv_kegg_dr_genes.columns[2]] 
keggGeneSymbolCols=tsv_kegg_dr_genes[tsv_kegg_dr_genes.columns[1]] 

records = len(keggDrugCols)
#print (records[1])

dict_kegg_dr_gene={}

for i in range(0,records):
	None	

	drug = str(keggDrugCols[i])
	gene = str(keggGeneSymbolCols[i])

	if ( drug not in dict_kegg_dr_gene):
		dict_kegg_dr_gene[drug] = set()
	dict_kegg_dr_gene[drug].add(gene)

			 
#	dict_pubChem_Drug_Bank[pubChemId]= drugBank



tsv_kegg_genes_dis = pandas.read_csv('download/Kegg/kegg-gene-dis.csv')

keggDiseaseCols=tsv_kegg_genes_dis[tsv_kegg_genes_dis.columns[2]] 
keggGeneSymbolCols=tsv_kegg_genes_dis[tsv_kegg_genes_dis.columns[1]] 

records = len(keggDiseaseCols)
#print (records[1])

dict_kegg_gene_dis={}

for i in range(0,records):
	None	

	disease = str(keggDiseaseCols[i])
	gene = str(keggGeneSymbolCols[i])
#	dict_pubChem_Drug_Bank[pubChemId]= drugBank

	if ( disease not in dict_kegg_gene_dis):
		dict_kegg_gene_dis[disease] = set()
	dict_kegg_gene_dis[disease].add(gene)


#print (dict_kegg_dr_gene)

#print (dict_kegg_gene_dis)


'''for x in dict_kegg_dr_gene.keys():
	print (x,dict_kegg_dr_gene[x])


for x in dict_kegg_gene_dis.keys():
	print (x,dict_kegg_gene_dis[x])



'''

kegg_dict_drugs=load_obj('kegg_dict_drugs')
kegg_dict_disease=load_obj('kegg_dict_disease')

#for x in kegg_dict_disease.keys():
#	print (x,kegg_dict_disease[x])


for disease in dict_kegg_gene_dis.keys():	
	
	#print (disease)
	#print(disease,dict_kegg_gene_dis[disease]['PATHWAY'])
	if(disease in 	kegg_dict_disease ):
#		print (kegg_dict_disease[disease]['GENE'])
#		print('+',dict_kegg_gene_dis[disease])
				
		(kegg_dict_disease[disease]['GENE'])=(kegg_dict_disease[disease]['GENE'])|((dict_kegg_gene_dis[disease]))


for drug in dict_kegg_dr_gene.keys():	
#	print ('+',drug)
	if (drug in kegg_dict_drugs):
		#print(drug,kegg_dict_drugs[drug]['PATHWAY'])
		kegg_dict_drugs[drug]['GENE'] =kegg_dict_drugs[drug]['GENE'] | ((dict_kegg_dr_gene[drug])		)
			


for x in kegg_dict_disease.keys():
	#if( kegg_dict_disease[x]['DBID'] is not None):
		if(len(kegg_dict_disease[x]['PATHWAY']) >0  or len(kegg_dict_disease[x]['GENE']) >0 ):		
			None
			print (x,kegg_dict_disease[x])



toDelete =[]
for x in kegg_dict_drugs.keys():
	if( kegg_dict_drugs[x]['DBID'] is not None):
		if(len(kegg_dict_drugs[x]['PATHWAY']) >0  or len(kegg_dict_drugs[x]['GENE']) >0 ):
			None			
			print (x,kegg_dict_drugs[x])
	else:
		toDelete.append(x)

for d in toDelete:
	del (kegg_dict_drugs[d])

save_obj(kegg_dict_drugs,'kegg_dict_drug_features')
save_obj(kegg_dict_disease,'kegg_dict_disease_features')




