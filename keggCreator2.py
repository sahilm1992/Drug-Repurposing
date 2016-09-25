#OM

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
		dict_kegg_dr_gene[disease] = set()
	dict_kegg_dr_gene[disease].add(gene)










