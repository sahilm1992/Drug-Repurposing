import pandas


tsv_pubChem_Drug_Bank = pandas.read_table('data/pubchem-mapping.tsv')

#tsv_pubChem_Drug_Bank= tsv_pubChem_Drug_Bank

drugBankIdCols=tsv_pubChem_Drug_Bank[tsv_pubChem_Drug_Bank.columns[0]] 
pubChemIdBankCols=tsv_pubChem_Drug_Bank[tsv_pubChem_Drug_Bank.columns[1]] 

records = len( tsv_pubChem_Drug_Bank)
#print (records[1])

dict_pubChem_Drug_Bank={}

for i in range(0,records):
	None	
	drugBank = str(drugBankIdCols[i])
	pubChemId = str(pubChemIdBankCols[i])
	dict_pubChem_Drug_Bank[pubChemId]= drugBank

#print (dict_pubChem_Drug_Bank)

#for x in dict_pubChem_Drug_Bank.keys():
#	print (x, dict_pubChem_Drug_Bank[x])


drug_Features={}
drug_Features['GENES']=set()
drug_Features['SE']=set()


tsv_drug_protein = pandas.read_table('data/proteins.tsv')

drug_Ids = tsv_drug_protein[tsv_drug_protein.columns[0]] 
gene_uniprot_id = tsv_drug_protein[tsv_drug_protein.columns[2]] 





