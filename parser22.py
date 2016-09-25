import pandas


import pickle
def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)




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
#drug_Features['GENES']=set()
#drug_Features['SE']=set()


tsv_drug_protein = pandas.read_table('data/proteins.tsv')

drug_Ids = tsv_drug_protein[tsv_drug_protein.columns[0]] 
gene_uniprot_id = tsv_drug_protein[tsv_drug_protein.columns[2]] 

for i in range(0,len(drug_Ids)):
	None	
	drugBank = str(drug_Ids[i])
	uniProt = str(gene_uniprot_id[i])
	if(drugBank not in drug_Features):
		drug_Features[drugBank]= {}
		drug_Features[drugBank]['GENES']= set()
		drug_Features[drugBank]['SE']= set()
			
	drug_Features[drugBank]['GENES'].add(uniProt)
	#drug_Features[drugBank]['SE'].add(uniProt)
	
	#dict_pubChem_Drug_Bank[pubChemId]= drugBank

'''for x in drug_Features.keys():
	if(len(drug_Features[x]['GENES'])>50):
		print (x, drug_Features[x])
'''


#drug SIDE EFFECT

tsv_drugpubChem_SE = pandas.read_table('data/meddra_all_indications.tsv')

drug_pubchemIds = tsv_drugpubChem_SE[tsv_drugpubChem_SE.columns[0]] 
se_Id = tsv_drugpubChem_SE[tsv_drugpubChem_SE.columns[1]] 


#print (dict_pubChem_Drug_Bank)

for i in range(0,len(drug_pubchemIds)):
#	print (i)
	pubChemId = str(drug_pubchemIds[i])
	pubChemId = pubChemId[3:] 	# CID 100000085
	pubChemId_int = int(pubChemId) # Subtract 1000
	pubChemId_int = str(pubChemId_int - 100000000)
	seId = str(se_Id[i])
	#print (pubChemId,pubChemId_int)
	if (pubChemId_int in dict_pubChem_Drug_Bank):
		#print (pubChemId_int)	
		dbId = 	str(dict_pubChem_Drug_Bank[pubChemId_int])
#		print (dbId,seId)
		if(dbId not in drug_Features):
			drug_Features[dbId]= {}
			drug_Features[dbId]['GENES']= set()
			drug_Features[dbId]['SE']= set()
	
		
		drug_Features[dbId]['SE'].add(seId)
#		print (drug_Features[dbId],seId)
		
		
'''
for x in drug_Features.keys():
	if(len(drug_Features[x]['SE']) >30  and len(drug_Features[x]['GENES'] ) >30 ):
		print (x, drug_Features[x])
'''

save_obj(dict_pubChem_Drug_Bank,'dict_pubChem_Drug_Bank')
save_obj(drug_Features,'drug_Features')





#print (drug_Features)
	

