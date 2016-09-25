#OM


import pickle
def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)



def removeJunk(string):
	string = string.replace('\n','')
	string = string.strip()
	return string
strEntry="ENTRY       "

strName="NAME        "
strPathway="PATHWAY     "
strGene="GENE        "
strDrug="DRUG        "

disease_file="download/Kegg/disease"



f = open(disease_file)

dict_disease={}
skip=0
line = f.readline()
count =0

while line:
#   print line
	if ( strEntry in line) :
		entryLineList= line.split(strEntry)
		entryName = entryLineList[1].split('                      ')[0]
		if(entryName not in dict_disease):
			dict_disease[entryName] = {}
			dict_disease[entryName]['PATHWAY'] =set()
			dict_disease[entryName]['GENE'] =set()						
			dict_disease[entryName]['DRUG'] =set()						
		print (entryName)
	elif (strPathway in line ):
		pathwaySplit = line.split(strPathway)[1]
		pathwayId=pathwaySplit.split('  ')[0]
		pathwayId=pathwayId.split('(')[0]
		line = f.readline()
		print ('\t',pathwayId)		
		dict_disease[entryName]['PATHWAY'].add(removeJunk(pathwayId))
		while('            hsa' in line):
			pathwaySplit = line.split('            ')[1]
			pathwayId=pathwaySplit.split('  ')[0]
			pathwayId=pathwayId.split('(')[0]
			line = f.readline()		
			count=count+1
			dict_disease[entryName]['PATHWAY'].add(removeJunk(pathwayId)			)
			print ('\t', ':',pathwayId)
		continue
	elif (strGene in line):
		geneSplit = line.split(strGene)[1]
		geneId=geneSplit#geneSplit.split('(')[0]
		dict_disease[entryName]['GENE'].add(removeJunk(geneId)	)
		print ('\t',geneId)
		line = f.readline()		
		while('            ' in line):
			geneSplit = line.split('            ')[1]
			geneId=geneSplit#geneSplit.split('(')[0]
			line = f.readline()
			count=count+1
			print ('\t',':',geneId)
			dict_disease[entryName]['GENE'].add(removeJunk(geneId)	)
		continue
	elif (strDrug in line):
		drugValue = line.split(strDrug)[1]
		drugId=drugValue.split('  ')[0]
		dict_disease[entryName]['DRUG'].add(removeJunk(drugId)	)
		print ('\t',drugId)
		line = f.readline()		
		while('            ' in line):
			drugValue = line.split('            ')[1]
			drugId=drugValue.split('  ')[0]
			line = f.readline()
			count=count+1
			print ('\t',':',drugId)
			dict_disease[entryName]['DRUG'].add(removeJunk(drugId)	)
		continue
		
	line = f.readline()
#	print (line)
f.close()


drug_file="download/Kegg/drug"
f = open(drug_file)

dict_drugs={}
skip=0
line = f.readline()
count =0

strEntry="ENTRY       "

strName="NAME        "
strPathway="  PATHWAY   "
#strGene="GENE        "
strTarget="TARGET      "
strDbLink="DBLINKS     "
tempEntry=""
while line :
#   print line
	if ( strEntry in line) :
		entryLineList= line.split(strEntry)
		entryName = entryLineList[1].split('                      ')[0]
		print (entryName)
		tempEntry=entryName	
		if(entryName not in dict_drugs):
			dict_drugs[entryName] = {}
			dict_drugs[entryName]['PATHWAY'] =set()						
	elif (strPathway in line ):
		pathwaySplit = line.split(strPathway)[1]
		pathwayId=pathwaySplit.split('  ')[0]
		pathwayId=pathwayId.split('(')[0]
		line = f.readline()
		print ('\t',pathwayId)		
		dict_drugs[entryName]['PATHWAY'].add(removeJunk(pathwayId))
		while('            hsa' in line):
			pathwaySplit = line.split('            ')[1]
			pathwayId=pathwaySplit.split('  ')[0]
			pathwayId=pathwayId.split('(')[0]
			line = f.readline()		
			count=count+1
			dict_drugs[entryName]['PATHWAY'].add(removeJunk(pathwayId))
			print ('\t', ':',pathwayId)
		continue
	elif (strTarget in line):
		geneSplit = line.split(strTarget)[1]
		geneId=geneSplit#geneSplit.split('(')[0]
		print ('\t',geneId)
		line = f.readline()		
		while('            ' in line):
			geneSplit = line.split('            ')[1]
			geneId=geneSplit#geneSplit.split('(')[0]
			line = f.readline()
			count=count+1
			print ('\t',':',geneId)
		continue
	elif (strDbLink in line):
#		print ("ya")
		dbLinks = line.split(strDbLink)[1]
		dbLinks=dbLinks#geneSplit.split('(')[0d]
		
		if('DrugBank' in dbLinks):
			dict_drugs[entryName]['DBID'] =removeJunk(dbLinks)
			print ('\t',':',dbLinks)
	
		line = f.readline()		
		while('            ' in line):
			dbLinks = line.split('            ')[1]
			dbLinks=dbLinks#geneSplit.split('(')[0]
			line = f.readline()
			count=count+1
			if('DrugBank' in dbLinks):
				dict_drugs[entryName]['DBID'] =removeJunk(dbLinks)
				print ('\t',':',dbLinks)
	
		continue
				
	line = f.readline()
#	print (line)

save_obj(dict_drugs,'kegg_dict_drugs')
save_obj(dict_disease,'kegg_dict_disease')

for dr in dict_disease:
	print (dr,dict_disease[dr])

print (len(dict_disease))
print (len(dict_drugs))
#print (dict_drugs)
#print (dict_disease)

f.close()
