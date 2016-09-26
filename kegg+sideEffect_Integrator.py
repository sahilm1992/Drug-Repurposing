#OM

import pickle

def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


kegg_dict_drugs_features=load_obj('kegg_dict_drug_features')
kegg_dict_disease_features=load_obj('kegg_dict_disease_features')




dict_sideEffectFeatures=load_obj('drug_Features')

for x in kegg_dict_drugs_features.keys():
#	print (x,kegg_dict_drugs_features[x])
	DBID = kegg_dict_drugs_features[x]['DBID']	
	if (DBID in dict_sideEffectFeatures and len(dict_sideEffectFeatures[DBID]['SE']) >0 ):
		kegg_dict_drugs_features[x]['SE'].add(str(dict_sideEffectFeatures[DBID]['SE'] ))
		#print (x, dict_sideEffectFeatures[DBID])	


for x in kegg_dict_drugs_features.keys():
	print (x,kegg_dict_drugs_features[x])


print ('#kegg_dict_drugs',len(kegg_dict_drugs_features))

print ('#kegg_dict_disease',len(kegg_dict_disease_features))


save_obj(kegg_dict_drugs_features,'kegg_dict_drug_features+SE')
save_obj(kegg_dict_disease_features,'kegg_dict_disease_features+SE')


#load drugbank ids

