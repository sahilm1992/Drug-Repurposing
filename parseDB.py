import os
import csv
import gzip
import collections
import re
import io
import json
import xml.etree.ElementTree as ET

import requests
import pandas

from requests.auth import HTTPProxyAuth

proxyDict = { 
          'http'  : 'http://172.16.115.215:8888', 
          'https' : 'https://172.16.115.215:8888'
        }
auth = HTTPProxyAuth('smuser', 'Technical1_')
#response =  requests.get("http://git.dhimmel.com/uniprot/data/map/GeneID.tsv.gz", proxies=proxyDict, auth=auth)

#text = io.TextIOWrapper(gzip.GzipFile(fileobj=response.raw))
#print (text)
#r =
#print (r)


print ("WELCOME")
xml_path = os.path.join('download', 'drugbank.xml.gz')
with gzip.open(xml_path) as xml_file:
    tree = ET.parse(xml_file)
root = tree.getroot()
print ("WELCOME2")
ns = '{http://www.drugbank.ca}'
inchikey_template = "{ns}calculated-properties/{ns}property[{ns}kind='InChIKey']/{ns}value"
inchi_template = "{ns}calculated-properties/{ns}property[{ns}kind='InChI']/{ns}value"

rows = list()
for i, drug in enumerate(root):
    row = collections.OrderedDict()
    #print ()
    assert drug.tag == ns + 'drug'
    row['type'] = drug.get('type')
    row['drugbank_id'] = drug.findtext(ns + "drugbank-id[@primary='true']")
    row['name'] = drug.findtext(ns + "name")
    row['description'] = drug.findtext(ns + "description")
    row['groups'] = [group.text for group in
        drug.findall("{ns}groups/{ns}group".format(ns = ns))]
    row['atc_codes'] = [code.get('code') for code in
        drug.findall("{ns}atc-codes/{ns}atc-code".format(ns = ns))]
    row['categories'] = [x.findtext(ns + 'category') for x in
        drug.findall("{ns}categories/{ns}category".format(ns = ns))]
    row['inchi'] = drug.findtext(inchi_template.format(ns = ns))
    row['inchikey'] = drug.findtext(inchikey_template.format(ns = ns))

   # calcuProperty = [calculated_property_list_type.get('value') for calculated_property_list_type in 
    root2=drug.findall("{ns}calculated-properties/{ns}property".format(ns = ns))
    smilesList=[]
    for x in root2:
        flag =0 
        if(len(x) >0):
            kindValue =x.findtext(ns+'kind') 
               #print (len(kindValue),kindValue)

            if((kindValue)  == "SMILES"):
                   # print ('ya')
                    SMILES= x.findtext(ns+'value')
                    #print (SMILES)   
                    flag=1
                    smilesList.append(SMILES  )
                   # print ( row['SMILES']) 

            if(flag==0):
                    smilesList.append('')  
                    
   # Add drug aliases
    row['SMILES']=smilesList
    print (row['SMILES'])
    aliases = {
        elem.text for elem in 
        drug.findall("{ns}international-brands/{ns}international-brand".format(ns = ns)) +
        drug.findall("{ns}synonyms/{ns}synonym[@language='English']".format(ns = ns)) +
        drug.findall("{ns}international-brands/{ns}international-brand".format(ns = ns)) +
        drug.findall("{ns}products/{ns}product/{ns}name".format(ns = ns))

    }
    aliases.add(row['name'])
    row['aliases'] = sorted(aliases)

    rows.append(row)

alias_dict = {row['drugbank_id']: row['aliases'] for row in rows}
with open('./data/aliases.json', 'w') as fp:
    json.dump(alias_dict, fp, indent=2, sort_keys=True)

def collapse_list_values(row):
    for key, value in row.items():
        if isinstance(value, list):
            row[key] = '|'.join(value)
    return row

rows = list(map(collapse_list_values, rows))

columns = ['drugbank_id', 'name', 'type', 'groups', 'atc_codes', 'categories', 'inchikey', 'inchi', 'description','SMILES']
drugbank_df = pandas.DataFrame.from_dict(rows)[columns]
drugbank_df.head()


drugbank_slim_df = drugbank_df[
    drugbank_df.groups.map(lambda x: 'approved' in x) &
    drugbank_df.inchi.map(lambda x: x is not None) &
    drugbank_df.type.map(lambda x: x == 'small molecule')
]
drugbank_slim_df.head()




# write drugbank tsv
path = os.path.join('data', 'drugbank.tsv')
drugbank_df.to_csv(path, sep='\t', index=False,encoding='utf-8')

# write slim drugbank tsv
path = os.path.join('data', 'drugbank-slim.tsv')
drugbank_slim_df.to_csv(path, sep='\t', index=False,encoding='utf-8')

protein_rows = list()
for i, drug in enumerate(root):
    drugbank_id = drug.findtext(ns + "drugbank-id[@primary='true']")
    for category in ['target', 'enzyme', 'carrier', 'transporter']:
        proteins = drug.findall('{ns}{cat}s/{ns}{cat}'.format(ns=ns, cat=category))
        for protein in proteins:
            row = {'drugbank_id': drugbank_id, 'category': category}
            row['organism'] = protein.findtext('{}organism'.format(ns))
            row['known_action'] = protein.findtext('{}known-action'.format(ns))
            actions = protein.findall('{ns}actions/{ns}action'.format(ns=ns))
            row['actions'] = '|'.join(action.text for action in actions)
            uniprot_ids = [polypep.text for polypep in protein.findall(
                "{ns}polypeptide/{ns}external-identifiers/{ns}external-identifier[{ns}resource='UniProtKB']/{ns}identifier".format(ns=ns))]            
            if len(uniprot_ids) != 1:
                continue
            row['uniprot_id'] = uniprot_ids[0]
            ref_text = protein.findtext("{ns}references[@format='textile']".format(ns=ns))
            pmids = re.findall(r'pubmed/([0-9]+)', ref_text)
            row['pubmed_ids'] = '|'.join(pmids)
            protein_rows.append(row)

protein_df = pandas.DataFrame.from_dict(protein_rows)

# Read our uniprot to entrez_gene mapping
response =  requests.get("http://git.dhimmel.com/uniprot/data/map/GeneID.tsv.gz", proxies=proxyDict, auth=auth,stream=True)
print (response)

text = io.TextIOWrapper(gzip.GzipFile(fileobj=response.raw))
print (text)

uniprot_df = pandas.read_table(text, engine='python')
uniprot_df.rename(columns={'uniprot': 'uniprot_id', 'GeneID': 'entrez_gene_id'}, inplace=True)

# merge uniprot mapping with protein_dfp
entrez_df = protein_df.merge(uniprot_df, how='inner')

columns = ['drugbank_id', 'category', 'uniprot_id', 'entrez_gene_id', 'organism',
           'known_action', 'actions', 'pubmed_ids']
entrez_df = entrez_df[columns]

path = os.path.join('data', 'proteins.tsv')
entrez_df.to_csv(path, sep='\t', index=False,encoding='utf-8')



len(set(entrez_df.entrez_gene_id))
len(set(entrez_df.drugbank_id))



