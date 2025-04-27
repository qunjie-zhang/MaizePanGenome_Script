import os,json,time
import requests
from tqdm import tqdm

'''
通过 API 将有关基因信息推送到数据库中
'''

def log(msg):
    print ('[{}] {}'.format(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime()),msg))


gene_info_json_path = 'gene_info.json'
dataset_id = 7
api = 'http://12.12.12.138/api/haplotype/geneAdd'



with open(gene_info_json_path,'r',encoding='utf8') as gene_info_json:
    gene_info_json = json.load(gene_info_json)

s = requests.session()

for gene_id,gene_info in gene_info_json.items():

    params = {
        'dataset_id': dataset_id,
        'gene_id': gene_id,
        'chr': gene_info['chr'],
        'start': gene_info['start'],
        'end': gene_info['end'],
        'check':'1'
    }
    res = s.post(api,params=params).text

    try:
        res = json.loads(res)
    except:
        log('Error: json.loads failed')
        log(res)
        exit()

    if res['status'] == True:
        log(f'[SUCCES]: {gene_id}')
    else:
        log(f'[FAILED]: {gene_id} MSG: {res["msg"]}')


