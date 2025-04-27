
import json

'''
用于生成数据集定义信息配置文件
'''

d = {
    'chr_vcf':'vcfs',
    'gff': '../modify_MSU.IGDBv1.Allset.gff',
    'acc': '../acc.csv',
    'pheo': '../pheo.csv',
}

with open('dataset.json','w',encoding='utf8') as f:
    f.write(json.dumps(d,indent=4))
