import os,re,json

'''
该脚本用于从GFF文件中提取基因相关信息。
包括基因的名称，所在染色体位置，开始结束位置等。

该脚本需手动更改 gff 文件名中 基因名称提取的正则表达式

'''

gff_path = 'Zea_mays.Zm-B73-REFERENCE-NAM-5.0.57_v2.gff3'

with open(gff_path,'r',encoding='utf8') as gff_file:
    gff_file = gff_file.read().splitlines()
    
gene_dict = dict()
chr_list = list()

# 获取基因ID正则表达式
regex_geneid = re.compile('ID=gene:(.*?);')

for line in gff_file:
    if line.startswith('#') or line.isspace() or line == '':continue
    line = line.split('\t')
    if line[2] == 'gene':
        gene_id = regex_geneid.search(line[8]).group(1)
        try:
            chr = line[0]
            start = int(line[3])
            end = int(line[4])
        except:
            print('Error: chr or start or end is not int')
            print(line)
            exit(1)
        if chr not in chr_list:
            chr_list.append(chr)

        if gene_id in gene_dict:
            print('Error: gene_id repeat')
        gene_dict[gene_id] = {
            'chr':chr,
            'start':start,
            'end':end
        }
with open('gene_info.json','w',encoding='utf8') as gene_info_file:
    gene_info_file.write(json.dumps(gene_dict,indent=4,ensure_ascii=False))
print('Done.')