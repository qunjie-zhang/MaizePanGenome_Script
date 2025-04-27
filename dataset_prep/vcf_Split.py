#!/usr/bin/env python
#-*- coding: utf-8 -*-

'''
 VCF文件按照染色体进行拆分，并用染色体名称命名。

python {self}.py -i input.vcf -o output

'''

import os,time
from tqdm import tqdm
from argparse import ArgumentParser

arg = ArgumentParser(description='Description - SCAU Bio_Lab')
arg.add_argument("-i",
                 "--input",
                 required=True,
                 help="input VCF file path")
arg.add_argument("-o",
                 "--output",
                 help="output path")

args = arg.parse_args()

input_vcf_path = args.input
output_vcf_path = args.output

# input_vcf_path = 'mini.vcf'
# output_vcf_path = 'output'


def llog(msg):
    msg = f'{time.asctime()} - {msg}'
    print(msg)

llog('正在读入VCF文件')
with open(input_vcf_path,'r',encoding='utf8') as vcf_file:
    vcf_file = vcf_file.read().splitlines()

llog('正在对VCF数据进行分类')

output_dict = dict()
for line in tqdm(vcf_file,desc='Process:',unit='line',ascii=True):
    # #开头行 为公共行 所有文件添加这部分内容
    if line.startswith('#'):
        if 'public' not in output_dict:
            output_dict['public'] = list()
            output_dict['public'].append(line)
            continue
        else:
            output_dict['public'].append(line)
            continue

    chr = line.split()[0]
    if chr not in output_dict:
        output_dict[chr] = list()
        output_dict[chr].append(line)
    else:
        output_dict[chr].append(line)

llog('VCF数据分类完成')

llog('正在输出分类文件')
# 如果输出路径不存在则创建该路径
if os.path.exists(output_vcf_path) == False:
    os.makedirs(output_vcf_path)

for chr_name,chr_context in tqdm(output_dict.items(),unit='file',ascii=True,desc='正在输出'):
    if chr_name == 'public':continue
    tqdm.write(f'{output_vcf_path}/{chr_name}.vcf')
    with open(f'{output_vcf_path}/{chr_name}.vcf','w',encoding='utf8') as out_file:
        # 写入公共表头数据
        for i in output_dict['public']:
            out_file.write(i+'\n')
        for i in chr_context:
            out_file.write(i+'\n')
llog('Done!')
