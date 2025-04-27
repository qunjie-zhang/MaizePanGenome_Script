from tqdm import tqdm
from argparse import ArgumentParser
'''
用于将图泛产生的VCF文件转换为标准VCF文件
主要将后面数据内容的 单个字符替换为 1/1 等形式
'''
arg = ArgumentParser(description='Description - SCAU Bio_Lab')
arg.add_argument("-i",
                 "--input",
                 required=True,
                 help="input file path")
arg.add_argument("-o",
                 "--output",
                 help="output file path")
args = arg.parse_args()
input = args.input
output = args.output

with open(input,'r',encoding='utf8') as vcf_file:
    vcf_file = vcf_file.read().splitlines()

with open(output,'w',encoding='utf8') as out_file:
    for line in tqdm(vcf_file,desc='Processing VCF file'):
        if line.startswith('#'):
            out_file.write(line + '\n')
            continue
        line = line.split()
        for index, i in enumerate(line[9:]):
            line[index + 9] = f'{i}/{i}'
        out_file.write('\t'.join(line) + '\n')
print('Done!')