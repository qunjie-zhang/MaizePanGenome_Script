from tqdm import tqdm

'''
VCF 文件瘦身
去除多余的表头行，仅保留版本信息，格式信息以及列名。
INFO列缩减为 .

'''

from argparse import ArgumentParser

arg = ArgumentParser(description='Description - SCAU Bio_Lab')
arg.add_argument("-i",
                 "--input",
                 required=True,
                 help="input file path")
arg.add_argument("-o",
                 "--output",
                 help="output file path")

args = arg.parse_args()

file = args.input
output = args.output


with open(file,'r',encoding='utf8') as vcf_file:
    vcf_file = vcf_file.read().splitlines()

title_line = None
info_index = None

with open(output,'w',encoding='utf8') as output_file:

    for line in tqdm(vcf_file,desc='Processing VCF file'):

        if line.startswith('#CHROM'):
            output_file.write(line + '\n')
            title_line = line.split()
            info_index = title_line.index('INFO')
            continue
        if line.startswith('##fileformat'):
            output_file.write(line + '\n')
            continue
        if line.startswith('##FORMAT'):
            output_file.write(line + '\n')
            continue
        if line.startswith('#'):
            continue

        if title_line == None:
            print('Error: No title line found!')
            break

        line = line.split()
        line[info_index] = '.'
        output_file.write('\t'.join(line) + '\n')

print('Done!')