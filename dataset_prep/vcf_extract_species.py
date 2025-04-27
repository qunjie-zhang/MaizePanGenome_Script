
# 从VCF文件中提取指定物种信

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

# vcf_file_path = 'mini.vcf'
vcf_file_path = file

with open(vcf_file_path,'r',encoding='utf8') as vcf_file:
    vcf_file = vcf_file.read().splitlines()

# 所需物种名称列表 42
species_list = ["PH207","A632","OH43","Huangzaosi","Dan340","Ye478","zheng58","Chang72","Jing92","Jing724","Xu178","Mo17","Qi319"]
# 获取标题行
title_line = None
for line in vcf_file:
    if line.startswith('#CHROM'):
        title_line = line
        break
# 标题行 分割列表
title_line = title_line.split()
# 物种位置索引
species_index = list()


for species in species_list:
    if species not in title_line:
        print('Error: 未找到物种',species)
        exit()
    else:
        species_index.append(title_line.index(species))


with open(output,'w',encoding='utf8') as output_file:
    # 按行处理vcf文件
    for line in vcf_file:
        # 处理标题行
        if line.startswith('#CHROM'):
            line = line.split()
            # 临时列表,用于存放当前行编辑数据
            t = list()
            t = t + line[0:9]
            for index in species_index:
                t.append(line[index])
            output_file.write('\t'.join(t) + '\n')
            continue
        # 注释行直接添加
        if line.startswith('#'):
            output_file.write(line + '\n')
            continue

        # 数据行进行分割
        line = line.split()
        t = list()
        t = t + line[0:9]
        for index in species_index:
            t.append(line[index])
        output_file.write('\t'.join(t) + '\n')

print('Done!')
