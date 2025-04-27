

gff_file_path = 'Nip.gff3'

output_file_path = 'output.gff3'

'''
用于处理GFF文件中 行首以及行尾 多余的空格或者tab
尤其是首行 GFF 版本信息，该行不能缺少或出现多余的字符，会导致 genehapR 中的 import_gff 函数异常
'''

with open(gff_file_path,'r',encoding='utf8') as gff_file:
    gff_file = gff_file.read().splitlines()


with open(output_file_path,'w',encoding='utf8') as output_file:

    for line in gff_file:
        line = line.strip()
        output_file.write(line+'\n')