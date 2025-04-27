from argparse import ArgumentParser


'''
用于提取指定VCF文件区段到指定位置

python vcf_extract_location.py -v mini.vcf -o mini_extract.vcf -c chr02 -s 1000 -e 2000
'''


arg = ArgumentParser(description='Description - SCAU Bio_Lab')
arg.add_argument("-v",
                 "--vcf",
                 required=True,)
arg.add_argument("-o",
                 "--output",
                 )
arg.add_argument("-c",
                 "--chr",
                 )
arg.add_argument("-s",
                 "--start",
                 )

arg.add_argument("-e",
                 "--end",
                 )
args = arg.parse_args()

vcf_file_path = args.vcf
output_file_path = args.output
chr = args.chr
start = args.start
end = args.end

with open(vcf_file_path,'r',encoding='utf8') as vcf_file:
    vcf_file = vcf_file.read().splitlines()

with open(output_file_path,'w',encoding='utf8') as output_file:
    for line in vcf_file:
        if line.startswith('#'):
            output_file.write(line + '\n')
            continue
        line = line.split()
        if line[0] == chr and int(line[1]) >= int(start) and int(line[1]) <= int(end):
            output_file.write('\t'.join(line) + '\n')

print('Done!')