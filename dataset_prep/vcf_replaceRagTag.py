from tqdm import tqdm

vcf_file_path = 'rice43_graph.vcf'

vcf_file_output_path = 'rice43_graph.v2.vcf'

with open(vcf_file_path,'r',encoding='utf8') as vcf_file:
    vcf_file = vcf_file.read().splitlines()

with open(vcf_file_output_path,'w',encoding='utf8') as output_file:

    for line in tqdm(vcf_file,desc='Processing:',unit='line',ascii=True):

        if line.startswith('#'):
            output_file.write(line+'\n')
            continue

        line = line.split()
        line[0] = line[0].replace('_RagTag','')

        output_file.write('\t'.join(line)+'\n')

print('Done!')
