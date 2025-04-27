from tqdm import tqdm

vcf_file_path = 'm350_new_from_zlx.vcf'

vcf_file_output_path = 'm350_new_from_zlx.v2.vcf'

map = {
    '1':'chr01',
    '2':'chr02',
    '3':'chr03',
    '4':'chr04',
    '5':'chr05',
    '6':'chr06',
    '7':'chr07',
    '8':'chr08',
    '9':'chr09',
    '10':'chr10',
}


with open(vcf_file_path,'r',encoding='utf8') as vcf_file:
    vcf_file = vcf_file.read().splitlines()

with open(vcf_file_output_path,'w',encoding='utf8') as output_file:

    for line in tqdm(vcf_file,desc='Processing:',unit='line',ascii=True):

        if line.startswith('#'):
            for _ in ['##fileformat','##FILTER','#CHROM']:
                if line.startswith(_):
                    output_file.write(line+'\n')
                    break
            continue

        line = line.split()
        line[0] = map.get(line[0],line[0])

        output_file.write('\t'.join(line)+'\n')

print('Done!')
