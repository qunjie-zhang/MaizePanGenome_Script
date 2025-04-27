
# 用于gff文件中 染色体名称的批量替换
from tqdm import tqdm

gff_file_path= 'Zea_mays.Zm-B73-REFERENCE-NAM-5.0.57.gff3'
output_file_path = 'Zea_mays.Zm-B73-REFERENCE-NAM-5.0.57_v2.gff3'


with open(gff_file_path, 'r') as gff_file:
    gff_file = gff_file.read().splitlines()

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

with open(output_file_path,'w',encoding='utf8') as output:
    output.write('##gff-version 3\n')
    for line in tqdm(gff_file,desc='Processing'):
        if line.startswith('#'):
            output.write(line+'\n')
            continue
        line = line.split()
        line[0] = map.get(line[0],line[0])
        output.write('\t'.join(line)+'\n')
print('Done')