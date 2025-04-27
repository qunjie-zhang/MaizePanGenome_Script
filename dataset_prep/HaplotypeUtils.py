import pydoc
import time,os
from tqdm import tqdm
import re,json
import requests

from multiprocessing import Pool

'''
数据分析平台 haplotype 模块 数据集自动构建流程
 '''
class HaplotypeUtils:

    # 初始化操作
    def __init__(self,**kwargs):
        pass

    '''
    VCF 文件按照染色体进行拆分
    '''
    def vcf_split(self,vcf_path,vcfs_dir):
        with open(vcf_path, 'r', encoding='utf8') as vcf_file:
            vcf_file = vcf_file.read().splitlines()

        # 用于输出的VCF文件
        output_dict = dict()
        # 解析VCF文件按照染色体进行拆分保存在字典中
        for line in tqdm(vcf_file, desc='正在识别:', unit='line', ascii=True):
            # 注释行进行统一存储在 public 标识内
            if line.startswith('#'):
                if 'public' not in output_dict:output_dict['public'] = list()
                output_dict['public'].append(line)
                continue
            # 染色体名称标识
            chr = line.split()[0]
            if chr not in output_dict:
                output_dict[chr] = list()
                output_dict[chr].append(line)
            else:
                output_dict[chr].append(line)

        # 输出文件
        if os.path.exists(vcfs_dir) == False:
            os.makedirs(vcfs_dir)
        for chr_name, chr_context in tqdm(output_dict.items(), unit='file', ascii=True, desc='正在输出'):
            if chr_name == 'public': continue
            with open(f'{vcfs_dir}/{chr_name}.vcf', 'w', encoding='utf8') as out_file:
                # 写入公共表头数据
                for i in output_dict['public']:
                    out_file.write(i + '\n')
                for i in chr_context:
                    out_file.write(i + '\n')
        return True

    '''
    VCF 文件瘦身
    - 删除多余表头信息
    - INFO列(py索引7) 信息重置为 .
    '''
    def vcf_slim(self,input_vcf_path,output_vcf_path):

        with open(input_vcf_path,'r',encoding='utf8') as input_vcf_file:
            input_vcf_file = input_vcf_file.read().splitlines()

        with open(output_vcf_path,'w',encoding='utf8') as output_vcf_file:
            # 关键表头字段
            head_keyword_list = ['#CHROM','##fileformat','##FORMAT']
            for input_vcf_file_line in input_vcf_file:
                # 注释行信息处理
                if input_vcf_file_line.startswith('#'):
                    # 将关键表头信息写入文件
                    for keyword in head_keyword_list:
                        if input_vcf_file_line.startswith(keyword):
                            output_vcf_file.write(input_vcf_file_line+'\n')
                            break
                    continue

                input_vcf_file_line = input_vcf_file_line.split()

                # 将INFO字段中的内容进行处理
                input_vcf_file_line[7] = '.'
                output_vcf_file.write('\t'.join(input_vcf_file_line)+'\n')
        return True

    # 按照给定物种列表输出VCF文件
    def vcf_extract_species(self,input_vcf_path,output_vcf_path,species_list=[]):
        with open(input_vcf_path,'r',encoding='utf8') as input_vcf_file:
            input_vcf_file = input_vcf_file.read().splitlines()
        with open(output_vcf_path,'w',encoding='utf8') as output_vcf_file:

            # 定义标题行属性
            title_line = None

            for input_vcf_file_line in input_vcf_file:
                if input_vcf_file_line.startswith('#CHROM'):
                    title_line = input_vcf_file_line
                    break

            if title_line == None:
                print('VCF文件缺失表头标题行！ -> #CHROM...')
                exit()

            title_line = title_line.split()

            if len(title_line) <= 9:
                print('VCF文件表头标题行字段不足！ -> #CHROM...')
                print(title_line)
                exit()

            # 所需物种在vcf文件中的列索引
            species_index_list = list()

            for species in species_list:
                if species not in title_line:
                    print(species,'物种不存在于VCF文件中！')
                    continue
                species_index_list.append(title_line.index(species))

            for input_vcf_file_line in input_vcf_file:
                # 标题行单独处理一遍
                if input_vcf_file_line.startswith('#CHROM'):
                    input_vcf_file_line = input_vcf_file_line.split()
                    output_line = input_vcf_file_line[:9]
                    for index in species_index_list:
                        output_line.append(input_vcf_file_line[index])
                    output_vcf_file.write('\t'.join(output_line) + '\n')
                    continue
                elif input_vcf_file_line.startswith('#'):
                    output_vcf_file.write(input_vcf_file_line+'\n')
                    continue

                input_vcf_file_line = input_vcf_file_line.split()
                output_line = input_vcf_file_line[:9]

                # 将所需索引位置信息写入输出信息
                for index in species_index_list:
                    output_line.append(input_vcf_file_line[index])
                output_vcf_file.write('\t'.join(output_line)+'\n')
        return True

    '''
    修正VCF文件格式
    - 将图泛VCF文件模拟为标准VCF文件格式。
    '''
    def vcf_graph2stand(self,input_vcf_path,output_vcf_path):
        with open(input_vcf_path,'r',encoding='utf8') as input_vcf_file:
            input_vcf_file = input_vcf_file.read().splitlines()
        with open(output_vcf_path, 'w', encoding='utf8') as out_file:
            for line in input_vcf_file:
                if line.startswith('#'):
                    out_file.write(line + '\n')
                    continue
                line = line.split()
                for index, i in enumerate(line[9:]):
                    line[index + 9] = f'{i}/{i}'
                out_file.write('\t'.join(line) + '\n')
        return True

    '''
    Gff 文件过滤清理
    - 清除前后多余的空格或TAB。
    '''
    def gff_clean(self,input_gff_path,output_gff_path):
        with open(input_gff_path, 'r', encoding='utf8') as gff_file:
            gff_file = gff_file.read().splitlines()

        with open(output_gff_path, 'w', encoding='utf8') as output_file:
            for line in gff_file:
                line = line.strip()
                output_file.write(line + '\n')
        return True

    '''
    从GFF文件中提取基因信息到json文件
    '''
    def gff_getGeneList(self,input_gff_path,output_json_path):

        with open(input_gff_path, 'r', encoding='utf8') as gff_file:
            gff_file = gff_file.read().splitlines()

        # 根据第一个gene行出现的信息判断本次过滤基因ID使用的正则表达式
        for line in gff_file:
            if line.startswith('#') or line.isspace() or line == '': continue
            line = line.split()
            if line[2] == 'gene':
                if line[8].startswith('ID=gene:'):
                    regex_geneid = re.compile('ID=gene:(.*?);')
                    break
                elif line[8].startswith('ID='):
                    regex_geneid = re.compile('ID=(.*?);')
                    break
                else:
                    print('未注册的基因ID格式！')
                    exit()

        # 用于存储基因信息
        gene_dict = dict()

        for line in gff_file:
            if line.startswith('#') or line.isspace() or line == '': continue
            line = line.split()
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

                if gene_id in gene_dict:
                    print('Error: gene_id repeat')
                gene_dict[gene_id] = {
                    'chr': chr,
                    'start': start,
                    'end': end
                }

        with open(output_json_path, 'w', encoding='utf8') as output_file:
            json.dump(gene_dict, output_file, indent=4, ensure_ascii=False)
        return True


    '''
    将数据集信息上传到数据库
    '''
    def website_createDataset(self,url,dataset_name,dataset_path):
        # 请求数据集创建接口
        response = requests.post(url, data={
            'dataset_name': dataset_name,
            'dataset_path': dataset_path
        })
        return response.json()


    '''
    将数据集的基因信息写入数据库
    '''
    def website_createDatasetGene(self,url,dataset_id,gene_id,chr,start,end,check=True):
        response = requests.post(url, data={
            'dataset_id': dataset_id,
            'gene_id': gene_id,
            'chr': chr,
            'start': start,
            'end': end,
            'check': check
        })
        return response.json()

    '''
    # 初始化数据集信息
    '''
    def website_initDataset(self,vcfs_path,gff_path,acc_path,pheo_path,output_path):
        d = {
            'chr_vcf': 'Maize_138/',
            'gff': 'Zea_mays.Zm-B73-REFERENCE-NAM-5.0.57_v2.gff3',
            'acc': 'acc.csv',
            'pheo': 'pheo.csv',
        }
        with open(output_path, 'w', encoding='utf8') as output_file:
            json.dump(d, output_file, indent=4, ensure_ascii=False)


    # TEST
    def test(self,msg):
        return (msg)


if __name__ == "__main__":
    haplotype_utils = HaplotypeUtils()

    # # 定义原始VCF文件
    vcf_raw_file_path = 'syri_chrall.vcf'


    #### 1 VCF文件拆分
    vcf_split_dir = '01.vcf_split'
    haplotype_utils.vcf_split(vcf_raw_file_path,vcf_split_dir)

    ### 2 VCF文件瘦身
    pool = Pool(16)
    vcf_slime_dir = '02.vcfq_slim'
    if not os.path.exists(vcf_slime_dir):os.makedirs(vcf_slime_dir)

    for file in os.listdir(vcf_split_dir):
        input_vcf_path = f'{vcf_split_dir}/{file}'
        output_vcf_path = f'{vcf_slime_dir}/{file}'
        pool.apply_async(haplotype_utils.vcf_slim,args=(input_vcf_path,output_vcf_path))
    pool.close()
    pool.join()


    # ### 3 Graph VCF 2 Stand VCF
    # vcf_stand_dir = '03.vcf_stand'
    # if not os.path.exists(vcf_stand_dir):os.makedirs(vcf_stand_dir)

    # pool = Pool(16)
    # for file in os.listdir(vcf_slime_dir):
    #     input_vcf_path = f'{vcf_slime_dir}/{file}'
    #     output_vcf_path = f'{vcf_stand_dir}/{file}'
    #     pool.apply_async(haplotype_utils.vcf_graph2stand,args=(input_vcf_path,output_vcf_path))
    # pool.close()
    # pool.join()



    # ### 4 仅提取对应物种
    # species_list = [
    #     '02428',
    #     'ZS97',
    # ]
    # vcf_species_dir = '04.vcf_species'
    # if not os.path.exists(vcf_species_dir):os.makedirs(vcf_species_dir)
    # pool = Pool(16)
    # for file in os.listdir(vcf_stand_dir):
    #     input_vcf_path = f'{vcf_stand_dir}/{file}'
    #     output_vcf_path = f'{vcf_species_dir}/{file}'
    #     pool.apply_async(haplotype_utils.vcf_extract_species,args=(input_vcf_path,output_vcf_path,species_list))
    # pool.close()
    # pool.join()


