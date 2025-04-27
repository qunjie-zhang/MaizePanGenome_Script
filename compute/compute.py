
import time,json,os
import subprocess
from BioinfodataPlatformWeb import BioinfodataPlatformWeb


if __name__ == '__main__':

    # 网站交互工具类
    bpw = BioinfodataPlatformWeb(url='http://12.12.12.138')
    app = 'haplotype'

    while True:
        # 生产参数
        job_info = bpw.getOneJob(app)

        if job_info == None:
            time.sleep(10)
            continue

        # 系统参数
        job_id = job_info['id']
        job_title = job_info['title']
        create_time = job_info['create_time']
        start_time = job_info['start_time']
        end_time = job_info['end_time']
        user_id = job_info['user_id']
        parameter = job_info['parameter']
        status = job_info['status']
        # 作业参数
        chr = parameter['chr']

        start = parameter['start']
        end = parameter['end']

        pheo_file = parameter['pheo_file']
        acc_file = parameter['acc_file']

        prefix = parameter['prefix']
        gene_id = parameter['gene_id']
        na_drop = parameter['na_drop']
        dataset_path = parameter['dataset_path']
        custom_pos = parameter['custom_pos']
        hetero_remove = parameter['hetero_remove']
        filter_vcf_mode = parameter['filter_vcf_mode']
        filter_vcf_type = parameter['filter_vcf_type']
        output_path = parameter['output_path']

        range_file = parameter['range_file']
        # 物种信息列表
        species_list = list()
        # 位置区段列表
        locations_list = list()
        # 位置区段最大最小值
        locations_max = 0
        locations_min = 0


        # 任务日志打印并保存文件
        def log(msg):
            res = f'{time.asctime()}: {msg}'
            with open(os.path.join(output_path,'log.txt'),'a',encoding='utf8') as log_file:
                log_file.write(res + '\n')
            print(res)

        print('-----------------------------------------------------------------')
        log('接收任务：' + job_title)
        bpw.set_job_status(job_id,'running')

        # range_file 文件加载
        if range_file != None:
            log('已上传自定义位置信息配置文件')
            with open(os.path.join(output_path,range_file),'r',encoding='utf8') as range_file:
                range_file = json.load(range_file)
                species_list = range_file['species']
                chr = range_file['chr']

                locations_list = range_file['location']
                # 多个位置区段信息平展，用于获取最大值最小值
                locations_expand = [item for sublist in locations_list for item in sublist]
                locations_max = max(locations_expand)
                locations_min = min(locations_expand)
        else:
            # 若不存在range_file 则直接将 start end 参数写入位置列表
            # 若不存在 start end 则任务失败
            log('未上传自定义位置配置文件')
            if start != None and end != None:
                locations_list.append([start,end])
                locations_max = end
                locations_min = start
            else:
                log('位置信息不完整,任务失败')
                bpw.set_job_status(job_id,'failed')
                continue


        acc_file = os.path.join(output_path,parameter['acc_file']) if parameter['acc_file'] != None else None
        pheo_file = os.path.join(output_path,parameter['pheo_file']) if parameter['pheo_file'] != None else None
        range_file = os.path.join(output_path,parameter['range_file']) if parameter['range_file'] != None else None

        # 加载数据集配置文件
        dataset_info = os.path.join(dataset_path,'dataset.json')
        if not os.path.exists(dataset_path):
            log('数据集路径或配置文件不存在,任务失败')
            bpw.set_job_status(job_id,'failed')
            continue
        with open(dataset_info,'r',encoding='utf8') as dataset_info_file:
            dataset_info = json.load(dataset_info_file)


        vcf_file_path = os.path.join(dataset_path,dataset_info['chr_vcf'],f'{chr}.vcf')
        gff_file_path = os.path.join(dataset_path,dataset_info['gff'])
        acc_file_path = os.path.join(dataset_path,dataset_info['acc'] if acc_file == None else os.path.join(output_path,acc_file))
        pheo_file_path = os.path.join(dataset_path,dataset_info['pheo'] if pheo_file == None else os.path.join(output_path,pheo_file))

        with open(vcf_file_path,'r',encoding='utf8') as vcf_file:
            vcf_file = vcf_file.read().splitlines()

        # 获取标题行内容
        title_line = None
        for line in vcf_file:
            if line.startswith('#CHROM'):
                title_line = line.split()
                break

        if title_line == None:
            log('VCF表头信息获取失败')
            bpw.set_job_status(job_id,'failed')
            continue


        # VCF 所需列位置索引信息
        column_index_list = list()
        # 判断物种列表是否为空。若为空则使用全部物种信息
        if len(species_list) != 0:
            for species in species_list:
                if species in title_line:
                    column_index_list.append(title_line.index(species))
        else:
            # 若为空，则将全部列索引写入
            for i in range(len(title_line)):
                if i >= 9:
                    column_index_list.append(i)

        with open(os.path.join(output_path,'work.vcf'),'w',encoding='utf8') as work_vcf:
            work_vcf.write('##fileformat=VCFv4.2\n')
            # 处理表头行相关信息
            for line in vcf_file:
                if line.startswith("#"):
                    if line.startswith('#CHROM'):
                        line = line.split()
                        output_line = line[:9]
                        for index in column_index_list:
                            output_line.append(line[index])
                        # 将处理好的表头写入文件
                        work_vcf.write('\t'.join(output_line) + '\n')
                        break

            # 处理VCF每一行数据
            for line in vcf_file:
                if line.startswith('#'):continue
                line = line.split()
                # 前九列为固定列
                output_line = line[:9]

                # 非目标染色体信息
                if line[0] != chr:continue

                for locations in locations_list:
                    start = int(locations[0])
                    end = int(locations[1])
                    # 若当前染色体名称匹配 且 位置在范围内
                    if int(line[1]) >= start and int(line[1]) <= end:
                        for index in column_index_list:
                            output_line.append(line[index])
                        work_vcf.write('\t'.join(output_line) + '\n')
                        break
        log('VCF文件提取完成')

        # 将网站根目录映射到此
        cmd = [f'Rscript /biodataplatform/compute/hap.version2_ctl.v2.r',
               f'--geneid "{job_title if gene_id == None else gene_id}"',
               f'--chr {chr}',
               f'--start {locations_min}',
               f'--end {locations_max}',
               f'--happrefix {prefix}',
               f'--vcf {os.path.join(output_path,"work.vcf")}',
               f'--gff {gff_file_path}',
               f'--pheno {pheo_file_path}',
               f'--accinfo {acc_file_path}',
               f'--hetero_remove {hetero_remove}',
               f'--na_drop {na_drop}',
               f'--filter_vcf_mode {'both' if filter_vcf_mode else 'POS'}',
               f'--filter_vcf_type {','.join(filter_vcf_type)}',
               f'--output {output_path}']
        cmd = ' '.join(cmd)

        log('EXEC:' + cmd)

        try:
            t = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8",timeout=1800)
        except subprocess.TimeoutExpired:
            log('任务超时结束')
            bpw.set_job_status(job_id, 'failed')
            continue
        except Exception as e:
            log('任务异常:' + str(e))
            bpw.set_job_status(job_id, 'failed')
            continue

        # 写入本次任务日志文件
        log(t.stdout)

        if t.returncode == 0:
            log(f'[任务完成] ID:{job_id}  TITLE:{job_title} OUTPUT: {output_path}')
            bpw.set_job_status(job_id, 'success')
        else:
            log(f'[任务失败] ID:{job_id}  TITLE:{job_title} OUTPUT: {output_path}')
            bpw.set_job_status(job_id, 'failed')

        log('Task Done.')
