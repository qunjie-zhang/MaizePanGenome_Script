
from urllib.parse import urljoin
import json,requests,os
import datetime

class BioinfodataPlatformWeb:
    def __init__(self,url='http://127.0.0.1'):
        self.base_url = url
        # 更新作业信息接口
        self.job_updateStatus_url = urljoin(self.base_url,'/api/job/updateInfo')
        # 获取作业信息接口
        self.job_getOneJob_url = urljoin(self.base_url,'/api/job/getOneJob')

    # 版本号
    def version(self):
        return '1.0.0'

    # 获取一个任务
    def getOneJob(self,app):
        getOneJobUrl = urljoin(self.job_getOneJob_url,'?app='+app)
        res = requests.post(getOneJobUrl).json()
        if res['status'] == True:
            return res['data']
        else:
            return None

    # 更改任务状态
    def set_job_status(self,job_id,status):
        status_map = ('wait','running','success','failed')
        if status not in status_map:
            raise ValueError(f'status must be one of {status_map}')

        if status == 'running':
            requests.post(self.job_updateStatus_url, data={'job_id': job_id, 'key': 'start_time','value': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")})
            requests.post(self.job_updateStatus_url, data={'job_id': job_id, 'key': 'status', 'value': 'running'})

        if status == 'success':
            requests.post(self.job_updateStatus_url, data={'job_id': job_id, 'key': 'end_time','value': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")})
            requests.post(self.job_updateStatus_url, data={'job_id': job_id, 'key': 'status', 'value': 'success'})

        if status == 'failed':
            requests.post(self.job_updateStatus_url, data={'job_id': job_id, 'key': 'end_time','value': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")})
            requests.post(self.job_updateStatus_url, data={'job_id': job_id, 'key': 'status', 'value': 'failed'})

        if status == 'wait':
            requests.post(self.job_updateStatus_url, data={'job_id': job_id, 'key': 'end_time','value': None})
            requests.post(self.job_updateStatus_url, data={'job_id': job_id, 'key': 'start_time', 'value': None})
            requests.post(self.job_updateStatus_url, data={'job_id': job_id, 'key': 'status', 'value': 'wait'})

if __name__ == '__main__':
    self = BioinfodataPlatformWeb()
    print(self.version())




