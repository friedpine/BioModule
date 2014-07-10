import cPickle as pickle
import re,os,time
import subprocess

class geodata(dict):
	def __init__(self):
		self['sample'] = []
	def download_geo(self,**kwargs):
		soft = open(kwargs['soft'])
		count = 0
		cwd = os.getcwd()
		while 1:
			line = soft.readline()
			if not line:
				break
			if re.match('\^SAMPLE',line):
				GSM = re.split('\s',line)[2]
			if re.match('\!Sample_title',line):
				title = re.split('\s',line)[2]
			if re.match('\!Sample_supplementary_file',line) and re.findall('SRX',line):
				link = re.split('\s',line)[2]
				name = kwargs['names'][count]
				self['sample'].append(name)
				self[name] = {}
				self[name]['link'] = link
				self[name]['title'] = title
				self[name]['GSM'] = GSM
				self[name]['fold'] = cwd+'/'+re.split('//',link)[1]
				self[name]['files'] = []
				if os.path.exists(self[name]['fold']):
					for i in os.listdir(self[name]['fold']):
						print self[name]['fold']+'/'+i
						files = os.listdir(self[name]['fold']+'/'+i)
						self[name]['files'].append(self[name]['fold']+'/'+i+'/'+files[0])
				print name,self[name]
				count += 1
		pickle.dump(self,open(kwargs['data_file'],'w'))
	def download_run(self,**kwargs):
		for i in kwargs['ranges']:
			if self[i]['files'] == []:
			        a = "wget -r "+self[i]['link']+" ./"
				print a
				child = subprocess.call(a,shell='True',stdout=subprocess.PIPE)
				if os.path.exists(self[name]['fold']):
                                        for i in os.listdir(self[name]['fold']):
                                                print self[name]['fold']+'/'+i
                                                files = os.listdir(self[name]['fold']+'/'+i)
                                                self[name]['files'].append(self[name]['fold']+'/'+i+'/'+files[0])
			pickle.dump(self,open(kwargs['data_file'],'w'))
