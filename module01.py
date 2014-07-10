#!/usr/bin/python

from __future__ import division
import re
import matplotlib
from matplotlib.mlab import PCA
#matplotlib.use('Cairo')
import matplotlib.pyplot as plt
#import matplotlib_venn
#from matplotlib_venn import venn3,venn3_circles,venn2,venn2_circles 
import subprocess
import cPickle as pickle
import time
import numpy as np

class expression(dict):
	def __init__(self,id_input,name_input,counts):
		self['id'] = id_input
		self['sample'] = name_input
		self['sample_count'] = counts
		self['genelists'] = []
		self['sampinfo'] = {}
		self['geneinfo'] = {}
		self['gene_det'] = {}
		self['maprate'] = {}
	def initialize_express_data(self,db,total_counts):
		self['db'] = db
		self['sample_number'] = total_counts
		self['genelist'] =  db['array_gene']
		for i in db['array_gene']:
			self['geneinfo'][i] = {}
			self['geneinfo'][i]['length'] = self['db']['genes'][i]['length']
			self['geneinfo'][i]['count'] = [0]*total_counts
			self['geneinfo'][i]['RPKM'] = [0]*total_counts
			self['geneinfo'][i]['FPKM'] = [0]*total_counts
			self['geneinfo'][i]['logRPKM'] = [0]*total_counts
			self['geneinfo'][i]['logFPKM'] = [0]*total_counts
			self['sampinfo'] = {}
			self['geneinfo'][i]["detect_number"] = 0
			self['geneinfo'][i]["detect_cond"] = ""
	def integrate_expr_file(self,files,gene_classes,counts_colume,sample_order,g_or_t):
		sample = self['sample'][sample_order]		
		for i in range(0,len(gene_classes)):
			gene_class = gene_classes[i]
			file = open(files[i])
			for line in file:
				array = str.split(line)
				if g_or_t == 'G':
					gene = array[0]
				else:
					gene = self['db'][gene_class]['t'][array[0]]
				if gene in self['geneinfo']:
					try:
						self['geneinfo'][gene]['count'][sample_order] += int(array[counts_colume])
					except:
						print array[0],array,self['geneinfo'][gene]['count'],"NO_GENE_EXPRESSION_DATA"
		self['sampinfo'][sample] = {}
		self['sampinfo'][sample]["count"] = [] 
		self['sampinfo'][sample]["RPKM"] = []
		self['sampinfo'][sample]["FPKM"] = []
		self['sampinfo'][sample]["logRPKM"] = []
		self['sampinfo'][sample]["logFPKM"] = []
		self['sampinfo'][sample]["mapreads"] = 0
		for i in self['genelist']:
			self['sampinfo'][sample]["count"].append(self['geneinfo'][i]['count'][sample_order])
			self['sampinfo'][sample]["mapreads"] += self['geneinfo'][i]['count'][sample_order]
		for i in self['genelist']:
			self['sampinfo'][sample]["RPKM"].append(float(self['geneinfo'][i]['count'][sample_order]*1000*1000000/(self['sampinfo'][sample]["mapreads"]*self['geneinfo'][i]['length'])))
			self['geneinfo'][i]['RPKM'][sample_order] = float(self['geneinfo'][i]['count'][sample_order]*1000*1000000/(self['sampinfo'][sample]["mapreads"]*self['geneinfo'][i]['length']))
			try:
				self['geneinfo'][i]["logRPKM"][sample_order] = max(np.log10(self['geneinfo'][i]['RPKM'][sample_order]),-2)
			except:
				self['geneinfo'][i]["logRPKM"][sample_order] = -2
	def integrate_cuff_file(self,filename,sample_order):
		sample = self['sample'][sample_order]
		file = open(filename)
		file.readline()
		for line in file:
			array = str.split(line)
			gene = array[4]
			if gene in self['geneinfo']:
				self['geneinfo'][gene]['FPKM'][sample_order] = float(array[9])
				try:
					self['geneinfo'][gene]['logFPKM'][sample_order] = max(np.log10(array[9]),-2)
				except:
					self['geneinfo'][gene]["logFPKM"][sample_order] = -2
		for i in self['genelist']:
			self['sampinfo'][sample]["FPKM"].append(self['geneinfo'][i]["FPKM"][sample_order])
			self['sampinfo'][sample]["logFPKM"].append(self['geneinfo'][i]["logFPKM"][sample_order])
	def genes_detect(self,record,gene_class,normalize_way,threhold):
		tStart = time.time()
		self['gene_det'][record] = []
		for i in self['sample']:
			temp = 0
			for j in self['db']['lists_ids'][gene_class]:
				if self['sampinfo'][i][normalize_way][j] >= threhold:
					temp = temp+1
			self['gene_det'][record].append(temp)
		tEnd = time.time()
		print record,
		print "It cost %f sec" % (tEnd - tStart)
	def genes_detect_pair_wise(self,samplelists,listnames,geneclass,normal,threhold,required_numbers):
		if 'venn' not in self:
			self['venn'] = {}
		indexes = []
		for i in range(0,len(samplelists)):
			index = []
			for j in samplelists[i]:
				index.append(self['sample'].index(j))
			indexes.append(index)
		for i in indexes:
			id = indexes.index(i)
			self['venn'][listnames[id]] = []
			for j in self['db']['lists'][geneclass]:
				number = 0
				for k in i:
					if self['geneinfo'][j][normal][k] >= threhold:
						number += 1
				if number >= required_numbers[id]:
					self['venn'][listnames[id]].append(j)
	def sample_bin_genes_with_bulk_expression(self,gene_class,normal_way,record,legend_name,samples,colors):
                plt.figure(figsize=(4,4), dpi=300)
		count = 0
		for lists in samples:
			rec_std = record+'_std'
	                index = []
			bulk_index = self['sample'].index('G01')
	       	        for i in lists:
	                        index.append(self['sample'].index(i))
	               	self['bin'] = {}
	                self['bin'][record] = [0]*20
	                self['bin'][rec_std] = [0]*20
	                temp = [0]*20
	                xlabel = []
	                for i in range(0,20):
	                        temp[i] = []
	                        xlabel.append(0.2*i)
	                for i in self['db']['lists'][gene_class]:
	                        t = [self['geneinfo'][i][normal_way][j] for j in index]
				bulk = self['geneinfo'][i][normal_way][bulk_index]
				if bulk >= 0:
	                                bin_index = min(int(bulk/0.2),19)
	                                std = np.std(t,ddof=1)
	                                temp[bin_index].append(std)
	                for i in range(0,20):
	                        print len(temp[i])
	                        self['bin'][record][i] = np.average(temp[i])
	                        self['bin'][rec_std][i] = np.std(temp[i],ddof=1)
	                plt.plot(xlabel,self['bin'][record],linewidth=2,linestyle="-",color=colors[count],label=legend_name[count])
			count += 1
	        plt.legend(loc=1)
                plt.savefig('f.'+record+'.png')
                plt.clf()
	def random_selected_samples(self,inlist,counts,rec):
		while 1:
			out = []
			for i in inlist:
				if np.random.uniform(0,1)<float(counts/len(inlist)):
					out.append(i)
			if len(out) == 6:
				break
		self[rec] = out
	def gene_detection_ability_with_RPKM_bins(self,gene_class,normal_way,threhold_up,threhold_down,fig_name,legend_name,samples,colors):
		f = open('info.sample_RPKM_detected.txt','w')
		RPKM_bin_genes = [0]*20
		bin_width = threhold_up/20
		xlabel = []
		bulk_index = self['sample'].index('G01')
		for i in range(0,20):
			RPKM_bin_genes[i] = []
			xlabel.append(bin_width*i)
		for i in self['db']['lists'][gene_class]:
			bulk = self['geneinfo'][i][normal_way][bulk_index]
			if bulk >= 0:
				bin_index = min(int(bulk/bin_width),19)
				RPKM_bin_genes[bin_index].append(i)
		plt.figure(figsize=(10,10), dpi=300)
		count = 0 
		for sample in samples:
			if not isinstance(sample,list):
				print >>f,sample,
				detected_number = [0]*20
				detected_ratio = [0]*20
				print sample,self['sample']		
				index = self['sample'].index(sample)
				for i in range(0,20):
					for gene in RPKM_bin_genes[i]:
						if self['geneinfo'][gene][normal_way][index] >= threhold_down:
							detected_number[i] += 1
					print sample,i,len(RPKM_bin_genes[i]),detected_number[i]
					detected_ratio[i] = detected_number[i]/len(RPKM_bin_genes[i])
					print >>f,detected_ratio[i],
                       		plt.plot(xlabel,detected_ratio,linewidth=2,linestyle="-",color=colors[count],label=legend_name[count])
				print >>f,'##'
			else:
                                detected_ratios = [0]*len(sample)
                                detected_ratio_mean = [0]*20
				detected_ratio_std = [0]*20
				for samp in sample:
					index = self['sample'].index(samp)
					id = sample.index(samp)
					detected_number = [0]*20
					detected_ratio = [0]*20
					for i in range(0,20):
						for gene in RPKM_bin_genes[i]:
							if self['geneinfo'][gene][normal_way][index] >= threhold_down:
								detected_number[i] += 1
						detected_ratio[i] = detected_number[i]/len(RPKM_bin_genes[i])
					detected_ratios[id] = detected_ratio
				for i in range(0,20):
					temp = [detected_ratios[k][i] for k in range(0,len(sample))]
					detected_ratio_mean[i] = np.average(temp)
					detected_ratio_std[i] = np.std(temp,ddof=1)
				plt.errorbar(xlabel,detected_ratio_mean,yerr=detected_ratio_std,linewidth=2,linestyle="-",color=colors[count],label=legend_name[count])
				print >>f,xlabel,detected_ratio_mean
				print >>f,xlabel,detected_ratio_std	
				count += 1
		plt.axis([0,3.6,0,1])
		plt.legend(loc=4)
                plt.savefig('f.png.'+fig_name+'.png',dpi=300)
		plt.savefig('f.eps.'+fig_name+'.eps')
                plt.clf()
	def geneinfo_item_generator(self,process,samples,p):
		if process == 'identical':
			self[p[1]] = []
			for i in self['geneinfo']:
				self['geneinfo'][i][p[1]] = self['geneinfo'][i][p[0]][self['sample'].index(samples[0])]
				self[p[1]].append(self['geneinfo'][i][p[1]])
		if process == 'count':
			for i in self['geneinfo']:
				count = 0
				for j in samples:
					if self['geneinfo'][i][p[0]][self['sample'].index(j)] >= p[1]:
						count += 1
				self['geneinfo'][i][p[2]] = float(count)/len(samples)
		if process == 'stat':
			for i in self['geneinfo']:
				RPKM = []
				for j in samples:
					RPKM.append(self['geneinfo'][i][p[0]][self['sample'].index(j)])
				self['geneinfo'][i][p[1]] = np.average(RPKM)
				if np.average(RPKM)>0:
					self['geneinfo'][i][p[2]] = np.std(RPKM,ddof=1)/np.average(RPKM)
				else:
					self['geneinfo'][i][p[2]] = 0
		if process == 'log':
			for i in self['geneinfo']:
				self['geneinfo'][i][p[1]] = np.log10(self['geneinfo'][i][p[0]]+0.01)
	def sampinfo_generator(self,**kwargs):
		if kwargs['meth'] == 'sigma':
			self['sampinfo'][kwargs['rec']] = {}
			self['sampinfo'][kwargs['outliers']] = {}
			self['sampinfo'][kwargs['outliers']]['RPKM_mean_log'] = []
			self['sampinfo'][kwargs['outliers']]['RPKM_std_norm'] = []
			self['sampinfo'][kwargs['outliers']]['gene_lists'] = []
#			self['sampinfo'][kwargs['outliers']]['anno'] = {}
			x = np.array(self['sampinfo'][kwargs['sampleid']]['RPKM_mean_log'])
			y = np.array(self['sampinfo'][kwargs['sampleid']]['RPKM_std_norm'])
			genes = np.array(self['sampinfo'][kwargs['sampleid']]['gene_lists'])
			x_target = np.array(self['sampinfo'][kwargs['targetid']]['RPKM_mean_log'])
                        y_target = np.array(self['sampinfo'][kwargs['targetid']]['RPKM_std_norm'])
                        genes_target = np.array(self['sampinfo'][kwargs['targetid']]['gene_lists'])
			xdata = []
			ydata = []
			for i in np.arange(kwargs['range'][0],kwargs['range'][1],kwargs['step']):
				left = i-kwargs['flank']
				right = i+kwargs['flank']
				print left,right,np.average(y[np.logical_and(x<right,x>=left)])
                                up = np.average(y[np.logical_and(x<right,x>=left)])+kwargs['std']*np.std(y[np.logical_and(x<right,x>=left)])
				xdata.append(i)
				ydata.append(up)
				gene_in_range = genes_target[np.logical_and(x_target<i+0.5*kwargs['step'],x_target>=i-0.5*kwargs['step'])]			
				mean_in_range = x_target[np.logical_and(x_target<i+0.5*kwargs['step'],x_target>=i-0.5*kwargs['step'])]
				cv_in_range = y_target[np.logical_and(x_target<i+0.5*kwargs['step'],x_target>=i-0.5*kwargs['step'])]
				for j,gene in enumerate(gene_in_range):
					if cv_in_range[j]>up and gene not in self['sampinfo'][kwargs['outliers']]['gene_lists']:
						self['sampinfo'][kwargs['outliers']]['RPKM_mean_log'].append(mean_in_range[j])
						self['sampinfo'][kwargs['outliers']]['RPKM_std_norm'].append(cv_in_range[j])
						self['sampinfo'][kwargs['outliers']]['gene_lists'].append(gene)
			self['sampinfo'][kwargs['rec']]['RPKM_mean_log'] = xdata
			self['sampinfo'][kwargs['rec']]['RPKM_std_norm'] = ydata
			if 'file' in kwargs:
				thefile = open(kwargs['file'],'w')
				for i,gene in enumerate(self['sampinfo'][kwargs['outliers']]['gene_lists']):
					print >>thefile,gene,self['sampinfo'][kwargs['outliers']]['RPKM_mean_log'][i],self['sampinfo'][kwargs['outliers']]['RPKM_std_norm'][i]
				thefile.close()
		if kwargs['meth'] == 'subsample':
			self['sampinfo'][kwargs['rec']] = {}
			self['sampinfo'][kwargs['rec']]['RPKM_mean_log'] = []
                        self['sampinfo'][kwargs['rec']]['RPKM_std_norm'] = []
                        self['sampinfo'][kwargs['rec']]['gene_lists'] = kwargs['genes']
			for gene in kwargs['genes']:
				if gene in self['sampinfo'][kwargs['sampleid']]['gene_lists']:
					index = self['sampinfo'][kwargs['sampleid']]['gene_lists'].index(gene)
					self['sampinfo'][kwargs['rec']]['RPKM_mean_log'].append(self['sampinfo'][kwargs['sampleid']]['RPKM_mean_log'][index])
					self['sampinfo'][kwargs['rec']]['RPKM_std_norm'].append(self['sampinfo'][kwargs['sampleid']]['RPKM_std_norm'][index])
				else:
					print gene,"NOt_in_lists"
	def geneinfo_plot_scatter(self,geneclass,cat1,cat2,sizes,colors,config,figname):
		xdata = []
		ydata = []
		for i in self['db']['lists'][geneclass]:
			a = self['geneinfo'][i][cat1]
			b = self['geneinfo'][i][cat2]
			if 'xlog' in config and self['geneinfo'][i][cat1]<=config[2]:
				a = config[2]
			if 'ylog' in config and self['geneinfo'][i][cat2]<=config[4]:
				b = config[4]
			if 'CV' in config and self['geneinfo'][i][cat1] == 0:
				a = config[3]
			if 'CV' in config and self['geneinfo'][i][cat2] == 0:
                                b = config[5]
			xdata.append(a)
			ydata.append(b)
		if re.findall('map',colors):
                        plt.figure(figsize=(10,8), dpi=200)
                        cat = re.split(':',colors)[1]
                        dotcolors = []
                        for i in self['db']['lists'][geneclass]:
                                dotcolors.append(self['geneinfo'][i][cat])
		else:
                        plt.figure(figsize=(8,8), dpi=200)
                        dotcolors = [colors]*len(xdata)
		if 'SORT' in config:
			xdata = np.array(xdata)
			ydata = np.array(ydata)	
			orders = xdata.argsort()
			xdata = xdata[orders]
			ydata = ydata[orders]
			dotcolors = np.array(dotcolors)
			dotcolors = dotcolors[orders]
			print xdata,ydata
			plt.scatter(range(0,len(xdata)),ydata,c=dotcolors,edgecolors='none',s=sizes)
			plt.colorbar()
			plt.scatter(range(0,len(xdata)),xdata,c='black',edgecolors='none',s=sizes)
		else:
			plt.scatter(xdata,ydata,c=dotcolors,edgecolors='none',s=sizes)
		if 'xlog' in config:
                	plt.xscale('log')
               	if 'ylog' in config:
                	plt.yscale('log')
		plt.xlabel(cat1,fontsize=18)
		plt.ylabel(cat2,fontsize=18)
		plt.xlim(config[2],config[3])
		plt.ylim(config[4],config[5])
		if re.findall('map',colors):
			plt.colorbar()
		plt.savefig('./f.png.'+figname+'.png',dpi=400)
		plt.savefig('./f.eps.'+figname+'.eps')
		plt.clf()
	def geneinfo_bin_line_plot(self,geneclass,numbers,cat1,cat2,xmin,xmax,DIS_CONTINOUS,binsize,colors,legends,xlabel,ylabel,figname):
		plt.figure(figsize=(8,8), dpi=200)
		for n in range(numbers):
	                xdata = []
	                ydata = []
	                for i in self['db']['lists'][geneclass]:
	                        xdata.append(self['geneinfo'][i][cat1[n]])
	                        ydata.append(self['geneinfo'][i][cat2[n]])
			x = np.array(xdata)
			y = np.array(ydata)
			ymean = []
			ystd = []
			if DIS_CONTINOUS == 'DIS':
				xpos = range(xmin,xmax+1)
				xticks = [str(i) for i in range(xmin,xmax+1)]
				for i in range(xmin,xmax+1):
					ymean.append(np.average(y[x==i]))
					ystd.append(np.std(y[x==i],ddof=1))
			if DIS_CONTINOUS == 'CON':
				ranges = np.arange(xmin,xmax,binsize)
				xpos = []
				xticks = []
				for i in range(1,len(ranges)):
					left = ranges[i-1]
					right = ranges[i]
					ymean.append(np.average(y[np.logical_and(x<right,x>=left)]))
					ystd.append(np.std(y[np.logical_and(x<right,x>=left)]))
					xpos.append(left)
					xticks.append(str(round(left,2)))
			plt.errorbar(xpos,ymean,yerr=ystd,linewidth=2,linestyle="-",color=colors[n],label=legends[n])
		plt.xticks(xpos,xticks,rotation=45, fontsize=10)
		plt.legend(loc=2)
		plt.xlabel(xlabel,fontsize=20)
		plt.ylabel(ylabel,fontsize=20)
		plt.savefig('f.png.'+figname+'.png',dpi=400)
		plt.savefig('f.eps.'+figname+'.eps')
	def genes_variation_of_samples(self,samples,geneclass,normal_way,catname):
		samplesize = len(samples)
		if 'samplebin' not in self:
			self['samplebin'] = []
		if 1==1:
			self['samplebin'].append(catname)
			self['sampinfo'][catname] = {}
			self['sampinfo'][catname]['RPKM'] = []
			self['sampinfo'][catname]['logRPKM'] = []
			self['sampinfo'][catname]['gene_lists'] = []
			self['sampinfo'][catname]['RPKM_mean'] = []
			self['sampinfo'][catname]['RPKM_mean_log'] = []
			self['sampinfo'][catname]['logRPKM_mean'] = []
			self['sampinfo'][catname]['RPKM_std'] = []
			self['sampinfo'][catname]['RPKM_std_norm'] = []
			self['sampinfo'][catname]['RPKM_std_norm_samplesize'] = []
			self['sampinfo'][catname]['logRPKM_std'] = []
			self['sampinfo'][catname]['logRPKM_std_norm'] = []
		for gene in self['db']['lists'][geneclass]:
			RPKM = []
			logRPKM = []
			for i in samples:
				RPKM.append(self['geneinfo'][gene]['RPKM'][self['sample'].index(i)])
				logRPKM.append(self['geneinfo'][gene]['logRPKM'][self['sample'].index(i)])
			if np.average(RPKM)==0 or np.average(logRPKM)==0:
				continue
			self['sampinfo'][catname]['RPKM'].append(RPKM)
			self['sampinfo'][catname]['logRPKM'].append(logRPKM)
			self['sampinfo'][catname]['gene_lists'].append(gene)
			self['sampinfo'][catname]['RPKM_mean'].append(np.average(RPKM))
			self['sampinfo'][catname]['logRPKM_mean'].append(np.average(logRPKM))
			self['sampinfo'][catname]['RPKM_mean_log'].append(np.log10(np.average(RPKM)))
			self['sampinfo'][catname]['RPKM_std'].append(np.std(RPKM,ddof=1))
			self['sampinfo'][catname]['logRPKM_std'].append(np.std(logRPKM,ddof=1))
			self['sampinfo'][catname]['RPKM_std_norm'].append(np.std(RPKM,ddof=1)/np.average(RPKM))
			self['sampinfo'][catname]['RPKM_std_norm_samplesize'].append(np.std(RPKM,ddof=1)/(np.average(RPKM)*np.sqrt(samplesize)))
			self['sampinfo'][catname]['logRPKM_std_norm'].append(np.std(logRPKM,ddof=1)/np.average(logRPKM))
	def correlation(self,normalize,threhold,record,gene_class):
		self[record] = [[0]*len(self['sample']) for i in range(len(self['sample']))]
		for i in range(0,len(self['sample'])):
			print "CORR_CALCULATE",self['sample'][i],i
			self[record][i][i] = 1
			for j in range(i+1,len(self['sample'])):
				t1 = []
				t2 = []
				for k in self['db']['lists'][gene_class]:
					t1.append(self['geneinfo'][k][normalize][i])
					t2.append(self['geneinfo'][k][normalize][j])
				t1 = np.array(t1)	
				t2 = np.array(t2)
				log_t1 = np.log10(t1+0.001)
				log_t2 = np.log10(t2+0.001)
				log_t1[log_t1<threhold] = threhold
				log_t2[log_t2<threhold] = threhold
				corr = np.corrcoef(log_t1[(log_t1>threhold) | (log_t2>threhold)],log_t2[(log_t1>threhold) | (log_t2>threhold)])[1][0]
				self[record][i][j] = corr
				self[record][j][i] = corr
	def corr_unlog(self,normalize,threhold,record,gene_class):
                self[record] = [[0]*len(self['sample']) for i in range(len(self['sample']))]
                for i in range(0,len(self['sample'])):
                        print "CORR_CALCULATE_unlog",self['sample'][i],i
                        self[record][i][i] = 1
                        for j in range(i+1,len(self['sample'])):
				t1 = []
                                t2 = []
                                for k in self['db']['lists'][gene_class]:
                                        t1.append(self['geneinfo'][k][normalize][i])
                                        t2.append(self['geneinfo'][k][normalize][j])
                                t1 = np.array(t1)                            
                                t2 = np.array(t2)
                                corr = np.corrcoef(t1[(t1>threhold) | (t2>threhold)],t2[(t1>threhold) | (t2>threhold)])[1][0]
				self[record][i][j] = corr
                                self[record][j][i] = corr
	def corr_on_demand(self,genelists,threhold,cat1,cat2,rec):
		t1 = []
		t2 = []
		for k in genelists:
			t1.append(self['geneinfo'][k][cat1])
			t2.append(self['geneinfo'][k][cat2])
			print k,self['geneinfo'][k][cat1],self['geneinfo'][k][cat2]
		t1 = np.array(t1)
		t2 = np.array(t2)
		corr = np.corrcoef(t1[(t1>threhold) | (t2>threhold)],t2[(t1>threhold) | (t2>threhold)])[1][0]
		self[rec] = corr
		print corr
	def hist_of_selected_genes(self,genelist,samples,sample_groups,groupnames,normal_way,outfile):
                N = len(genelist)
                groupinfo = {}
                for i in range(0,len(groupnames)):
                        group_sample_orders = []
                        groupinfo[i] = {}
                        groupinfo[i]['mean'] = []
                        groupinfo[i]['std'] = []
                        for j in range(sample_groups.index(i),sample_groups.index(i)+sample_groups.count(i)):
                                group_sample_orders.append(self['sample'].index(samples[j]))
                        for gene in genelist:
                                genes_exp = []
                                for sample in group_sample_orders:
                                        genes_exp.append(self['geneinfo'][gene][normal_way][sample])
                                groupinfo[i]['mean'].append(np.average(genes_exp))
                                groupinfo[i]['std'].append(np.std(genes_exp,ddof=1))
                print groupinfo
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ind = np.arange(N)
                width = 0.25
                bar_c = ['blue','red','green']
                error_c = ['black','black','black']
                for i in range(0,len(groupnames)):
                        ax.bar(ind+width*i,groupinfo[i]['mean'], width,color=bar_c[i],yerr=groupinfo[i]['std'],error_kw=dict(elinewidth=2,ecolor=error_c[i]),label=groupnames[i])
                ax.set_xlim(-width,len(ind)+width)
                ax.set_ylabel('RPKM')
                ax.set_title('MT-genes')
                xTickMarks = genelist
                ax.set_xticks(ind+width)
                xtickNames = ax.set_xticklabels(xTickMarks)
                plt.setp(xtickNames, rotation=45, fontsize=10)
                ax.legend(loc=2)
                plt.savefig(outfile)
                plt.clf()
	def hist_of_gene_for_samples(self,samples,gene,normal_way,config,figname):
		express = np.array([self['geneinfo'][gene][normal_way][self['sample'].index(k)] for k in samples])
		expression = express
		ind = np.arange(len(samples))
		names = np.array(samples)
		if 'SORT' in config:
			expression = express[express.argsort()]
			names = names[express.argsort()]
		width = 0.4
		plt.figure()
		plt.bar(ind,expression,0.4)
		plt.ylabel(normal_way)
		plt.xticks(ind+width,names,rotation=45, fontsize=10)
		plt.savefig('f.png.'+figname+'.png')
		plt.savefig('f.eps.'+figname+'.eps')
                plt.clf()
		plt.hist(express)
		plt.savefig('f.png.hist'+figname+'.png')
                plt.savefig('f.eps.hist'+figname+'.eps')
	def venn_digram(self,names,figname):
		import matplotlib_venn
		from matplotlib_venn import venn3,venn3_circles,venn2,venn2_circles
		plt.figure(figsize=(4,4))
		if len(names)==2:
			set1 = set(self['venn'][names[0]])
			set2 = set(self['venn'][names[1]])
			venn2([set1, set2], (names[0],names[1]))
			venn2_circles([set1, set2])
		if len(names)==3:
                        set1 = set(self['venn'][names[0]])
                        set2 = set(self['venn'][names[1]])
			set3 = set(self['venn'][names[2]])
			venn3([set1, set2, set3], (names[0],names[1],names[2]))
			venn3_circles([set1, set2, set3])
		plt.savefig('f.png.'+figname+'.png')
		plt.savefig('f.eps.'+figname+'.eps')
                plt.clf()
	def heatmap_corr(self,cat,size,figname):
		plt.figure(figsize=(size*1.1,size))
		heat = []
		for i in range(0,self['sample_count']):
                        heat.append(self[cat][i])
		plt.imshow(heat,cmap=plt.cm.jet,interpolation='nearest')
		plt.xticks(range(0,len(self['sample'])),self['sample'],rotation=45,fontsize=9)
		plt.colorbar()
		plt.savefig('f.eps.'+figname+'.eps')
		plt.savefig('f.png.'+figname+'.png')
		plt.clf()
	def writing_gene_expression_files(self,geneclass,normal_way,filenames):
		f = open("info.genes."+filenames+".txt",'w')
		print >>f,'genes','\t',"\t".join(n for n in self['sample'])
		for i in self['db']['lists'][geneclass]:
        	        print >>f,i,'\t',"\t".join(str(n) for n in self['geneinfo'][i][normal_way])
	        f.close
	def writing_gene_info(self,genes,normal_way,samples,file):
		f = open(file,'w')
	#	print >>f,'genes','\t',"\t".join(n for n in samples)
		index = []
		for i in samples:
			index.append(self['sample'].index(i))
		for i in genes:
			print >>f,i,'\t',"\t".join(str(self['geneinfo'][i][normal_way][n]) for n in index)
	def scat_genes_express(self,threhold,pairs_numbers,cat1,cat2,colors,figname):
                plt.figure(figsize=(8,8), dpi=200)
		for i in range(0,pairs_numbers):
			plt.scatter(cat1[i],cat2[i],s=1.5,edgecolors='none',facecolors=colors[i])
		plt.axis([threhold,4,threhold,4])
                plt.savefig('./f.png.corr.'+figname+'.png')
     		plt.savefig('./f.eps.corr.'+figname+'.eps')
	        plt.clf()
	def scat_samples_infos(self,samples,cat1,cat2,colors,sizes,alphas,legends,ymax,figname):
		plt.figure(figsize=(8,8), dpi=200)
		for sample in samples:
			if 'anno' not in self['sampinfo'][sample]:
				i = samples.index(sample)
				plt.scatter(self['sampinfo'][sample][cat1],self['sampinfo'][sample][cat2],s=sizes[i],edgecolors='none',facecolors=colors[i],alpha=alphas[i],label=legends[i])
			else:
				i = samples.index(sample)
                                plt.scatter(self['sampinfo'][sample][cat1],self['sampinfo'][sample][cat2],s=sizes[i],edgecolors='none',facecolors=colors[i],alpha=alphas[i],label=legends[i])
                                for j,txt in enumerate(self['sampinfo'][sample]['gene_lists']):
                                        print txt,self['sampinfo'][sample]['RPKM_mean_log'][j],self['sampinfo'][sample]['RPKM_std_norm'][j]
                                        plt.annotate(txt,(self['sampinfo'][sample]['RPKM_mean_log'][j],self['sampinfo'][sample]['RPKM_std_norm'][j]),fontsize=2)
		plt.axis([0,5,0,ymax])
		plt.legend(loc=1,fontsize=25)
		plt.xlabel('log10(average RPKM)',fontsize=18)
		plt.ylabel('CV normalized standard deviation',fontsize=18)
                plt.savefig('./f.png.'+figname+'.png')
		plt.savefig('./f.eps.'+figname+'.eps')
		plt.clf()
	def reads_mapped_to_ERCC(self,cat,rec):
		self[rec] = []
		for index,sample in enumerate(self['sample']):
			count = 0
			for gene in self['db']['lists'][cat]:
				count += self['geneinfo'][gene]['count'][index]
			self[rec].append(count)
	def pca_(self,samples,normal_way,rec):
		index = []
		genes = []
		data = []
		self[rec] = {}
		for i in samples:
			index.append(self['sample'].index(i))
		for i in self['db']['lists']['ref']:
			express = []
			for j in index:
				express.append(self['geneinfo'][i][normal_way][j])
			if np.sum(express) >0:	
				genes.append(i)
				data.append(express)
		self[rec]['data'] = data
		data = np.transpose(np.array(data))
		results = PCA(data)
		self[rec]['wt'] = results.Wt
		self[rec]['fracs'] = results.fracs
		self[rec]['Y'] = results.Y
		
			
 
			
