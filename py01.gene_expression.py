#!/usr/bin/python

from __future__ import division
import sys
sys.path.append("..")
import module.module01 as m01
import subprocess
import cPickle as pickle
import time
import string
import matplotlib.pyplot as plt
import os

db = pickle.load(open("hg19_Ref_Gen.dat","rb"))

#FUNc:integrate_expr_file(self,files,counts_colume,sample_order,g_or_t)
#FUNC:integrate_cuff_file(self,filename,sample_order)
#FUNC:genes_detect(self,record,normalize_way,threhold)
#FUNC:correlation(self,normalize,record,gene_class)
#FUNC:integrate_expr_file(self,files,gene_class,counts_colume,sample_order,g_or_t)
#FUNC:genes_detect(self,record,gene_class,normalize_way,threhold)
#FUNC:hist_of_selected_genes(self,genelist,samples,sample_groups,groupnames,normal_way,outfile)
#LISTS:['both', 'ref', 'gen_anti', 'ens_only', 'gen_lincRNA', 'gen_coding', 'ref_lincRNA', 'gen']

try:
	cells = pickle.load(open("Data.dat","rb"))
except:
	n = 35
	ids = ["SRR944285","SRR944290","SRR944293","SRR944295","SRR944298","SRR944300","SRR944302","SRR944304","SRR944309","SRR944311","SRR944313","SRR944315","SRR944316","D03","D04","D11","D12","D13","D14","D15","D16","D17","D18","D19","D03_SE","D04_SE","D11_SE","D12_SE","D13_SE","D14_SE","D15_SE","D16_SE","D17_SE","D18_SE","D19_SE"]
	names = ids
	cells = m01.expression(ids,names,n)
	cells.initialize_express_data(db,n)
	for i in range(0,n):
		cells.integrate_expr_file(['./00.data/Sum.refseq.'+names[i]+'.txt','./00.data/Sum.Gencode.'+names[i]+'.txt'],['ref','gen'],1,i,'T')
		cells.integrate_cuff_file('./00.data/Ensg.'+names[i]+'.txt',i)
	cells.correlation("RPKM","corr_RPKM_log","ref")
	cells.corr_unlog("RPKM","corr_RPKM_unlog","ref")
	cells.correlation("FPKM","corr_FPKM_log","ref")
        cells.corr_unlog("FPKM","corr_FPKM_unlog","ref")
	cells.genes_detect('RPKM_0.1',"ref",'RPKM',0.1)
	cells.genes_detect('FPKM_0.1',"ref",'FPKM',0.1)
	cells.genes_detect('RPKM1',"ref",'RPKM',1)
	cells.genes_detect('FPKM_1',"ref",'FPKM',1)
	cells.genes_detect('count2',"ref",'count',2)
	cells.genes_detect('gen_linc_RPKM_0.1',"gen_lincRNA",'RPKM',0.1)
	cells.genes_detect('gen_linc_FPKM_0.1',"gen_lincRNA",'FPKM',0.1)
	cells.genes_detect('gen_anti_RPKM_0.1',"gen_anti",'RPKM',0.1)
        cells.genes_detect('gen_anti_FPKM_0.1',"gen_anti",'FPKM',0.1)
	pickle.dump(cells,open("Data.dat","wb"),True)
	f = open('Gene_log_RPKM_FPKM.txt','w')
	for i in cells['geneinfo']:
		print >>f,i,'\t',"	".join(str(n) for n in cells['geneinfo'][i]['logRPKM'])
	f.close
	f = open('Gene_count_RPKM_FPKM.txt','w')
	for i in cells['geneinfo']:
	        print >>f,i,'\t',"\t".join(str(n) for n in cells['geneinfo'][i]['count']),"\t".join(str(n) for n in cells['geneinfo'][i]['RPKM'])
	f.close
if not os.path.exists('info.txt'):
	f = open('info.txt','w')
	print >>f,"corr_RPKM_log"
	for i in cells['corr_RPKM_log']:
	        print >>f,'\t',"\t".join(str(n) for n in cells['corr_RPKM_log'][i])
	print >>f,"corr_FPKM_log"
	for i in cells['corr_FPKM_log']:
	        print >>f,'\t',"\t".join(str(n) for n in cells['corr_FPKM_log'][i])
	print >>f,"corr_RPKM_unlog"
	for i in cells['corr_RPKM_unlog']:
	        print >>f,'\t',"\t".join(str(n) for n in cells['corr_RPKM_unlog'][i])
	print >>f,"corr_FPKM_unlog"
	for i in cells['corr_FPKM_unlog']:
	        print >>f,'\t',"\t".join(str(n) for n in cells['corr_FPKM_unlog'][i])
	print >>f,"Mapreads",",".join(str(cells['sampinfo'][n]["mapreads"]) for n in cells['sample'])
	print >>f,"RPKM_0.1",cells['gene_det']['RPKM_0.1']
	print >>f,"FPKM_0.1",cells['gene_det']['FPKM_0.1']
	print >>f,"count2",cells['gene_det']['count2']
	print >>f,"RPKM_1",cells['gene_det']['RPKM1']
	print >>f,"FPKM_1",cells['gene_det']['FPKM_1']
	print >>f,'gen_linc_RPKM_0.1',cells['gene_det']['gen_linc_RPKM_0.1']
	print >>f,'gen_linc_FPKM_0.1',cells['gene_det']['gen_linc_FPKM_0.1']
	print >>f,'gen_anti_RPKM_0.1',cells['gene_det']['gen_anti_RPKM_0.1']
	print >>f,'gen_anti_FPKM_0.1',cells['gene_det']['gen_anti_FPKM_0.1']
	f.close

#FUNC:mapreads_count_genes_number(self,normalize_way,color,legends)
normalize_way=['RPKM_0.1','RPKM1']
color=['r','b']
legends=normalize_way
cells.mapreads_count_genes_number(normalize_way,color,legends)

if not os.path.exists('MT-genes_RPKM.png'):
	samples = ["SRR944285","SRR944290","SRR944293","SRR944295","SRR944298","SRR944300","SRR944302","SRR944304","SRR944309","SRR944311","SRR944313","SRR944315","SRR944316","D03","D04","D11","D12","D14","D15","D16","D17","D18","D19"]
	sample_group = [0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1]
	MT_genes = ["MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND5","MT-ND6","MT-CYB","MT-ATP8","MT-CO1","MT-CO2","MT-ND4L","MT-ATP6","MT-CO3"]
	cells.hist_of_selected_genes(MT_genes,samples,sample_group,['smart2','new'],"FPKM","MT-genes_FPKM.png")
	cells.hist_of_selected_genes(MT_genes,samples,sample_group,['smart2','new'],"RPKM","MT-genes_RPKM.png")

#FUNC:random_select_reads_subsampling(self,sample,cut_reads)
#cells.random_select_reads_subsampling('SRR944285',[2,4,6,8])

