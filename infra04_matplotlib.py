import re,sys,os
import matplotlib
matplotlib.use('Cairo')
import matplotlib.pyplot as plt
import subprocess
import cPickle as pickle
import time
import numpy as np

def array_density(array,width,height,xmax,ymax,bin_number,record):
	plt.figure(figsize=(width,height), dpi=200)
	xc = len(array)
	yc = len(array[0])
	count = 1
	print record,xc,yc
	bins = np.linspace(0,xmax,bin_number)
        for i in range(xc):
		for j in range(yc):
			ax = plt.subplot(xc,yc,count)
			ax.hist(array[i][j],bins,label='A'+str(i+1)+'_B'+str(j+1))
			plt.xlim(0,xmax)
			plt.ylim(0,ymax)
			plt.legend()
			count += 1
	plt.savefig('f.'+record+'.png',dpi=200)
def hist_dict(data,key_name,width,height,xmax,bin_number,record):
	plt.figure(figsize=(width,height), dpi=200)
	t = []
	bins = np.linspace(0,xmax,bin_number)
	for i in data:
		t.append(data[i][key_name])
	plt.hist(t,bins)
	plt.xlim(0,xmax)
	plt.savefig('f.'+record+'.png',dpi=200)
def hist_array(datas,labels,width,height,xmax,ymax,bin_number,record):
	plt.figure(figsize=(width,height), dpi=200)
	count = 1
	bins = np.linspace(0,xmax,bin_number)
	for rec in labels:
		ax = plt.subplot(len(labels),1,count)
		ax.hist(datas[rec],bins)
		plt.xlim(0,xmax)
		plt.ylim(0,ymax)
		count += 1
		plt.title(rec)
	plt.savefig('f.'+record+'.png',dpi=200)
def density_array(datas,labels,width,height,xmax,ymax,bin_number,record):
	plt.figure(figsize=(width,height), dpi=200)
	count = 1
	bins = np.linspace(0,xmax,bin_number)
	for rec in labels:
		ax = plt.subplot(len(labels),1,count)
		ax.density(datas[rec],bins)
		plt.xlim(0,xmax)
		count += 1
		plt.title(rec)
	plt.savefig('f.'+record+'.png',dpi=200)
def scatters_sets(datas,labels,width,height,xmax,ymax,size,record):
	plt.figure(figsize=(width,height), dpi=200)
	count = 1
	for rec in labels:
		ax = plt.subplot(len(labels),1,count)
		ax.scatter(datas[rec][0],datas[rec][1],s=size)
		plt.xlim(0,xmax)
		plt.ylim(0,ymax)
		count += 1
		plt.title(rec)
	plt.savefig('f.'+record+'.png',dpi=200)
def scatters_layers(datas,labels,width,height,xmax,ymax,sizes,colors,record):
	plt.figure(figsize=(width,height), dpi=200)
	count = 0
	for data in datas:
		plt.scatter(data[0],data[1],s=sizes[count],color=colors[count],label=labels[count])
		plt.legend(loc = 2)
		count += 1
	plt.xlim(0,xmax)
	plt.ylim(0,ymax)
	plt.savefig('f.'+record+'.png',dpi=200)

