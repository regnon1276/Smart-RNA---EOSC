#!/usr/bin/env python3

from __future__ import print_function
import networkx as nx
import pickle as pkl
import json
import sys
import random
from networkx.readwrite import json_graph
import multiprocessing as mp
from threading import Lock
import pickle as pkl

#LIBRARIES FOR GA
import sklearn.svm

#OTHER LIB
import binascii
import socket
import numpy as np
from numpngw import write_png
import time 
import select
import os
import struct
from pathlib import Path
import atexit
import signal
import glob
from time import sleep
import time
import matplotlib.pyplot as plt
import matplotlib
import gzip
import gensim
import logging
from tkinter import *
from tkinter import ttk
import pandas as pd
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import askopenfilenames
import tkinter.simpledialog
import tkinter.messagebox
#from tkhtmlview import HTMLLabel
import webbrowser
import Bio
from Bio import SeqIO
import Bio; print(Bio.__version__)
from Bio import Align
from Bio import SeqIO
from Bio import AlignIO
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Align import AlignInfo
from Bio import SubsMat
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

url = 'https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Severe%20acute%20respiratory%20syndrome%20coronavirus%202,%20taxid:2697049'
def OpenUrl():
	webbrowser.open(url)

import itertools

from scipy import linalg

from sklearn import mixture



matplotlib.use("TkAgg")
import string
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.widgets  import RectangleSelector
import warnings

from datetime import datetime
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import style

style.use('ggplot')

sys.setrecursionlimit(10000)
warnings.filterwarnings("ignore", category=Warning)
HOST = '192.168.1.10'  # The server's hostname or IP address
PORT = 7	       # The port used by the server

ans = True

soc= socket.socket(socket.AF_INET, socket.SOCK_STREAM) 	
soc.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)	
#soc.connect((HOST, PORT))

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
mypath = os.path.join(project_root, 'Recorder_acq')
root = Tk(  )
mypath = os.path.join(project_root, '')
knowledge_loaded=0
knowledge_slice=0
files_training=''
files_testing=''
files_verif=''
np.set_printoptions(threshold=sys.maxsize)
knowledge_loaded_from_server=0
train_GENSIM_REF_file=''
test_GENSIM_TEST_file=''
counts_test =dict()
counts_ref=dict()
prot_sim_result=dict()
np_counts_ref=np.zeros((4,2))
np_counts_test=np.zeros((4,2))
DICT_REF=np.zeros((15000, 500), dtype=object )	
DICT_TEST=np.zeros((15000, 500), dtype=object )	
NB_WORD_TEST=0
NB_WORD_REF=0
line_item=0
word_ref=np.zeros((20000), dtype=object )
word_test=np.zeros((20000), dtype=object )
sentence_cleaned_ref=''
sentence_cleaned_test=''
ANNOTATION_MATRIX=np.zeros((3000, 270), dtype=np.uint8) 
category=0
RNA_AI_RECO=0

logging.basicConfig(
    format='%(asctime)s : %(levelname)s : %(message)s',
    level=logging.INFO)

###################################################
# 		VARIABLE NN LEARNING		  #
###################################################
neuron_max=512 #512
feature_number=256 #256
knowledge_from_simu=np.zeros((9+neuron_max*261), dtype=np.dtype(object))	


feature_line=np.zeros((3000,270))
def show_file_contents(input_file):
    with gzip.open(input_file, 'rb') as f:
        for i, line in enumerate(f):
            print(line)
            break

##############################################################
#		A MUTATION ANALYSIS 			     #
#							     #	
#							     #	
# 							     #	
#							     #	
#						             #
##############################################################

def read_input(input_file):
    """This method reads the input file which is in gzip format"""
    logging.info("reading file {0}...this may take a while".format(input_file))
    with gzip.open(input_file, 'rb') as f:
        for i, line in enumerate(f):

            if (i % 10000 == 0):
                logging.info("read {0} reviews".format(i))
            # do some pre-processing and return list of words for each review
            # text
            yield gensim.utils.simple_preprocess(line)
##############################################################
#		OPEN AND TRANSFORM REFERENCE RNA FILE	     #
##############################################################
def Open_REF_RNA_FILE():
	global train_GENSIM_REF_file
	global counts_ref
	global DICT_REF
	global NB_WORD_REF
	global word_ref		
	global sentence_cleaned_ref
	size_word=3
	index_word=0
#Open multiple files

	files = askopenfilenames(parent=root,initialdir=mypath+'Smart_RNA/dataset/REFERENCE DATA'
	                           ,title='Choose files',multiple = 1)
	files_for_Word2Vect = root.tk.splitlist(files)
	feature_file_for_Word2Vect_number_of_raw=0	
	tot_raw=0
	raw=np.zeros(len(files_for_Word2Vect))
	#print(len(files_for_Word2Vect))
	data_final=''
	for i in range (len(files_for_Word2Vect)):
	
		

#OPEN THE FILE DOWNLOADED FROM THE DATA BAS
		with open (files_for_Word2Vect[i], "r") as myfile:
		    data=myfile.read()
		data_raw= re.sub(r'(.{1})(?!$)','\\1 ', data[96:])

#SPLIT THE RNA IN WORDS OF X LETTERS			
		if (size_word==3):		
			data= re.sub(r'(.{3})(?!$)','\\1 ', data[96:])	
		if (size_word==4):		
			data= re.sub(r'(.{4})(?!$)','\\1 ', data[96:])	
		if (size_word==5):		
			data= re.sub(r'(.{5})(?!$)','\\1 ', data[96:])	
		if (size_word==6):		
			data= re.sub(r'(.{6})(?!$)','\\1 ', data[96:])	
		if (size_word==7):		
			data= re.sub(r'(.{7})(?!$)','\\1 ', data[96:])	
		if (size_word==8):		
			data= re.sub(r'(.{8})(?!$)','\\1 ', data[96:])	
		if (size_word==9):		
			data= re.sub(r'(.{9})(?!$)','\\1 ', data[96:])	
		if (size_word==10):		
			data= re.sub(r'(.{10})(?!$)','\\1 ', data[96:])	


#CHANGE THE TYPE TO LOWER CASE
		data=data.lower()
		NB_WORD_REF=int(len(data[96:])/size_word)
#IF MORE THAN ONE FILE SELECTED CONCAT THE FILES

		if len(files_for_Word2Vect)>1:		
			data_final=data+data_final 	
		if len(files_for_Word2Vect)==1:		
			data_final=data 	
#SAVE THE FILE TO A GZIP FORMAT COMPATIBLE WITH GENSIM
	if len(files_for_Word2Vect)==1:		
		train_GENSIM_REF_file=files_for_Word2Vect[i] +'transformed.gz'		
	if len(files_for_Word2Vect)>1:		
		train_GENSIM_REF_file=mypath+'/Smart_RNA/dataset/REFERENCE DATA/CONCAT_FILE_transformed.gz'			
	with gzip.open(train_GENSIM_REF_file, 'w') as f:
	    	f.write(data_final.encode())
	#CALCULATE STATISTICS
	word_count_ref(data_raw)
	sentence=data
	words = sentence.split(' ')
	Word_OK=0
	u=0
	sentence_cleaned=''
	sentence_cleaned_ref=''
	for k,w in enumerate(words):
		Word_OK=0				
		w=w.replace("\r","")
		w=w.replace("\n","")
		for j in range (len(w)):			
			if w[j]=='a' or w[j]=='t' or w[j]=='u' or w[j]=='g' or w[j]=='c':
				Word_OK=Word_OK+1 			
		if Word_OK==size_word:			
			DICT_REF[u+1,0]=w
			sentence_cleaned_ref=sentence_cleaned_ref+w
			sentence_cleaned=sentence_cleaned+' '+w
			u=u+1		
			word=w			
			word_ref[index_word+1]=w
			index_word=index_word+1			
	words_cleaned = sentence_cleaned.split(' ')
	
	for k,w in enumerate(words_cleaned):
		word=w		
		positions = [ i for i,w in enumerate(words_cleaned) if w == word ]		
		DICT_REF[k,1:len(positions)+1]=positions

##############################################################
#		OPEN AND TRANSOFRM TEST RNA FILE	     #
##############################################################
def Open_TEST_RNA_FILE():
	global counts_test
	global DICT_TEST
	global train_GENSIM_TEST_file
	global NB_WORD_TEST
	global word_test
	global sentence_cleaned_test
	size_word=3
	index_word=0
	files = askopenfilenames(parent=root,initialdir=mypath+'/Smart_RNA/dataset/TEST DATA'
	                           ,title='Choose files',multiple = 1)
	files_for_Word2Vect = root.tk.splitlist(files)
	print(files_for_Word2Vect)
	feature_file_for_Word2Vect_number_of_raw=0	
	tot_raw=0
	raw=np.zeros(len(files_for_Word2Vect))
	#print(len(files_for_Word2Vect))
	data_final=''
	for i in range (len(files_for_Word2Vect)):
#OPEN THE FILE DOWNLOADED FROM THE DATA BASE
		with open (files_for_Word2Vect[i], "r") as myfile:
		    data=myfile.read()
		data_raw=re.sub(r'(.{1})(?!$)','\\1 ', data[96:])
#SPLIT THE RNA IN WORDS OF X LETTERS			
		if (size_word==3):		
			data= re.sub(r'(.{3})(?!$)','\\1 ', data[96:])	
		if (size_word==4):		
			data= re.sub(r'(.{4})(?!$)','\\1 ', data[96:])	
		if (size_word==5):		
			data= re.sub(r'(.{5})(?!$)','\\1 ', data[96:])	
		if (size_word==6):		
			data= re.sub(r'(.{6})(?!$)','\\1 ', data[96:])	
		if (size_word==7):		
			data= re.sub(r'(.{7})(?!$)','\\1 ', data[96:])	
		if (size_word==8):		
			data= re.sub(r'(.{8})(?!$)','\\1 ', data[96:])	
		if (size_word==9):		
			data= re.sub(r'(.{9})(?!$)','\\1 ', data[96:])	
		if (size_word==10):		
			data= re.sub(r'(.{10})(?!$)','\\1 ', data[96:])	
#CHANGE THE TYPE TO LOWER CASE
		data=data.lower()
		NB_WORD_TEST=int(len(data[96:])/size_word)

#IF MORE THAN ONE FILE SELECTED CONCAT THE FILES
		if len(files_for_Word2Vect)>1:		
			data_final=data+data_final 	
		if len(files_for_Word2Vect)==1:		
			data_final=data 	
#SAVE THE FILE TO A GZIP FORMAT COMPATIBLE WITH GENSIM
	if len(files_for_Word2Vect)==1:		
		train_GENSIM_TEST_file=files_for_Word2Vect[i] +'transformed.gz'		
	if len(files_for_Word2Vect)>1:		
		train_GENSIM_TEST_file=mypath+'/Smart_RNA/dataset/TEST DATA/CONCAT_FILE_transformed.gz'			
	with gzip.open(train_GENSIM_TEST_file, 'w') as f:
	    	f.write(data_final.encode())
	sentence=data
	words = sentence.split(' ')
	Word_OK=0
	u=0
	sentence_cleaned=''
	sentence_cleaned_test=''	
	for k,w in enumerate(words):
		Word_OK=0				
		w=w.replace("\r","")
		w=w.replace("\n","")
		for j in range (len(w)):			
			if w[j]=='a' or w[j]=='t' or w[j]=='u' or w[j]=='g' or w[j]=='c':
				Word_OK=Word_OK+1 			
		if Word_OK==size_word:			
			DICT_TEST[u+1,0]=w
			sentence_cleaned_test=sentence_cleaned_test+w
			sentence_cleaned=sentence_cleaned+' '+w
			u=u+1		
			word=w			
			word_test[index_word+1]=w
			index_word=index_word+1			
	words_cleaned = sentence_cleaned.split(' ')

	for k,w in enumerate(words_cleaned):
		word=w		
		positions = [ i for i,w in enumerate(words_cleaned) if w == word ]		
		DICT_TEST[k,1:len(positions)+1]=positions

#CALCULATE STATISTICS
	word_count_test(data_raw)

	Generate_report()

##############################################################################
#		OPEN AND TRANSOFRM BATCH TEST RNA FILE FOR TRAINING SET	     #
##############################################################################
def Open_TEST_BATCH_RNA_FILE():
	global counts_test
	global DICT_TEST
	global train_GENSIM_TEST_file
	global NB_WORD_TEST
	global word_test
	global sentence_cleaned_test
	global index_file_in_test
	global category
	global line_item	
	size_word=4
	index_word=0
#Open multiple files
	line_item=0
	files = askopenfilenames(parent=root,initialdir=mypath+'/Smart_RNA/dataset/'
	                           ,title='Choose files',multiple = 1)
	files_for_training_set = root.tk.splitlist(files)
	files_for_training_set_of_raw=0	
	tot_raw=0
	raw=np.zeros(len(files_for_training_set))
	data_final=''
	for i in range (len(files_for_training_set)):
		category=set_file_category(files_for_training_set[i])
		print(category)		
		index_file_in_test=i
#OPEN THE FILE DOWNLOADED FROM THE DATA BASE
		with open (files_for_training_set[i], "r") as myfile:
		    data=myfile.read()
		data_raw=re.sub(r'(.{1})(?!$)','\\1 ', data[96:])
#SPLIT THE RNA IN WORDS OF X LETTERS			
		if (size_word==3):		
			data= re.sub(r'(.{3})(?!$)','\\1 ', data[96:])	
		if (size_word==4):		
			data= re.sub(r'(.{4})(?!$)','\\1 ', data[96:])	
		if (size_word==5):		
			data= re.sub(r'(.{5})(?!$)','\\1 ', data[96:])	
		if (size_word==6):		
			data= re.sub(r'(.{6})(?!$)','\\1 ', data[96:])	
		if (size_word==7):		
			data= re.sub(r'(.{7})(?!$)','\\1 ', data[96:])	
		if (size_word==8):		
			data= re.sub(r'(.{8})(?!$)','\\1 ', data[96:])	
		if (size_word==9):		
			data= re.sub(r'(.{9})(?!$)','\\1 ', data[96:])	
		if (size_word==10):		
			data= re.sub(r'(.{10})(?!$)','\\1 ', data[96:])	
#CHANGE THE TYPE TO LOWER CASE
		data=data.lower()
		NB_WORD_TEST=int(len(data[96:])/size_word)
#IF MORE THAN ONE FILE SELECTED CONCAT THE FILES
		#if len(files_for_training_set)>1:		
		#	data_final=data+data_final 	
		#if len(files_for_training_set)==1:		
		data_final=data 	
		sentence=data
		words = sentence.split(' ')
		Word_OK=0
		u=0
		sentence_cleaned=''
		sentence_cleaned_test=''
		index_word=0	
		for k,w in enumerate(words):
			Word_OK=0				
			w=w.replace("\r","")
			w=w.replace("\n","")
			for j in range (len(w)):			
				if w[j]=='a' or w[j]=='t' or w[j]=='u' or w[j]=='g' or w[j]=='c':
					Word_OK=Word_OK+1 			
			if Word_OK==size_word:			
				DICT_TEST[u+1,0]=w
				sentence_cleaned_test=sentence_cleaned_test+w
				sentence_cleaned=sentence_cleaned+' '+w
				u=u+1		
				word=w			
				word_test[index_word+1]=w
				index_word=index_word+1			
		words_cleaned = sentence_cleaned.split(' ')
		for k,w in enumerate(words_cleaned):
			word=w		
			positions = [ i for i,w in enumerate(words_cleaned) if w == word ]		
			DICT_TEST[k,1:len(positions)+1]=positions
		chain_similarity()
	#RCE_NN_Learning()


def set_file_category(file_name):
	category=0	
	print(file_name)
	if "US" in file_name and "HU" in file_name and "1" in file_name:
		print('US HUMAN PERIOD 1')
		category=1
	if "US" in file_name and "HU" in file_name and "2" in file_name:
		print('US HUMAN PERIOD 2')
		category=2
	if "US" in file_name and "HU" in file_name and "3" in file_name:
		print('US HUMAN PERIOD 3')
		category=3
	if "US" in file_name and "HU" in file_name and "4" in file_name:
		print('US HUMAN PERIOD 4')
		category=4

	if "EU" in file_name and "HU" in file_name and "1" in file_name:
		print('EU HUMAN PERIOD 1')
		category=5
	if "EU" in file_name and "HU" in file_name and "2" in file_name:
		print('EU HUMAN PERIOD 2')
		category=6
	if "EU" in file_name and "HU" in file_name and "3" in file_name:
		print('EU HUMAN PERIOD 3')
		category=7
	if "EU" in file_name and "HU" in file_name and "4" in file_name:
		print('EU HUMAN PERIOD 4')
		category=8

	if "CHINA" in file_name and "HU" in file_name and "1" in file_name:
		print('CHINA HUMAN PERIOD 1')
		category=9
	if "CHINA" in file_name and "HU" in file_name and "2" in file_name:
		print('CHINA HUMAN PERIOD 2')
		category=10
	if "CHINA" in file_name and "HU" in file_name and "3" in file_name:
		print('CHINA HUMAN PERIOD 3')
		category=11
	if "CHINA" in file_name and "HU" in file_name and "4" in file_name:
		print('CHINA HUMAN PERIOD 4')
		category=12

	if "SA" in file_name and "HU" in file_name and "1" in file_name:
		print('SOUTH AMERICA HUMAN PERIOD 1')
		category=13
	if "SA" in file_name and "HU" in file_name and "2" in file_name:
		print('SOUTH AMERICA HUMAN PERIOD 2')
		category=14
	if "SA" in file_name and "HU" in file_name and "3" in file_name:
		print('SOUTH AMERICA HUMAN PERIOD 3')
		category=15
	if "SA" in file_name and "HU" in file_name and "4" in file_name:
		print('SOUTH AMERICA HUMAN PERIOD 4')
		category=16

	if "ASIA" in file_name and "HU" in file_name and "1" in file_name:
		print('ASIA HUMAN PERIOD 1')
		category=17
	if "ASIA" in file_name and "HU" in file_name and "2" in file_name:
		print('ASIA HUMAN PERIOD 2')
		category=18
	if "ASIA" in file_name and "HU" in file_name and "3" in file_name:
		print('ASIA HUMAN PERIOD 3')
		category=19
	if "ASIA" in file_name and "HU" in file_name and "4" in file_name:
		print('ASIA HUMAN PERIOD 4')
		category=20

	if "AF" in file_name and "HU" in file_name and "1" in file_name:
		print('AF HUMAN PERIOD 1')
		category=21
	if "AF" in file_name and "HU" in file_name and "2" in file_name:
		print('AF HUMAN PERIOD 2')
		category=22
	if "AF" in file_name and "HU" in file_name and "3" in file_name:
		print('AF HUMAN PERIOD 3')
		category=23
	if "AF" in file_name and "HU" in file_name and "4" in file_name:
		print('AF HUMAN PERIOD 4')
		category=24

	if "IND" in file_name and "HU" in file_name and "1" in file_name:
		print('INDIA HUMAN PERIOD 1')
		category=25
	if "IND" in file_name and "HU" in file_name and "2" in file_name:
		print('INDIA HUMAN PERIOD 2')
		category=26
	if "IND" in file_name and "HU" in file_name and "3" in file_name:
		print('INDIA HUMAN PERIOD 3')
		category=27
	if "IND" in file_name and "HU" in file_name and "4" in file_name:
		print('INDIA HUMAN PERIOD 4')
		category=28
	if "HU" in file_name:
		print('HUMAN')
	else:
		print('ANIMAL')
		category=100

	return category 
def set_chain_similarity():
	print('set chain similarity')
	#Similarity_progressbar.grid(row=6, column=1)
	#Similarity_progressbar.update()
	#time.sleep(0.1)
	#Similarity_progressbar.update()		
	chain_similarity()
def chain_similarity():
	global category
	global DICT_TEST
	global DICT_REF
	global NB_WORD_TEST
	global NB_WORD_REF
	global word_ref
	global word_test
	global sentence_cleaned_test
	global sentence_cleaned_ref
	global index_file_in_test
	global line_item
	global feature_line
	global ANNOTATION_MATRIX	
	global RNA_AI_RECO

	similarity=np.zeros((10000))
	Seq_TEST_PREV=np.zeros((10000), dtype=object )								
	Seq_REF_PREV=np.zeros((10000), dtype=object )								
	Seq_REF_NEXT=np.zeros((10000), dtype=object )
	Seq_TEST_NEXT=np.zeros((10000), dtype=object )
	sim_seq_test=np.zeros((10000), dtype=object )
	sim_seq_ref=np.zeros((10000), dtype=object )
	diff_seq_ref=np.zeros((10000,10), dtype=object )
	diff_seq_test=np.zeros((10000,10), dtype=object )	
	diff_seq_test_neg=np.zeros((10000,10), dtype=object )
	diff_seq_ref_neg=np.zeros((10000,10), dtype=object )	
	MAT_WORD_ALREADY_TESTED=np.zeros((10000,1), dtype=object )		
	POS_DIFF_REF=0
	POS_DIFF_TEST=0		
	similarity_count=5
	Max_similarity=6
	Target_similarity_ratio=50
	count_hih_similarity_rate=0	
	Max_Glitch=3000
	print('chain similarity')
	COUNT_NEG=0
	COUNT_POS=0
	currentValue=1
	divisions=1
	maxValue=100
	Similarity_progressbar["value"]=currentValue
	Similarity_progressbar["maximum"]=maxValue
	aligner = Align.PairwiseAligner()
	aligner.mode = 'local'
	print('*****************************************')
	print('			USER DATA			')
	print('*****************************************')
	seq1= sentence_cleaned_test.lower()
	seq2= sentence_cleaned_ref.lower()
	nb_subseq_seq1=int(len(seq1)/100)
	nb_subseq_seq2=int(len(seq2)/100)
	if nb_subseq_seq2>=nb_subseq_seq1:
		nb_subseq=nb_subseq_seq1
	if nb_subseq_seq1>nb_subseq_seq2:
		nb_subseq=nb_subseq_seq2
	
	score=np.zeros((nb_subseq))
	for i in range (nb_subseq):
		score[i] = aligner.score(seq1[i*100:i*100+100], seq2[i*100:i*100+100])
	maxi=0	
	for i in range (nb_subseq):
		if score[i]>maxi:
			maxi=score[i]
			index_max=i

	alignments = aligner.align(seq1,seq2)
	for i in range(1):	
		align_str=str(alignments[i])
		pos_mut=np.zeros((int(len(align_str)/3)))
		offset=int(2*len(align_str)/3)
		match=0
		Frame_shift=0
		Missense=0
		j=0		
		for index in range (int(len(align_str)/3)):
			if align_str[index]=='-' or align_str[index+offset]=='-':
				#print('Frameshift Mutation')
				Frame_shift=Frame_shift+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]==align_str[index+offset]:
				#print('match')
				match=match+1 
			if align_str[index]=='a' and align_str[index+offset]=='t':
				print('Missense Mutation A -> T')
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='a' and align_str[index+offset]=='u':
				print('Missense Mutation A -> U')
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='a' and align_str[index+offset]=='g':
				print('Missense Mutation A -> G')
				Missense=Missense+1		
				pos_mut[j]=index
				j=j+1
			if align_str[index]=='a' and align_str[index+offset]=='c':
				print('Missense Mutation A -> C')
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1
			
			if align_str[index]=='t' and align_str[index+offset]=='a':
				print('Missense Mutation T -> A')
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='t' and align_str[index+offset]=='u':
				print('Missense Mutation T -> U')
				Missense=Missense+1		
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='t' and align_str[index+offset]=='g':
				print('Missense Mutation T -> G')
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='t' and align_str[index+offset]=='c':
				print('Missense Mutation T -> C')
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='u' and align_str[index+offset]=='a':
				print('Missense Mutation U -> A')
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='u' and align_str[index+offset]=='t':
				print('Missense Mutation U -> T')		
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='u' and align_str[index+offset]=='g':
				print('Missense Mutation U -> G')
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='u' and align_str[index+offset]=='c':
				print('Missense Mutation U -> C')
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='g' and align_str[index+offset]=='a':
				print('Missense Mutation G -> A')
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='g' and align_str[index+offset]=='t':
				print('Missense Mutation G -> T')
				Missense=Missense+1
				pos_mut[i]=index
				i=i+1

			if align_str[index]=='g' and align_str[index+offset]=='u':
				print('Missense Mutation G -> U')
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='g' and align_str[index+offset]=='c':
				print('Missense Mutation G -> C')
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='c' and align_str[index+offset]=='a':
				print('Missense Mutation C -> A')

				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='c' and align_str[index+offset]=='t':
				print('Missense Mutation C -> T')
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='c' and align_str[index+offset]=='u':
				print('Missense Mutation C -> U')
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

			if align_str[index]=='c' and align_str[index+offset]=='g':
				print('Missense Mutation C -> G')
				Missense=Missense+1
				pos_mut[j]=index
				j=j+1

		for index in range (int(len(align_str)/3)-2):
			if align_str[index]=='u' and align_str[index+1]=='a' and align_str[index+2]=='g':
				print('Stop sequence detected')
		for index in range (int(len(align_str)/3)-2):
			if align_str[index]=='u' and align_str[index+1]=='a' and align_str[index+2]=='a':
				print('Stop sequence detected')
		for index in range (int(len(align_str)/3)-2):
			if align_str[index]=='u' and align_str[index+1]=='g' and align_str[index+2]=='a':
				print('Stop sequence detected')
		
		for index in range (int(len(align_str)/3)-2):
			if align_str[index+offset]=='u' and align_str[index+offset+1]=='a' and align_str[index+offset+2]=='g':
				print('Stop sequence detected')
		for index in range (int(len(align_str)/3)-2):
			if align_str[index+offset]=='u' and align_str[index+offset+1]=='a' and align_str[index++offset+2]=='a':
				print('Stop sequence detected')
		for index in range (int(len(align_str)/3)-2):
			if align_str[index+offset]=='u' and align_str[index+offset+1]=='g' and align_str[index+offset+2]=='a':
				print('Stop sequence detected')
		
		print('missense Mutation',Missense)
		print('Frame Shift Mutation',Frame_shift)
		print('match', match)

		homogeneity=100*match/int(len(align_str)/3)
		print('Homogeneity :' + str(homogeneity)+'%')
		#print(align_str[1+ int(2*len(align_str)/3)])
		print("Done")
	#for alignment in alignments:
	#	print(alignment)
		
		text_report.insert(INSERT, "Homogeneity "+ str(homogeneity)+"%" + '\n')
		text_report.insert(INSERT, "Missense Mutation "+ str(Missense) + '\n')
		text_report.insert(INSERT, "Frame Shift Mutation "+ str(Frame_shift) + '\n')
		text_report.insert(INSERT, "Nucleotide Match "+ str(match) + '\n')
		text_report.insert(INSERT, '\n')
	Result_mat=pd.DataFrame(columns=['word_tested', 'Position_TEST', 'Position_REF', 'COUNT_POS','COUNT_NEG', 'Size_Sequence', 'Similarity_Ratio (No Gitch)', 'Similarity_Ratio'])
	word_already_present=0
	pd_index=0
	np.random.seed(0)
	import matplotlib.pyplot as pltt
	X_tot=pos_mut
	nb_mutation=j
	from sklearn.mixture import GaussianMixture
	random_state = np.random.RandomState(seed=1)
	X =  np.concatenate([X_tot[0:60]]).reshape(-1, 1)
	X1=X_tot[0:nb_mutation].reshape(-1, 1)
	N = np.arange(1, 11)
	models = [None for i in range(len(N))]

	for i in range(len(N)):
		models[i] = GaussianMixture(N[i]).fit(X)
	

	# compute the AIC and the BIC
	AIC = [m.aic(X) for m in models]
	BIC = [m.bic(X) for m in models]
	fig = plt.figure(figsize=(8, 6))
	fig.subplots_adjust(left=0.12, right=0.97,
		            bottom=0.21, top=0.9, wspace=0.5)
	ax = fig.add_subplot(111)
	M_best = models[np.argmin(AIC)]
	
	x = np.linspace(0, 100, 6000)
	logprob = M_best.score_samples(x.reshape(-1, 1))
	responsibilities = M_best.predict_proba(x.reshape(-1, 1))
	pdf = np.exp(logprob)
	pdf_individual = responsibilities * pdf[:, np.newaxis]

	ax.hist(X, 30, density=True, histtype='stepfilled', alpha=0.4)
	from scipy.stats import kde
	if nb_mutation>0:	
		density = kde.gaussian_kde(X_tot[0:nb_mutation]) # x: list of price
		xgrid = np.linspace(0, int(len(align_str)/3), nb_mutation)   
		plt.plot(xgrid, 5*density(xgrid))
	ax.set_xlabel('$x$')
	ax.set_ylabel('$p(x)$')
	scatter1 = FigureCanvasTkAgg(fig, frame16) 
	scatter1._tkcanvas.grid(row =0, column =0)
	ax.legend() 
	scatter1.draw()

	from sklearn.cluster import MeanShift, estimate_bandwidth

	cluster_index=0
	cluster=0
	item=0
	sum_pos=0
	if nb_mutation>0:	
		cluster=np.zeros((nb_mutation))
		item_in_cluster=np.zeros((nb_mutation))
		cluster_center=np.zeros((nb_mutation))
		for i in range (nb_mutation):
			if X_tot[i]-X_tot[i-1]<500:
				cluster[i]=cluster_index
				item_in_cluster[cluster_index]=item_in_cluster[cluster_index]+1
				sum_pos=sum_pos+X_tot[i]
				cluster_center[cluster_index]=(sum_pos)/item_in_cluster[cluster_index]
			else:
				sum_pos=0		
				cluster_index=cluster_index+1
				cluster[i]=cluster_index
				item_in_cluster[cluster_index]=item_in_cluster[cluster_index]+1
				sum_pos=sum_pos+X_tot[i]
				cluster_center[cluster_index]=(sum_pos)/item_in_cluster[cluster_index]

				item=0
		cluster_center_pos_matrix=np.zeros((255))
		cluster_mut_type_matrix=np.zeros((255))
		cluster_quantity_matrix=np.zeros((255))
		for i in range (cluster_index+1):
		#	print('cluster:', i)
		#	print('center cluster :',cluster_center[i])
		#	print('Nb element in cluster :',item_in_cluster[i])
			cluster_center_pos_matrix[int(cluster_center[i]/500)]=int(cluster_center[i])
			cluster_quantity_matrix[int(cluster_center[i]/500)]=int(item_in_cluster[i])		
		set_cat=2
		feature_line[line_item,10:265]=cluster_center_pos_matrix[:]
		feature_line[line_item+1,10:265]=cluster_quantity_matrix[:] 
		#feature_line[line_item+2,10:265]=cluster_mut_type_matrix[:] 			
		feature_line[line_item,0]=1   #index
		feature_line[line_item+1,0]=1   #index
		#feature_line[line_item+2,0]=1   #index
			
		feature_line[line_item,6]=1   #Context	
		feature_line[line_item+1,6]=2  #Context	
		#feature_line[line_item+2,6]=3   #Context	
		#line_item=line_item+3
		
		if RNA_AI_RECO==0: 
			feature_line[line_item,5]=category   #Category
			feature_line[line_item+1,5]=category   #Category
			#feature_line[line_item+2,5]=category  #Category
			ANNOTATION_MATRIX=feature_line

			title=mypath+'Smart_RNA/FEATURE FILE/training_set'

		if RNA_AI_RECO==1: 
			title=mypath+'Smart_RNA/FEATURE FILE/testing_set'
		#	Reco_by_feature_file(title)			
		np.savetxt(title, feature_line,'%f',delimiter=',')
		line_item=line_item+2



def Set_RNA_AI_RECO():
	global RNA_AI_RECO
	print('AI RECO')
	RNA_AI_RECO=1
	Open_TEST_BATCH_RNA_FILE()

def Set_RNA_AI_TRAINING():
	global RNA_AI_RECO
	print('AI TRAINING')
	RNA_AI_RECO=0
	Open_TEST_BATCH_RNA_FILE()		
		
def Set_RNA_AI_RECO_FILE():
	global RNA_AI_RECO
	print('AI RECO')
	files = askopenfilename(parent=root,initialdir=mypath+'Smart_RNA/FEATURE FILE/'
                           ,title='Choose files',multiple = 1)		
	files_for_test = root.tk.splitlist(files)
	test_file=files_for_test[0]	
	Reco_by_feature_file(test_file)
	

def Set_RNA_AI_TRAINING_FILE():
	global ANNOTATION_MATRIX
	global RNA_AI_RECO
	print('AI TRAINING')
	files = askopenfilenames(parent=root,initialdir=mypath+'Smart_RNA/FEATURE FILE/'
                           ,title='Choose files',multiple = 1)
	files_for_test = root.tk.splitlist(files)
	print ("list of feature files =",files_for_test)	
	print(len(files_for_test))
	for i in range (len(files_for_test)):
		Input_MATRIX=np.loadtxt(files_for_test[i],delimiter=',')
	normalized_matrix=Normalize(Input_MATRIX)
	ANNOTATION_MATRIX=normalized_matrix
	RCE_NN_Learning()

def Normalize(Matrix):
	max_cont1=0
	max_cont2=0

	for i in range (Matrix.shape[0]):
		for j in range(256):
			if Matrix[i, 6]==1:
				if Matrix[i, 10+j]>max_cont1:
					max_cont1=Matrix[i, 10+j]
			if Matrix[i, 6]==2:
				if Matrix[i, 10+j]>max_cont2:
					max_cont2=Matrix[i, 10+j]
	for i in range (Matrix.shape[0]):
		for j in range(256):
			if Matrix[i, 6]==1:
				Matrix[i, 10+j]=int(255*Matrix[i, 10+j]/max_cont1)
			if Matrix[i, 6]==2:
				Matrix[i, 10+j]=int(255*Matrix[i, 10+j]/max_cont2)
	return(Matrix)
	

def Generate_report():
	global counts_test	
	global counts_ref
	total_nucleotide_test=counts_test['A']+counts_test['C']+counts_test['T']+counts_test['G']
	total_nucleotide_ref= counts_ref['A']+counts_ref['C']+counts_ref['T']+counts_ref['G']
	text_report.insert(INSERT, "number of Nucleotide in the reference sequence   :"+ str(total_nucleotide_test) + '\n')
	text_report.insert(INSERT, "A :"+ str(counts_ref['A']) + '\n')
	text_report.insert(INSERT, "C :"+ str(counts_ref['C']) + '\n')
	text_report.insert(INSERT, "T :"+ str(counts_ref['T']) + '\n')
	text_report.insert(INSERT, "G :"+ str(counts_ref['G']) + '\n')
	text_report.insert(INSERT, "number of Nucleotide in the test sequence        :"+ str(total_nucleotide_ref) + '\n')
	text_report.insert(INSERT, "A :"+ str(counts_test['A']) + '\n')
	text_report.insert(INSERT, "C :"+ str(counts_test['C']) + '\n')
	text_report.insert(INSERT, "T :"+ str(counts_test['T']) + '\n')
	text_report.insert(INSERT, "G :"+ str(counts_test['G']) + '\n')	


def Generate_report_similarity():
	global counts_test	
	global counts_ref
	total_nucleotide_test=counts_test['A']+counts_test['C']+counts_test['T']+counts_test['G']
	total_nucleotide_ref= counts_ref['A']+counts_ref['C']+counts_ref['T']+counts_ref['G']
	text_report.insert(INSERT, "number of Nucleotide in the reference sequence   :"+ str(total_nucleotide_test) + '\n')
	text_report.insert(INSERT, "A :"+ str(counts_ref['A']) + '\n')
	text_report.insert(INSERT, "C :"+ str(counts_ref['C']) + '\n')
	text_report.insert(INSERT, "T :"+ str(counts_ref['T']) + '\n')
	text_report.insert(INSERT, "G :"+ str(counts_ref['G']) + '\n')
	text_report.insert(INSERT, "number of Nucleotide in the test sequence        :"+ str(total_nucleotide_ref) + '\n')
	text_report.insert(INSERT, "A :"+ str(counts_test['A']) + '\n')
	text_report.insert(INSERT, "C :"+ str(counts_test['C']) + '\n')
	text_report.insert(INSERT, "T :"+ str(counts_test['T']) + '\n')
	text_report.insert(INSERT, "G :"+ str(counts_test['G']) + '\n')	
											 
def word_count_test(str):
	global counts_test
	global np_counts_test	
	
	
	words=str.split()
	for word in words:
		if word=='A' or word=='C' or word=='T' or word=='G':			
			if word in counts_test:
				counts_test[word]=0
	for word in words:
		if word=='A' or word=='C' or word=='T' or word=='G':			
			if word in counts_test:
				counts_test[word]+=1
			else:
				counts_test[word]=0
	
	data = list(counts_test.items())
	np_counts_test = np.array(data)	
	
def word_count_ref(str):
	global counts_ref
	global np_counts_ref	
	
		
	words=str.split()
	for word in words:
		if word=='A' or word=='C' or word=='T' or word=='G':			
			if word in counts_ref:
				counts_ref[word]=0
		
	for word in words:
		if word=='A' or word=='C' or word=='T' or word=='G':			
			if word in counts_ref:
				counts_ref[word]+=1
			else:
				counts_ref[word]=0
	
	print(counts_ref)	
	keys=counts_ref.keys()
	values=counts_ref.values()
	data = list(counts_ref.items())
	np_counts_ref = np.array(data)	

def occurence_chart():
	global counts_test
	global counts_ref
	global np_counts_ref
	global np_counts_test
	values=counts_ref.values()
	figure1 = plt.Figure(figsize=(8,6))
	figure1, (ax, ax2) = plt.subplots(2, 1)	
	#ax = figure1.add_subplot(111)
	scatter1 = FigureCanvasTkAgg(figure1,frame16) 
	scatter1._tkcanvas.grid(row =0, column =0)
	ax.legend() 
	labels=[str(np_counts_ref[0,0]),str(np_counts_ref[1,0]), str(np_counts_ref[2,0]), str(np_counts_ref[3,0])]	
	vals=[int(np_counts_ref[0,1]),int(np_counts_ref[1,1]), int(np_counts_ref[2,1]), int(np_counts_ref[3,1])]	
	ax.bar(labels,vals)
	# Call the function above. All the magic happens there.
	add_value_labels(ax)		
	ax.set_title('Occurence Nucleotide REFERENCE data')
	ax.set_ylim(0, 10000)
	

	ax2.legend() 
	values=counts_test.values()
	labels=[str(np_counts_test[0,0]),str(np_counts_test[1,0]), str(np_counts_test[2,0]), str(np_counts_test[3,0])]	
	vals=[int(np_counts_test[0,1]),int(np_counts_test[1,1]), int(np_counts_test[2,1]), int(np_counts_test[3,1])]	
	
	#ax = figure1.add_subplot(112)
	#scatter1 = FigureCanvasTkAgg(figure1, frame16) 
	#scatter1._tkcanvas.grid(row =1, column =0)
	x = np.arange(4)		
	ax2.bar(labels,vals)
	ax2.set_xlabel('nucleotide')
	ax2.set_title('Occurence Nucleotide TEST data')
	#ax.set_xlim(0, 4)
	ax2.set_ylim(0, 10000)
	# Call the function above. All the magic happens there.
	add_value_labels(ax2)
	figure1.tight_layout()
	scatter1.draw()

def add_value_labels(ax, spacing=5):
    """Add labels to the end of each bar in a bar chart.

    Arguments:
        ax (matplotlib.axes.Axes): The matplotlib object containing the axes
            of the plot to annotate.
        spacing (int): The distance between the labels and the bars.
    """

    # For each bar: Place a label
    for rect in ax.patches:
        # Get X and Y placement of label from rect.
        y_value = rect.get_height()
        x_value = rect.get_x() + rect.get_width() / 2

        # Number of points between bar and label. Change to your liking.
        space = spacing
        # Vertical alignment for positive values
        va = 'bottom'

        # If value of bar is negative: Place label below bar
        if y_value < 0:
            # Invert space to place label below
            space *= -1
            # Vertically align label at top
            va = 'top'

        # Use Y value as label and format number with one decimal place
        label = "{:.1f}".format(y_value)

        # Create annotation
        ax.annotate(
            label,                      # Use `label` as label
            (x_value, y_value),         # Place label at end of the bar
            xytext=(0, space),          # Vertically shift label by `space`
            textcoords="offset points", # Interpret `xytext` as offset in points
            ha='center',                # Horizontally center label
            va=va)                      # Vertically align label differently for
                                        # positive and negative values.





##############################################################
#		TRAIN GENSIM FOR REF AND TEST RNA FILE	     #
##############################################################

def Train_GENSIM_REF_TEST_RNA_FILE():
#Setting data for GENSIM
	size_word=3
	feature_number=256
	size_window=8
	epochs_number=300
	global train_GENSIM_REF_file
	global train_GENSIM_TEST_file
#READ CREATED REF DOCUMENT WITH GENSIM AND PRE-TREAT THE DOCUMENT TO CREATE WORD TOKEN	
	doc_REF = list(read_input(train_GENSIM_REF_file))
	logging.info("Done reading data file for REF data")
	print(doc_REF)
#TRAIN THE GENSIN WORD2VECT WITH GENSIM AND PRE-TREAT THE DOCUMENT TO CREATE WORD TOKEN FOR REF
	model_REF = gensim.models.Word2Vec(doc_REF,size=feature_number,window=size_window,min_count=1,workers=10)
	model_REF.train(doc_REF, total_examples=len(doc_REF), epochs=epochs_number)
#SAVE VECTORS CREATED BY GENSIM FOR REF
	model_REF.wv.save(mypath+'Smart_RNA/dataset/REFERENCE DATA/vectors/REF_Vect_default')
#READ CREATED TEST DOCUMENT WITH GENSIM AND PRE-TREAT THE DOCUMENT TO CREATE WORD TOKEN	
	doc_TEST = list(read_input(train_GENSIM_TEST_file))
	logging.info("Done reading data file for TEST data")
	print(doc_TEST)
#TRAIN THE GENSIN WORD2VECT WITH GENSIM AND PRE-TREAT THE DOCUMENT TO CREATE WORD TOKEN
	model_TEST = gensim.models.Word2Vec(doc_TEST,size=feature_number,window=size_window,min_count=1,workers=10)
	model_TEST.train(doc_REF, total_examples=len(doc_TEST), epochs=epochs_number)
#SAVE VECTORS CREATED BY GENSIM
	model_TEST.wv.save(mypath+'Smart_RNA/dataset/TEST DATA/vectors/TEST_Vect_default')
	w1 = "cac"
	print("Most similar to {0}".format(w1), model_REF.wv.most_similar(positive=w1))
	print(
	"Most similar to {0}".format(w1),
		model_REF.wv.most_similar(
		positive=w1,
		topn=15))
	print("Similarity between 'ttt' and 'gat'",model_REF.wv.similarity(w1="ttt", w2="gat"))
	print('word vector cac:', model_REF["cac"])

	w1 = "cac"
	print("Most similar to {0}".format(w1), model_TEST.wv.most_similar(positive=w1))
	print(
	"Most similar to {0}".format(w1),
		model_TEST.wv.most_similar(
		positive=w1,
		topn=15))
	print("Similarity between 'ttt' and 'gat'",model_TEST.wv.similarity(w1="ttt", w2="gat"))
	print('word vector cac:', model_TEST["cac"])


logging.basicConfig(
    format='%(asctime)s : %(levelname)s : %(message)s',
    level=logging.INFO)

	
def progress(currentValue):
	Train_progressbar["value"]=currentValue
	GA_progressbar["value"]=currentValue	
	Similarity_progressbar["value"]=currentValue
def Set_Reco_NM():
	frame1.grid_remove()
	frame3.grid_remove()
	frame2.grid_remove()
	frame4.grid_remove()
	frame5.grid_remove()
	frame6.grid_remove()
	frame7.grid_remove()
	frame8.grid_remove()
	frame9.grid_remove()	
	frame10.grid_remove()
	frame11.grid_remove()
	frame12.grid_remove()
	frame13.grid_remove()
	frame14.grid_remove()	
	frame15.grid_remove()
	frame16.grid_remove()
	frame17.grid_remove()
	frame18.grid_remove()
	frame19.grid_remove()
	frame20.grid_remove()
	frame30.grid_remove()
	Frame_Set_Ga.grid_remove()
	frame21.grid(row =1, column =0)	
	frame22.grid(row =1, column =4)
	frame23.grid(row =1, column =5)
	Button_Start_Reco_NM.grid(row =0, column =0)
	Button_Stop_Reco_NM.grid(row =1, column =0)
	Button_Reco_by_feature.grid(row =2, column =0)
	Button_Set_Maxif_NM.grid(row =3, column =0)
	


##############################################################
#		B SPECIES PREDICTION 			     #
#		PROTEIN PROTEIN INTERACTION		     #	
#							     #	
# 							     #	
#							     #	
#						             #
##############################################################
path_sentence=''
def set_protein_comp():
	
	protein_seq_homogeneity()
	
def protein_seq_homogeneity():
	global prot_sim_result
	similarity=np.zeros((10000))
	Seq_TEST_PREV=np.zeros((10000), dtype=object )								
	Seq_REF_PREV=np.zeros((10000), dtype=object )								
	Seq_REF_NEXT=np.zeros((10000), dtype=object )
	Seq_TEST_NEXT=np.zeros((10000), dtype=object )
	sim_seq_test=np.zeros((10000), dtype=object )
	sim_seq_ref=np.zeros((10000), dtype=object )
	diff_seq_ref=np.zeros((10000,10), dtype=object )
	diff_seq_test=np.zeros((10000,10), dtype=object )	
	diff_seq_test_neg=np.zeros((10000,10), dtype=object )
	diff_seq_ref_neg=np.zeros((10000,10), dtype=object )	
	MAT_WORD_ALREADY_TESTED=np.zeros((10000,1), dtype=object )		
	POS_DIFF_REF=0
	POS_DIFF_TEST=0		
	similarity_count=5
	Max_similarity=6
	Target_similarity_ratio=50
	count_hih_similarity_rate=0	
	Max_Glitch=3000
	print('chain similarity')
	COUNT_NEG=0
	COUNT_POS=0
	currentValue=1
	divisions=1
	maxValue=100
	Similarity_progressbar["value"]=currentValue
	Similarity_progressbar["maximum"]=maxValue
	aligner = Align.PairwiseAligner()
	aligner.mode = 'local'
	print('*****************************************')
	print('			USER DATA			')
	print('*****************************************')
	files = askopenfilename(parent=root,initialdir=mypath+'Smart_RNA/dataset/',title='Choose file containing protein sequence to be tested',multiple = 1)		
	files_for_test = root.tk.splitlist(files)
	test_prot=files_for_test[0]	
	with open(test_prot, 'r', encoding = "ISO-8859-1") as protfile:
		protein_seq=protfile.read()
	print(protein_seq)
	protein_seq = ''.join(filter(str.isalnum, protein_seq))
	text_report.insert(INSERT, "Protein sequence under test :" + protein_seq + '\n')
	column={'# protein_xref_1','protein_xref_2','alternative_identifiers_1','alternative_identifiers_2','protein_alias_1','protein_alias_2','detection_method','author_name				#pmid','protein_taxid_1','protein_taxid_2','interaction_type','source_database_id','database_identifier','confidence','protein_xref_1_unique','protein_xref_2_unique','protein_taxid_1_cat',	'protein_taxid_2_cat',	'protein_taxid_1_name						#protein_taxid_2_name','protein_seq1','protein_seq2','source_database','protein_xref_1_display_id','protein_xref_2_display_id'}
	df = pd.read_csv(mypath+'Smart_RNA/dataset/PH_Virus/hpidb2.mitab_plus.txt',  sep='\t', encoding = "ISO-8859-1", header=(0)) 
	seq2=str(protein_seq)
	idx = df.shape[0]
	ind=0
	for k in range(1,idx):
		seq1=df.loc[k, 'protein_seq1']
		if len(seq1)>len(seq2):
			min_len=len(seq2)
		else:
			min_len=len(seq1)
		
		seq1=seq1.lower()
		seq2= seq2.lower()
		maxi=0	
		alignments = aligner.align(seq1,seq2)
		for i in range(1):	
			align_str=str(alignments[i])
			offset=int(2*len(align_str)/3)
			match=0
			Frame_shift=0
			Missense=0
			j=0		
			for index in range (int(len(align_str)/3)):
				if align_str[index]==align_str[index+offset]:
					#print('match')
					match=match+1 
			homogeneity=100*match/int(min_len+10)
						
			if homogeneity>95.8:
				prot_sim_result[ind]=df.loc[k, 'protein_xref_2_unique']			
				if "NCBI_ACC:" in prot_sim_result[ind]:				
					prot_sim_result[ind]=prot_sim_result[ind].replace("NCBI_ACC:", "")
				if "UNIPROT_AC:" in prot_sim_result[ind]:				
					prot_sim_result[ind]=prot_sim_result[ind].replace("UNIPROT_AC:", "")


								
				print('Homogeneity :' + str(homogeneity)+'% ' + str(prot_sim_result[ind]))
				text_report.insert(INSERT, "Homogeneity :" + str(homogeneity)+'% ' + str(prot_sim_result[ind]) + '\n')
				ind=ind+1	
	print('done')	
	
#		text_report.insert(INSERT, "Alignment ID "+ str(i) + '\n')
#		text_report.insert(INSERT, "Homogeneity "+ str(Missense) + '\n')
	
#		text_report.insert(INSERT, "Missense Mutation "+ str(Missense) + '\n')
#		text_report.insert(INSERT, "Frame Shift Mutation "+ str(Frame_shift) + '\n')
#		text_report.insert(INSERT, "Nucleotide Match "+ str(match) + '\n')
#		text_report.insert(INSERT, '\n')

#Global Variable Definition
Virus_Axiom_File_Name =''
Species_Axiom_File_Name =''

def Create_Species_Pheno_FILE():
	print("compute vectors for species phenotype data")
	Handle_Species_axioms_files()
def Create_Virus_Pheno_FILE():
	print("compute vectors for Virus phenotype data")
	Handle_Virus_axioms_files()


def Handle_Virus_axioms_files():
	global Virus_Axiom_File_Name

#Open multiple files
	files = askopenfilenames(parent=root,initialdir=mypath+'Smart_RNA/dataset/PH_Virus/'
	                           ,title='Choose files',multiple = 1)
	Virus_axiom = root.tk.splitlist(files)
	raw=np.zeros(len(Virus_axiom))
	data_final=''
	with open(mypath+'Smart_RNA/dataset/PH_Virus/Axiom_Virus', 'w') as outfile:
		for i in range (len(Virus_axiom)):
			print('file ', i)
			with open(Virus_axiom[i]) as infile:
				for line in infile:
					outfile.write(line)

	Virus_Axiom_File_Name=mypath+'Smart_RNA/dataset/PH_Virus/Axiom_Virus'
def Handle_Species_axioms_files():
	global Species_Axiom_File_Name
#Open multiple files
	files = askopenfilenames(parent=root,initialdir=mypath+'Smart_RNA/dataset/PH_Species'
	                           ,title='Choose files',multiple = 1)
	species_axiom = root.tk.splitlist(files)
	raw=np.zeros(len(species_axiom))
	data_final=''
	with open(mypath+'Smart_RNA/dataset/PH_Species/Axiom_Species', 'w') as outfile:
		for i in range (len(species_axiom)):
			print('file ', i)
			with open(species_axiom[i]) as infile:
				for line in infile:
					outfile.write(line)
	Species_Axiom_File_Name=mypath+'Smart_RNA/dataset/PH_Species/Axiom_Species'

##############################################################
#	HUMAN/SPECIES PHENOTYPES - GENES ASSOCIATION	     #
#							     #
##############################################################

# generate the human disease and hp association

def Create_Species_Pheno_protein_asso():
	print('Create Species Pheno/protein association')
	dis_phe=dict()
	with open(mypath+'Smart_RNA/dataset/PH_Species/phenotype_annotation.tab',"r") as f:
		for line in f.readlines():
			data=line.split("\t")	
			if (data[0]=="OMIM")&(data[5][:4]=="OMIM"):
				try:
					dis_phe[data[5].strip()].append(data[4].strip())
				except:
					dis_phe[data[5].strip()]=[data[4].strip()]	
			
	print("the number of disease phenotypes ",len(dis_phe))		
	# obtain the mouse gene and mouse disease association
	# obtain the human gene and human disease association
	mgi_do=pd.read_table(mypath+'Smart_RNA/dataset/PH_Species/MGI_DO.rpt')
	human_gene_disease=[]
	for index in mgi_do.index:
		if (mgi_do.loc[index,"Common Organism Name"]=="human"):
			disease=mgi_do.loc[index,"OMIM IDs"]
			sub_dis=disease.split("|")
			human_gene=mgi_do.loc[index,"EntrezGene ID"]
			#print(disease,sub_dis, human_gene)
			if disease:
				if human_gene:
					for dis in sub_dis:
						human_gene_disease.append([int(human_gene),dis.strip()])	
	mouse_gene_disease=dict()
	for index in mgi_do.index:
		mouse_gene=mgi_do.loc[index,"Mouse MGI ID"]
		if str(mouse_gene)!="nan":
			disease=str(mgi_do.loc[index,"OMIM IDs"]).strip()
			sub_dis=disease.split("|")
			for dis in sub_dis:
				if dis !="nan":
					try:
						mouse_gene_disease[mouse_gene].append(dis.strip())
					except:
						mouse_gene_disease[mouse_gene]=[dis.strip()]
	## generate the dictionary from human to mouse and from mouse to human
	human_to_mouse=dict()
	mouse_to_human=dict()
	geneName_to_id=dict()
	with open(mypath+'Smart_RNA/dataset/PH_Species/HMD_HumanPhenotype.rpt','r') as f:
		for line in f.readlines():
			data=line.split("\t")

			gene_name=data[0].strip()
			
			human_id=data[1].strip()
			mouse_id=data[5].strip()
			human_to_mouse[human_id]=mouse_id
			mouse_to_human[mouse_id]=human_id
			geneName_to_id[gene_name]=human_id
		# obtain the human gene and gene function association
		gene_go_feature=dict()
		prot_go_feature=dict()
		gene_expression_name=set()
		i=0
	with open(mypath+'Smart_RNA/dataset/PH_Species/goa_human.gaf',"r") as f:
		for line in f.readlines():
			data=line.split("\t")
			gene_name=data[2].strip()
			evidence_score=data[6].strip()
			go_id=data[4].strip()
			prot_id=data[10].strip()
			if gene_name in geneName_to_id.keys():
				if not ((evidence_score=="IEA") or (evidence_score=="ND")):
					human_gene=geneName_to_id[gene_name]
					if (human_gene in mouse_to_human.values()):
						try:
							gene_go_feature[human_gene].append(go_id)
							prot_go_feature[human_gene].append(prot_id)
			
						except:
							gene_go_feature[human_gene]=[go_id]
							prot_go_feature[human_gene]=[prot_id]

						#print(prot_go_feature[i])
			i=i+1
		with open(mypath+'Smart_RNA/dataset/PH_Species/gen_prot.pkl', 'wb') as f:
			pkl.dump(prot_go_feature, f)

        

	mp_pheno=pd.read_table(mypath+'Smart_RNA/dataset/PH_Species/MGI_GenePheno.rpt',names=["Allelic_Composition","Allele_Symbol","Allele_ID","Genetic_Background","Mammalian_Phenotype_ID","PubMed_ID","MGI_Marker","MGI_Genotype_Accession"])
	gene_mp_feature=dict()
	for index in mp_pheno.index:
		mouse_gene=mp_pheno.loc[index,"MGI_Marker"].strip()
		mouse_pheno=mp_pheno.loc[index,"Mammalian_Phenotype_ID"].strip()
		if (mouse_gene in mouse_to_human.keys()):
			human_gene=str(int(mouse_to_human[mouse_gene]))
			try:
		    		gene_mp_feature[human_gene].append(mouse_pheno)
			except:
			    gene_mp_feature[human_gene]=[mouse_pheno]
		# generate the intersection data
		gene_intersection_feature=dict()
	##  updated one
	gene_mp_intersection_feature=dict()
	gene_go_intersection_feature=dict()
	##  updated one
	for data in gene_go_feature.keys():
		if data in gene_mp_feature.keys():
			features=set()
			go_intersection_features=set()
			mp_intersection_features=set()
			for value in gene_go_feature[data]:
				features.add(value)
			for value in gene_mp_feature[data]:
				features.add(value)
			for value in gene_go_feature[data]:
				go_intersection_features.add(value)
			for value in gene_mp_feature[data]:
				mp_intersection_features.add(value)
			gene_go_intersection_feature[data]=go_intersection_features
			gene_mp_intersection_feature[data]=mp_intersection_features
			gene_intersection_feature[data]=features
	gene_total=set()
	for data in gene_go_feature.keys():
		gene_total.add(data)
	for data in gene_mp_feature.keys():
		gene_total.add(data)
	gene_union_feature=dict()
	for data in gene_total:
		features = set()
		if data in gene_go_feature.keys():
			for value in gene_go_feature[data]:
				features.add(value)
		if data in gene_mp_feature.keys():
			for value in gene_mp_feature[data]:
				features.add(value)
	gene_union_feature[data] = features

	print("num of mp  ",str(len(gene_mp_feature)))
	print("num of go  ",str(len(gene_go_feature)))
	print("num of intersection    ",str(len(gene_intersection_feature)))
	print("number of go in intersection",str(len(gene_go_intersection_feature)))
	print("number of mp in intersection",str(len(gene_mp_intersection_feature)))
	print("num of union     ",str(len(gene_union_feature)))




#generate the gene and disease associations based on mouse gene and mouse diseases
#associate the features with gene and diseases
#write the gene and disease out according to the specified data format
	path=mypath+'Smart_RNA/dataset/PH_Species/'
	f = open(path+"mp_association.txt","w")
	f.write( str(gene_mp_feature) )
	f.close()
	f = open(path+"go_association.txt","w")
	f.write( str(gene_go_feature) )
	f.close()
	f = open(path+"intersection_association.txt","w")
	f.write( str(gene_intersection_feature) )
	f.close()
	f = open(path+"go_intersection_association.txt","w")
	f.write( str(gene_go_intersection_feature) )
	f.close()
	f = open(path+"mp_intersection_association.txt","w")
	f.write( str(gene_mp_intersection_feature) )
	f.close()
	f = open(path+"union_association.txt","w")
	f.write( str(gene_union_feature) )
	f.close()
	generate_features(mouse_gene_disease,gene_mp_feature,dis_phe,"mp")
	generate_features(mouse_gene_disease,gene_go_feature,dis_phe,"go")
	generate_features(mouse_gene_disease,gene_intersection_feature,dis_phe,"intersection")
	generate_features(mouse_gene_disease,gene_go_intersection_feature, dis_phe,"go_intersection")
	generate_features(mouse_gene_disease, gene_mp_intersection_feature,dis_phe,"mp_intersection")
	generate_features(mouse_gene_disease,gene_union_feature,dis_phe,"union")


def generate_features(mouse_gene_disease_association,gene_feature,disease_feature,type):
	file=open(mypath+'Smart_RNA/dataset/PH_Species/'+type+'_association.txt',"w")
	disease_gene=dict()
	count=0
	total_gene=set()
	for gene in gene_feature.keys():
		total_gene.add(gene)
		features=gene_feature[gene]
		for feature in features:
			feature=feature.replace(":","_")
			file.write(gene+" "+'<http://purl.obolibrary.org/obo/'+str(feature)+">"+"\n")
	for disease in disease_feature.keys():
		features=disease_feature[disease]
		for feature in features:
			feature = feature.replace(":", "_")
			file.write(disease+" "+'<http://purl.obolibrary.org/obo/'+str(feature)+">"+"\n")
	for key in mouse_gene_disease_association.keys():
		try:
			human_gene=mouse_to_human[key]
			diseases=mouse_gene_disease_association[key]
			for dis in diseases:
				if dis in disease_feature.keys():
					if human_gene in gene_feature.keys():
						count+=1
                        # d_features=disease_feature[dis]
                        # for data in d_features:
                        #     file.write(dis+' '+"<http://purl.obolibrary.org/obo/"+str(data)+">"+"\n")
                        # for data in g_features:
                        #     file.write(human_gene+" "+'<http://purl.obolibrary.org/obo/'+str(data)+">"+"\n")
						try:
							disease_gene[dis].append(human_gene)
						except:
							disease_gene[dis]=[human_gene]
		except:
			pass
	print(disease_gene)
	with open(mypath+'Smart_RNA/dataset/'+type+'_disease_gene.pkl',"wb") as f:
		pkl.dump(disease_gene,f)
	with open(mypath+'Smart_RNA/dataset/'+type+'_gene_set.pkl',"wb") as f:
		pkl.dump(total_gene,f)
	print("the num of gene disease association :",str(count))
	print(" the number of total gene",len(total_gene))
	file.close()



##############################################################
#	VIRUS PHENOTYPES - PATHOGEN ASSOCIATION		     #
#							     #
##############################################################
def Create_Virus_Pheno_protein_asso():
	temp_filtered_data=''
	count=0
	i=0
	k=0
	data = {'Pathogene': [],'Phenotype': [],}
	df = pd.DataFrame(data)
	print(df)
	print('Create Virus Pheno/protein association')
	with open(mypath+'Smart_RNA/dataset/PH_Virus/patho_pheno_withsymbols.nt','r') as f, open(mypath+'Smart_RNA/dataset/PH_Virus/patho_pheno_withsymbols_filt.nt','w') as ff :
		for line in f.readlines():
			data=line.split(" ")
			Pathogene_ID=data[1].strip()
			Phenotype_Pat=data[2].strip()
			if "http://purl.obolibrary.org/obo/" in Pathogene_ID:
				new_line=Phenotype_Pat+ " " + Pathogene_ID+"\n"
				ff.write(new_line)								
		# obtain the human gene and gene function association




##############################################################
#			GRAPH CREATION			     #
#							     #
##############################################################
def set_graph_virus():
	global outfile
	global path_sentence	
	files_asso = askopenfilenames(parent=root,initialdir=mypath+'Smart_RNA/dataset/PH_Virus/'
	                           ,title='Choose Virus association files',multiple = 1)
	print(files_asso[0])

	files_axiom = askopenfilenames(parent=root,initialdir=mypath+'Smart_RNA/dataset/PH_Virus/'
	                           ,title='Choose Virus background knowledge files',multiple = 1)	
	print(files_axiom[0])
	outfile=mypath+'Smart_RNA/dataset/PH_Virus/W2V_Virus'
	path_sentence=mypath+'Smart_RNA/dataset/PH_Virus/walks_virus.txt'
	#f = open(path_sentence,"w")
	#f.write('')
	#f.close()
	generate_graph(files_asso[0],files_axiom[0])
	

def set_graph_species():
	global outfile
	global path_sentence
	files_asso = askopenfilenames(parent=root,initialdir=mypath+'Smart_RNA/dataset/PH_Species/'
	                           ,title='Choose Virus association files',multiple = 1)
	print(files_asso)
	files_axiom = askopenfilenames(parent=root,initialdir=mypath+'Smart_RNA/dataset/PH_Species/'
	                           ,title='Choose Virus background knowledge files',multiple = 1)
	print(files_axiom)
	outfile=mypath+'Smart_RNA/dataset/PH_Species/W2V_Species'
	path_sentence=mypath+'Smart_RNA/dataset/PH_Species/walks_species.txt'
	#f = open(path_sentence,"w")
	#f.write('')
	#f.close()
	generate_graph(files_asso[0],files_axiom[0])
	df = pd.read_csv(files_asso[0],  sep=' ')
	df.columns = ['a', 'b']
	df_filtered=df["a"].drop_duplicates()
	df_filtered.to_csv(mypath+'Smart_RNA/dataset/PH_Species/filter_genes_list.csv')


def Train_W2V_V_S():
	print('Train W2V Virus/Species Phenotypes')

	filenames = [mypath+'Smart_RNA/dataset/PH_Virus/walks_virus.txt',mypath+'Smart_RNA/dataset/PH_Species/walks_species.txt']
	file1=mypath+'Smart_RNA/dataset/PH_Virus/walks_virus.txt'
	file2=mypath+'Smart_RNA/dataset/PH_Species/walks_species.txt'
	file3=mypath+'Smart_RNA/walks_concat.txt'
	command='cat '+ file1 + ' ' +file2+' '+  '> '+ file3		
	os.system(command)
	print("start to train the word2vec models")
	sentences=gensim.models.word2vec.LineSentence(mypath+'/Smart_RNA/walks_concat.txt')
	model=gensim.models.Word2Vec(sentences,sg=1, min_count=1, size=100, window=10,iter=30,workers=5)
	w1 = "<http://purl.obolibrary.org/obo/RO_0002200>"
	print("Most similar to {0}".format(w1), model.wv.most_similar(positive=w1))
	print(
	"Most similar to {0}".format(w1),
		model.wv.most_similar(
		positive=w1,
		topn=15))
	print("Similarity between '<http://purl.obolibrary.org/obo/RO_0002200>' and '<http://purl.obolibrary.org/obo/RO_0002556>'",model.wv.similarity(w1="<http://purl.obolibrary.org/obo/RO_0002200>", w2="<http://purl.obolibrary.org/obo/RO_0002556>"))
	print('word vector <http://purl.obolibrary.org/obo/RO_0002200>:', model["<http://purl.obolibrary.org/obo/RO_0002200>"])

def test_virus():
	global prot_sim_result
	a=dict()
	b=dict()


	print('Train W2V Virus/Species Phenotypes')
	filenames = [mypath+'Smart_RNA/dataset/PH_Virus/walks_virus.txt',mypath+'Smart_RNA/dataset/PH_Species/walks_species.txt']
	file1=mypath+'Smart_RNA/dataset/PH_Virus/walks_virus.txt'
	file2=mypath+'Smart_RNA/dataset/PH_Species/walks_species.txt'
	file3=mypath+'Smart_RNA/walks_concat.txt'
	command='cat '+ file1 + ' ' +file2+' '+  '> '+ file3		
	os.system(command)
	print("start to train the word2vec models")
	sentences=gensim.models.word2vec.LineSentence(mypath+'/Smart_RNA/walks_concat.txt')
	model=gensim.models.Word2Vec(sentences,sg=1, min_count=1, size=100, window=10,iter=5,workers=5)
	df = pd.read_csv(mypath+'Smart_RNA/dataset/PH_Species/filter_genes_list.csv',  sep=',')
	idx = df.shape[0]
	with open(mypath+'Smart_RNA/dataset/PH_Species/gen_prot.pkl', 'rb') as f:
		prot_gene=pkl.load(f)		
	previous=''
	for k in range(0,idx):
		gene_under_test=df.loc[k, 'a']
		if gene_under_test in model.wv.vocab:		
			result=model.wv.similarity(w1="<http://purl.obolibrary.org/obo/RO_0002200>",w2=gene_under_test)
			if result>0.6:
				#print('Similarity between <http://purl.obolibrary.org/obo/RO_0002200> ' +   gene_under_test)
				#print(result)
				#text_report.insert(INSERT, "Similarity between <http://purl.obolibrary.org/obo/RO_0002200> " +   gene_under_test + ": " + str(result)+ '\n')
				if gene_under_test in prot_gene.keys():
					#print(prot_gene[gene_under_test][0])

					text_report.insert(INSERT, "protein associated: " +   str(prot_gene[gene_under_test][0]) + '\n')
					for key in prot_sim_result:
						if prot_sim_result[key] in prot_gene[gene_under_test]:
							print('**********POTENTIAL MATCH**********')
							text_report.insert(INSERT, "Potential Match Protein: " +   str(prot_gene[gene_under_test]) + '\n')
										


				#else:
				#	text_report.insert(INSERT, "No protein associated"'\n')
#	else:
#			print('gene not in vocab')
		

class Stack(object):
    def __init__(self):
        self.items = []

    def is_empty(self):
        return self.items == []

    def peak(self):
        return self.items[len(self.items) - 1]

    def size(self):
        return len(self.items)

    def push(self, item):
        self.items.append(item)
    def pop(self):
        return self.items.pop()

def split_sentence(sentence):
    new_sentence=""
    for d in sentence:
        if d=="(":
            new_sentence+="( "
        elif d==")":
            new_sentence+=" )"
        else:
            new_sentence+=d
    return new_sentence


def test_stack(sentence):
    sentence = split_sentence(sentence)
    sentence =sentence.split(" ")
    stack = Stack()
    for data in sentence:
        if data=="":
            pass
        else:
            stack.push(data)
    print(stack.items)

    while stack.is_empty()== False:

        print(stack.pop())

def judge_or_and_inside(data):
    if ("and" in data) or ("or" in data):
        return True
    else:
        return False

def judge_restrictions(data):
    restrctions = ["exactly","min","max","some","only"]
    for d in data:
        if d in restrctions:
            return True
    return False

def judge_condition(stack):
    if stack.is_empty()==False:
        if stack.peak() !="(":

            return True
        else:
            return False


    else:
        return False

def convert_triple(data):
    restrictions = ["some","only","exactly","min","max"]

    result = []
    first_entity = data[0]
    stack = Stack()

    tag=False
    for entity in data[2:]:
        if entity ==")":
            tag=True
            temp_data =[]
            while stack.peak() !="(":
                temp_data.append(stack.pop())
            stack.pop()
            if judge_restrictions(temp_data):
                new_relation = temp_data[-1]
                tail_node = temp_data[0]
                new_node = new_relation+" "+tail_node
                stack.push(new_node)
            else:
                if (not stack.is_empty()):
                    if stack.peak() =="(" or stack.peak() =="and" or stack.peak()=="or":


                        for da in temp_data:
                            stack.push(da)


                    else:
                        new_temp_data=[]
                        while judge_condition(stack):
                            new_temp_data.append(stack.pop())
                        if new_temp_data!=[]:
                            new_relation = new_temp_data[-1]
                            for da in temp_data:
                                if da !="and" and da !="or":
                                    new_element = new_relation+" "+da
                                    stack.push(new_element)
                            # for da in new_temp_data:
                            #     if da !="and" and da !="or":
                            #         new_element = new_relation+" "+da
                            #         stack.push(new_element)
                        else:

                            for da in temp_data:
                                stack.push(da)

        else:
            stack.push(entity)

    if tag:
        final_axioms = []
        while(stack.is_empty()==False):
            final_axioms.append(stack.pop())

        for element in final_axioms:
            if (element in restrictions):
                end_node = final_axioms[0]
                new_relation = final_axioms[-1]
                axiom = new_relation+" "+end_node
                final_axioms=[axiom]
                break

        for axiom in final_axioms:
            if axiom!="and" and axiom !="or":
                axiom=first_entity+" "+axiom


                axiom=axiom.split(" ")
                result.append(axiom)


        return result
    else:
        final_axioms = []
        while(stack.is_empty()==False):
            axiom =stack.pop()
            if axiom !="and" and axiom!="or":
                final_axioms.append(axiom)


        end_node = final_axioms[0]
        new_relation=final_axioms[-1]
        axiom = new_relation+" "+end_node
        axiom=first_entity+" "+axiom


        axiom=axiom.split(" ")

        result.append(axiom)


        return result

def convert_graph(data):
    sentence = split_sentence(data)
    sentence = sentence.split(" ")
    if len(sentence) <3:
        pass
    elif len(sentence)==3:
        result =[[sentence[0], sentence[1], sentence[2]]]
        new_result=[]
        for da in result:
            new_result.append([da[0]," ".join(da[1:-1]),da[-1]])

        return new_result
    else:

        result = convert_triple(sentence)
        new_result=[]
        for da in result:
            new_result.append([da[0]," ".join(da[1:-1]),da[-1]])


        return new_result
def generate_graph(annotation,axiom_file):
	global outfile
	G = nx.Graph()


	with open(axiom_file, "r") as f:
		for line in f.readlines():
			result = convert_graph(line.strip())
			if result!=None: 
				for entities in result:
					G.add_edge(entities[0].strip(), entities[2].strip())
					G.edges[entities[0].strip(), entities[2].strip()]["type"] = entities[1].strip()
					G.nodes[entities[0].strip()]["val"] = False
					G.nodes[entities[2].strip()]["val"] = False
	with open(annotation, "r") as f:
		for line in f.readlines():
			entities = line.split()
			G.add_edge(entities[0].strip(), entities[1].strip())
			G.edges[entities[0].strip(), entities[1].strip()]["type"] = "HasAssociation"
			G.nodes[entities[0].strip()]["val"] = False
			G.nodes[entities[1].strip()]["val"] = False

	gene_node_vector(G,annotation,outfile)
	#return G


lock = Lock()

WALK_LEN=30
N_WALKS=100

global data_pairs
data_pairs = []



def run_random_walks(G, nodes, num_walks=N_WALKS):
	print("now we start random walk")
	pairs = []
	for count, node in enumerate(nodes):
		if count<20:
			print(count, node)
			if G.degree(node) == 0:
				continue
			for i in range(num_walks):
				curr_node = node
				walk_accumulate=[]
				for j in range(WALK_LEN):
					next_node = random.choice(list(G.neighbors(curr_node)))
					type_nodes = G.edges[curr_node, next_node]["type"]
					if curr_node ==node:
						walk_accumulate.append(curr_node)
					walk_accumulate.append(type_nodes)
					walk_accumulate.append(next_node)
					curr_node = next_node
				pairs.append(walk_accumulate)
			if count % 1000 == 0:
				print("Done walks for", count, "nodes")
	write_file(pairs)


def run_walk(nodes,G):
    global data_pairs

    number=5
    length = len(nodes) // number

    processes = [mp.Process(target=run_random_walks, args=(G, nodes[(index) * length:(index + 1) * length])) for index
                 in range(number-1)]
    processes.append(mp.Process(target=run_random_walks, args=(G, nodes[(number-1) * length:len(nodes) - 1])))

    for p in processes:
        p.start()
    for p in processes:
        p.join()
    print("finish the work here")


def write_file(pair):
	global path_sentence	
	with lock:
		with open(path_sentence, "a") as fp:
			for p in pair:
				for sub_p in p:
					fp.write(str(sub_p)+" ")
				fp.write("\n")


def gene_node_vector(graph, entity_list,outfile):
	nodes_set=set()
	with open(entity_list,"r") as f:
		for line in f.readlines():
			data = line.strip().split()
			print('data', data)			
			for da in data:
				nodes_set.add(da)
	
	nodes_G= [n for n in graph.nodes()]
	G = graph.subgraph(nodes_G)
	nodes= [n for n in nodes_set]
	run_walk(nodes,G)

#annotation="/home/brignon/Documents/Smart_RNA/dataset/mp_association.txt"
#axiom_file="/home/brignon/Documents/Smart_RNA/axiomsorig.lst"
#G = generate_graph(annotation,axiom_file)
#gene_node_vector(G,annotation,outfile)



def NM_Ready():
	global NM_OK
	time.sleep(0.005)	
	soc.send(b'NM_status')	
	ready = select.select([soc], [], [], 1)
	data=''	
	if ready[0]:
		data = soc.recv(8092)
	np_eth = np.fromstring(data, dtype=np.ubyte)
	if (np_eth.size==1):
		
		if (np_eth[0]==1):
			NM_OK=1							
		else:
			NM_OK=0			

def Set_Reco_by_feature():
	if var_radio.get()==1: #NM RECO SOUND
		path_F=mypath+'Recorder_acq/Recorder_Folder/04_Feature_Sound/'
		Reco_by_feature_file(path_F)
		
	if var_radio.get()==2: #NM RECO VIB
		path_F=mypath+'Recorder_acq/Recorder_Folder/03_Feature_Vibration/'
		Reco_by_feature_file(path_F)


def Reco_by_feature():
	global NM_OK
	files = askopenfilenames(parent=root,initialdir=mypath+'Recorder_acq/Recorder_Folder/03_Feature_Vibration/'
                           ,title='Choose files',multiple = 1)
	files_for_test = root.tk.splitlist(files)
	print ("list of feature files =",files_for_test)	
	feature_file_for_test_number_of_raw=0	
	tot_raw=0
	raw=np.zeros(len(files_for_test))
	print(len(files_for_test))
	for i in range (len(files_for_test)):
		feature_file=np.loadtxt(files_for_test[i],delimiter=',',skiprows=2)
		for j in range(feature_file.shape[0]):
			if feature_file[j,0]!=0:
				tot_raw=tot_raw+1		
		feature_file_test_number_of_raw=int(tot_raw)
		feature_file_test_number_of_col=feature_file.shape[1]
	feature_file_test=np.zeros((int(tot_raw),int(feature_file_test_number_of_col)))
	prec=0	
	for i in range (len(files_for_test)):
		feature_file=np.loadtxt(files_for_test[i],delimiter=',',skiprows=2)			
		raw=0		
		for j in range(feature_file.shape[0]):
			if feature_file[j,0]!=0:
				feature_file_test[prec,:]=feature_file[j,:]
				prec=prec+1			
	data_inputs=feature_file_test
	feature_line=np.zeros((1, data_inputs.shape[1]))


def Reco_by_feature_file(file_name):
	global NM_OK
	files_for_test = file_name
	print ("list of feature files =",files_for_test)	
	feature_file_for_test_number_of_raw=0	
	tot_raw=0
	#raw=np.zeros(len(files_for_test))
	feature_file=np.loadtxt(files_for_test,delimiter=',')
	feature_file=Normalize(feature_file)
		
	for j in range(feature_file.shape[0]):
		if feature_file[j,0]!=0:
			tot_raw=tot_raw+1		
	feature_file_test_number_of_raw=int(tot_raw)
	feature_file_test_number_of_col=feature_file.shape[1]
	feature_file_test=np.zeros((int(tot_raw),int(feature_file_test_number_of_col)))
	prec=0	
	#feature_file=np.loadtxt(files_for_test,delimiter=',')			
	raw=0		
	for j in range(feature_file.shape[0]):
		if feature_file[j,0]!=0:
			feature_file_test[prec,:]=feature_file[j,:]
			prec=prec+1			
	data_inputs=feature_file_test
	feature_line=np.zeros((1, data_inputs.shape[1]))
	reco_file=np.zeros((int(tot_raw)+1,int(feature_file_test_number_of_col)), dtype=np.dtype('U50',1))
	for i in range(data_inputs.shape[0]):
		u=0				
		feature_line[0,0:265]=data_inputs[i,0:265] #Feature set
		#feature_line[0,0]=data_inputs[i,0]   #Index	
		a=np.zeros((4*feature_line.size),np.dtype(object))
		comb= b''	
		Send_Feature_Vector=b'Reco_Feature'
		soc.send(Send_Feature_Vector)
		k=0
		while k<266:				
			FV=struct.pack('i', int(feature_line[0,k]))	
			comb=comb+FV
			a[0]=comb 	
			k=k+1	
		soc.send(a[0])
		NM_OK=0
		count=0		
		while (NM_OK!=1):
			count=count+1					
			if count<1000:
				time.sleep(0.005)	
				soc.send(b'NM_status')	
				ready = select.select([soc], [], [], 1)
				data=''	
				if ready[0]:
					data = soc.recv(8092)
				np_eth = np.fromstring(data, dtype=np.ubyte)
				
				if (np_eth.size==31):
				
					if (np_eth[0]==1):
						NM_OK=1			
					else:
						NM_OK=0	
			else:
				NM_OK=1
		print('NSR ', np_eth[2])
		print('CONT ', np_eth[1])
		cat=(np_eth[3]<<8) + np_eth[4]		
		print('CAT ', (np_eth[3]<<8) + np_eth[4])
		print('DIST ', (np_eth[5]<<8) + np_eth[6])
		set_decode_result(cat)	


		#	reco_file[i+1,0]=data_inputs[i,0] 
	#	reco_file[i+1,2]=np_eth[3] #CAT
	#	reco_file[i+1,3]=np_eth[1] #CONT
	#	reco_file[i+1,4]=np_eth[2] #NSR
	#	reco_file[i+1,5]=(np_eth[5]<<8) + np_eth[6] #DIST
	#reco_file[0,:]=leg_F_Sound_Reco[0,:]
	#reco_file[1:(int(tot_raw)+1),10:]=data_inputs[:,10:]
	#time_stamp= str(datetime.now())
	#print(reco_file) 	
	#title = path + 'reco_file_' + time_stamp										
	#np.savetxt(title, reco_file,delimiter=',',fmt='%s')	
					
def set_decode_result(CAT):
	category=0	
	message=''
	if CAT==1:
		message='US HUMAN PERIOD 1'
		print('US HUMAN PERIOD 1')
	if CAT==2:
		message='US HUMAN PERIOD 2'
		print('US HUMAN PERIOD 2')
	if CAT==3:		
		message='US HUMAN PERIOD 3'
		print('US HUMAN PERIOD 3')
	if CAT==4:
		message='US HUMAN PERIOD 4'
		print('US HUMAN PERIOD 4')
	if CAT==5:
		message='EU HUMAN PERIOD 1'
		print('EU HUMAN PERIOD 1')
	if CAT==6:
		message='EU HUMAN PERIOD 2'
		print('EU HUMAN PERIOD 2')
	if CAT==7:
		message='EU HUMAN PERIOD 3'
		print('EU HUMAN PERIOD 3')
	if CAT==8:
		message='EU HUMAN PERIOD 4'
		print('EU HUMAN PERIOD 4')
	if CAT==9:
		message='CHINA HUMAN PERIOD 1'
		print('CHINA HUMAN PERIOD 1')
	if CAT==10:
		message='CHINA HUMAN PERIOD 2'
		print('CHINA HUMAN PERIOD 2')
	if CAT==11:
		message='CHINA HUMAN PERIOD 3'
		print('CHINA HUMAN PERIOD 3')
	if CAT==12:
		message='CHINA HUMAN PERIOD 4'
		print('CHINA HUMAN PERIOD 4')
	if CAT==13:
		message='SOUTH AMERICA PERIOD 1'
		print('SOUTH AMERICA HUMAN PERIOD 1')
	if CAT==14:
		message='SOUTH AMERICA PERIOD 2'
		print('SOUTH AMERICA HUMAN PERIOD 2')
	if CAT==15:
		message='SOUTH AMERICA PERIOD 3'
		print('SOUTH AMERICA HUMAN PERIOD 3')
	if CAT==16:
		message='SOUTH AMERICA PERIOD 4'
		print('SOUTH AMERICA HUMAN PERIOD 4')
	if CAT==17:
		message='ASIA HUMAN PERIOD 1'
		print('ASIA HUMAN PERIOD 1')
	if CAT==18:
		message='ASIA HUMAN PERIOD 2'
		print('ASIA HUMAN PERIOD 2')
	if CAT==19:
		message='ASIA HUMAN PERIOD 3'
		print('ASIA HUMAN PERIOD 3')
	if CAT==20:
		message='ASIA HUMAN PERIOD 4'
		print('ASIA HUMAN PERIOD 4')
	if CAT==21:
		message='AF HUMAN PERIOD 1'
		print('AF HUMAN PERIOD 1')

	if CAT==22:
		message='AF HUMAN PERIOD 2'
		print('AF HUMAN PERIOD 2')

	if CAT==23:
		message='AF HUMAN PERIOD 3'
		print('AF HUMAN PERIOD 3')

	if CAT==24:
		message='AF HUMAN PERIOD 4'
		print('AF HUMAN PERIOD 4')	
	if CAT==25:
		message='INDIA HUMAN PERIOD 1'
		print('INDIA HUMAN PERIOD 1')
	if CAT==26:
		message='INDIA HUMAN PERIOD 2'
		print('INDIA HUMAN PERIOD 2')
	if CAT==27:
		message='INDIA HUMAN PERIOD 3'
		print('INDIA HUMAN PERIOD 3')
	if CAT==28:
		message='INDIA HUMAN PERIOD 4'
		print('INDIA HUMAN PERIOD 4')





	if CAT==100:
		message='ANIMAL'		
		print('ANIMAL')
	text_report.insert(INSERT, "Sequence Recognized as  "+ message + '\n')
	#return category 						

def NeuroMEM_config(Channel_Status, MAXIF_value, MINIF_value, knowledge_loaded):
	s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	s.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
	s.connect((HOST, PORT))
	s.send(b'INIT CHANNEL STATUS, MAXIF AND MINIF')
	time.sleep(1)	
	
	Channel_Status=int(Channel_Status)
	
	Channel_Status=struct.pack('i', Channel_Status)
	MAXIF_value=int(MAXIF_value)
	MAXIF_value=struct.pack('i', MAXIF_value)
	MINIF_value=int(MINIF_value)
	MINIF_value=struct.pack('i', MINIF_value)
	
	s.send(Channel_Status)
	time.sleep(0.1)	
	s.send(MAXIF_value)
	time.sleep(0.1)
	s.send(MINIF_value)
	s.setblocking(0)
	ready = select.select([s], [], [], 1)
	if ready[0]:
		data = s.recv(8092)
		np_eth = np.fromstring(data, dtype=np.ubyte)
	recorder_ready = np_eth[0]		
	knowledge_size =np_eth[1]
	if knowledge_size ==0:
		knowledge_loaded=0							
	if knowledge_size !=0:
		knowledge_loaded=1			
	return(knowledge_loaded)
		
	s.close()
	time.sleep(0.1)



def knowledge_status(socket, knowledge_loaded):
	knowledge_size=np.zeros(2)
	global knowledge_slice	
	soc.send(b'knowledge_status')	
	ready = select.select([soc], [], [], 1)
	data=''	
	if ready[0]:
		data = soc.recv(8092)
	np_eth = np.fromstring(data, dtype=np.ubyte)
	if (np_eth.size==2):
		print((np_eth[0]<<8) + np_eth[1], knowledge_slice)

		if ((np_eth[0]<<8) + np_eth[1])!=knowledge_slice:
			knowledge_loaded=0							
		else:
			knowledge_loaded=1			
	else:
		knowledge_loaded=0
	return(knowledge_loaded)




#set save_KB_knowledge()
def restore_knowledge_KB():
	global knowledge_from_simu
	global KB_context_req
	global KNOWLEDGE_EXPERT_MATRIX
	print('Restore Knowledge from the recorder')
	Knowledge_Restore = tkinter.simpledialog.askstring("Restore a Knowledge from the Recorder: ","Enter name of the Knowledge:")
	Knowledge_Restore_name = mypath + '/00_Knowledge/' + Knowledge_Restore+'.txt'
	Knowledge_expert_name = Knowledge_Restore_name + '_expert.txt'
	np.savetxt(Knowledge_Restore_name, knowledge_from_simu, fmt='%i', delimiter=',')
	print(KNOWLEDGE_EXPERT_MATRIX)	
		 
	if KB_context_req==1:	
		ROI_H=KB_ROI_H_set.get()
		ROI_W=KB_ROI_W_set.get()
		KNOWLEDGE_EXPERT_MATRIX[KB_context,0]=1
		KNOWLEDGE_EXPERT_MATRIX[KB_context,1]=ROI_W
		KNOWLEDGE_EXPERT_MATRIX[KB_context,2]=ROI_H
		KNOWLEDGE_EXPERT_MATRIX[KB_context,3]=1	
		print(KNOWLEDGE_EXPERT_MATRIX)
	np.savetxt(Knowledge_expert_name, KNOWLEDGE_EXPERT_MATRIX, fmt='%i', delimiter=',')	


##########################################################################################################
#				LEARNING ALGO RCE-NN							 #
##########################################################################################################
def RCE_NN_Learning():
	global ANNOTATION_MATRIX
	global feature_number	
	global neuron_max	
	global knowledge_from_simu
	inc1=5
	inc2=10
	alpha=0.9
	beta=0.8
	teta=2
	Rmax=1500
	Rmin=2#3000
	Rw=Rmax	
	index_neuron=0
	neuron_activated=0 

	Neuron=np.zeros((neuron_max, 260), dtype=np.int16) #CAT, CONT, DIST, FEATURES
	R=np.zeros((neuron_max), dtype=np.int16) #CAT, CONT, DIST, FEATURES
	Temp_R=np.zeros((neuron_max), dtype=np.int16) 
	f=np.zeros((neuron_max), dtype=np.uint16) #CAT, CONT, DIST, FEATURE
	CAT_neuron=np.zeros((neuron_max), dtype=np.int16)
	CONT_neuron=np.zeros((neuron_max), dtype=np.int16)
	u=np.zeros((260), dtype=np.uint16) #CAT, CONT, DIST, FEATURES
		
	k=0	
	somme=0
	k=0
	for l in range(ANNOTATION_MATRIX.shape[0]):
		if ANNOTATION_MATRIX[l,0]!=0:
			k=k+1	
	input_vector_number=k
	Input_Vector=np.zeros((input_vector_number, 260), dtype=np.uint16) #CAT, CONT, FEATURES
	DIST_x_c=np.zeros((input_vector_number, neuron_max), dtype=np.uint16) 
	CAT_Input_Vector=np.zeros((input_vector_number), dtype=np.uint16)
	CONT_Input_Vector=np.zeros((input_vector_number), dtype=np.int16)
	k=0	
	max_context=0
	for l in range(ANNOTATION_MATRIX.shape[0]):
		if ANNOTATION_MATRIX[l,0]!=0:
			
			for j in range(256):		
				Input_Vector[k,j]=ANNOTATION_MATRIX[l,j+10]
				
			CAT_Input_Vector[k]=ANNOTATION_MATRIX[l,5]
			CONT_Input_Vector[k]=ANNOTATION_MATRIX[l,6]
			
			if CONT_Input_Vector[k]>max_context:
				max_context=CONT_Input_Vector[k]
			k=k+1
	for i in range (neuron_max):	
		R[i]=Rmax
	increment=0
	print(max_context, input_vector_number, neuron_max, feature_number)
	while increment<1:
		for m in range(1, max_context+1,1):	
			print('**************** CONTEXT:', m)			
			for l in range (input_vector_number):	
				#print(m, CONT_Input_Vector[l])				
				if CONT_Input_Vector[l]==m:				
					Temp_index_neuron=0			
					neuron_activated=0
					for i in range (neuron_max):
						Neuron_OK=0
						for j in range (feature_number):
							if Neuron[i,j]!=0:
								Neuron_OK=1
						somme=0				
						for j in range (feature_number):
						#u[j]=(Input_Vector[l,j]-Neuron[i,j])*(Input_Vector[l,j]-Neuron[i,j])						
							u[j]=abs(Input_Vector[l,j]-Neuron[i,j])
							somme=somme+u[j]			
						DIST_x_c[l,i]=somme
						#DIST_x_c[l,i]=int(math.sqrt(somme))

						if DIST_x_c[l,i]==0:
							neuron_activated=1	
						if Neuron_OK==1 and DIST_x_c[l,i]<=R[i] and DIST_x_c[l,i]!=0 :
							neuron_activated=1
							if CAT_Input_Vector[l]==CAT_neuron[i]:
								if (DIST_x_c[l,i]>(R[i]+Rmin)/2):
									f[i]=f[i]+1
									print('inc 0: ',f[i]) 
								if (DIST_x_c[l,i]>Rmin) and (DIST_x_c[l,i]<(R[i]+Rmin)/2):
									f[i]=f[i]+inc1
									print('inc 1: ',f[i])	
								if (DIST_x_c[l,i]<=Rmin):
									f[i]=f[i]+inc2
									print('inc 2: ',f[i])
							print('AIF before', i, CAT_Input_Vector[l], R[i])						
							if CAT_Input_Vector[l]!=CAT_neuron[i]:
								R[i]=DIST_x_c[l,i]							
						#		print('different_cat... shrink')
								if (CAT_Input_Vector[l]!=0): #Do not recreate a vector in case of cat0
									Neuron[index_neuron, 0:feature_number]=Input_Vector[l, 0:feature_number]
									f[index_neuron]=0
									Temp_R[Temp_index_neuron]=DIST_x_c[l,i]
						#		print('New Neuron temp',Temp_R[Temp_index_neuron], CAT_neuron[Temp_index_neuron])
									Temp_index_neuron=Temp_index_neuron+1
							if CAT_Input_Vector[l]==0:
								print('AIF after', i, R[i])						
											
										
					if Temp_index_neuron!=0:
						#if Temp_R.min()!=0:
						Min_val=10000000				
						for g in range(Temp_index_neuron):
							if Temp_R[g]<Min_val:
								if Temp_R[g]!=0:
									Min_val=Temp_R[g]			
						R[index_neuron]=Min_val
						CONT_neuron[index_neuron]=CONT_Input_Vector[l]
						CAT_neuron[index_neuron]=CAT_Input_Vector[l]
						#print('new neuron with ',R[index_neuron], CAT_Input_Vector[l])
								
						index_neuron=index_neuron+1
						Temp_index_neuron=0
						Temp_R[:]=0	
					if (neuron_activated==0): # and DIST_x_c[l,i]!=0:
						if (CAT_Input_Vector[l]!=0): #Do not recreate a vector in case of cat0
						
							CAT_neuron[index_neuron]=CAT_Input_Vector[l]
							CONT_neuron[index_neuron]=CONT_Input_Vector[l]
							Neuron[index_neuron, 0:feature_number]=Input_Vector[l, 0:feature_number]
							f[index_neuron]=0
							R[index_neuron]=Rmax
							#print('New Neuron',l, index_neuron, R[index_neuron], CAT_neuron[index_neuron])
							index_neuron=index_neuron+1
										
					#print('CAT...')
				#print(l, CAT_neuron[l])
			#Rw=alpha*Rw
		increment=increment+1
	#	print(Neuron[0:index_neuron,0:feature_number])
			#Neuron_knowledge=np.unique(Neuron,axis=0)
			#print(Neuron_knowledge)	
	#	print(R[0:index_neuron])
	print('index_neuron', index_neuron)
	text_report.insert(INSERT, "Knowledge with "+ str(index_neuron) + " neuron(s) created!"+ '\n')
	knowledge_from_simu[0]=5892 #HEX 1704
	knowledge_from_simu[2]=index_neuron
	for i in range(neuron_max):
		knowledge_from_simu[4+i*260]=CONT_neuron[i]
		#knowledge_from_simu[4+1+i*260:4+1+feature_number+i*260]=Neuron[i,0:feature_number]
		knowledge_from_simu[4+1+i*260:4+1+feature_number+i*260]=Neuron[i,0:feature_number]
		knowledge_from_simu[4+1+i*260+feature_number]=R[i]
		knowledge_from_simu[4+2+i*260+feature_number]=Rmin
		knowledge_from_simu[4+3+i*260+feature_number]=CAT_neuron[i]
	#print(knowledge_from_simu)
	size_knowledge=16+4*index_neuron*(1+256+3)
		
	#print('size_knowledge=',size_knowledge)
	knowledge=np.zeros((size_knowledge),dtype=np.dtype(object))
	knowledge=knowledge_from_simu[0:4+4+i*261+feature_number+1]

	#print(knowledge)
	########################################################	
	#            SEND KNOWLEDGE FILE TO ETHERNET           #
	########################################################

	#s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	#s.connect((HOST, PORT))
	print('Initialisation of the knowledge........')
	knowledge_init=b'AABBCCDDEEFF'
	soc.send(knowledge_init)
	time.sleep(0.1)
	comb= b''	
	j=0	
	a=np.zeros((size_knowledge),dtype=np.dtype(object))	
	while j<4:
		comb=comb+struct.pack("I", knowledge[j])
		a[0]=comb 	
		j=j+1	
	
	soc.send(a[0])		
	u=0			
	print('Loading the knowledge to the Platform.........')
	time.sleep(0.1)
	comb= b''
	knowledge_slice=(size_knowledge-16)/1040	
	while u<knowledge_slice:
		off=4
		j=0		
		comb= b''		
		a[0]=struct.pack("I", knowledge[off+u*260])
		soc.send(a[0])		
		time.sleep(0.01)						
		comb= b''		
		j=0						
		while j<256:
			comb=comb+struct.pack("I", knowledge[off+1+j+u*260])
			a[0]=comb 	
			j=j+1
		soc.send(a[0])		
		time.sleep(0.01)		
		comb= b''		
		j=0								
		while j<3:
			comb=comb+struct.pack("I", knowledge[off+1+256+j+u*260])
			a[0]=comb 	
			j=j+1
		#print(comb)
		soc.send(a[0])		
		time.sleep(0.01)				
		u=u+1
		print('neuron loading', u)		
		time.sleep(0.5)

def set_track_mutation():
	frame_track_mutation.grid(row=0, column=0)
	text_report.grid(row=0, column=0)
	frame_Virus_Host.grid_remove()

def set_virus_host_detection():
	frame_track_mutation.grid_remove()
	frame_Virus_Host.grid(row=0, column=0)	
	text_report.grid(row=0, column=0)

root.wm_title("SMART RNA")
menu = Menu(root)
root.config(menu=menu)
file = Menu(menu)

file.add_command(label = 'Track Mutation', command = set_track_mutation)
file.add_command(label = 'Virus Host Detectoon', command = set_virus_host_detection)

file.add_command(label = 'Exit', command = lambda:exit())
mode = Menu(menu)
menu.add_cascade(label = 'File', menu = file)
frame16 = Frame(root)
root.image=plt.imread(mypath+'Recorder_acq/Recorder_Folder/'+'main.jpg')
fig = plt.figure(figsize=(8,6))
plt.axis('off')
im = plt.imshow(root.image)



frame_track_mutation= Frame(root)
frame_text_report= Frame(root)

frame_track_mutation.grid(row=0, column=0)

Button_REF_RNA=Button(frame_track_mutation,text='Select reference RNA', command=Open_REF_RNA_FILE,height = 2, width = 25)
Button_query_RNA=Button(frame_track_mutation,text='Select query RNA', command=Open_TEST_RNA_FILE,height = 2, width = 25)

Button_Raw_stat=Button(frame_track_mutation,text='Statistics', command=occurence_chart,height = 2, width = 25)
Button_sequence_similarity=Button(frame_track_mutation,text='Sequence Similarty', command=set_chain_similarity,height = 2, width = 25)

Button_DATA_BASE=Button(frame_track_mutation,text='Database', command=OpenUrl,height = 2, width = 25)
Button_FEATURE_TRAIN=Button(frame_track_mutation,text='CREATE FEATURE FOR TRAINING', command=Set_RNA_AI_TRAINING,height = 2, width = 25)
Button_FEATURE_TEST=Button(frame_track_mutation,text='CREATE FEATURE FOR TESTING', command=Set_RNA_AI_RECO,height = 2, width = 25)
Button_Train_AI_file=Button(frame_track_mutation,text='AI TRAINING', command=Set_RNA_AI_TRAINING_FILE,height = 2, width = 25)
Button_RECO_AI_file=Button(frame_track_mutation,text='AI RECOGNITION', command=Set_RNA_AI_RECO_FILE,height = 2, width = 25)






Label_Single_request = Label(frame_track_mutation, text="Single Request")
Label_Mutation_Tracking = Label(frame_track_mutation, text="Mutation Tracking")

		


Button_DATA_BASE.grid(row =0, column =1)
Button_REF_RNA.grid(row =0, column =2)
Label_Single_request.grid(row =1, column = 1, sticky =W)
Button_query_RNA.grid(row =2, column =1)
Button_sequence_similarity.grid(row =2, column =2)
#Button_Raw_stat.grid(row =2, column =3)


Label_Mutation_Tracking.grid(row =3, column = 1, sticky =W)
Button_FEATURE_TRAIN.grid(row =4, column =1)
Button_FEATURE_TEST.grid(row =4, column =2)
Button_Train_AI_file.grid(row =5, column =1)
Button_RECO_AI_file.grid(row =5, column =2)


text_report = Text(frame_text_report)
Scroll=Scrollbar(frame_text_report)


text_report.insert(INSERT, "WELCOME TO THE SMART RNA TOOLSET"+ '\n')
frame_text_report.grid(row=0, column=1)
text_report.grid(row=0, column=0)
Scroll.grid(row=0, column=1)
text_report 
text_report.config(yscrollcommand=Scroll.set)
Scroll.config(command=text_report.yview) 

frame_Virus_Host= Frame(root)
Label_Corpus_Creation = Label(frame_Virus_Host, text="Creation of the corpus")

Label_PPI = Label(frame_Virus_Host, text="Virus Host detection PPI")

Button_Virus_PH_Axioms=Button(frame_Virus_Host,text='Virus: Phenotypes Background Knowledge', command=Create_Virus_Pheno_FILE,height = 2, width = 35)
Button_Species_PH_Axioms=Button(frame_Virus_Host,text='Species: Background Knowledge', command=Create_Species_Pheno_FILE,height = 2, width = 35)
Button_Virus_PH_Prot_asso=Button(frame_Virus_Host,text='Virus: Association files', command=Create_Virus_Pheno_protein_asso,height = 2, width = 35)
Button_Species_PH_Prot_asso=Button(frame_Virus_Host,text='Species: Association files', command=Create_Species_Pheno_protein_asso,height = 2, width = 35)
Button_Virus_Graph=Button(frame_Virus_Host,text='Virus: Graph', command=set_graph_virus,height = 2, width = 35)
Button_Species_Graph=Button(frame_Virus_Host,text='species: Graph', command=set_graph_species,height = 2, width = 35)
Button_V_S_Word2Vect=Button(frame_Virus_Host,text='Train Word2Vec Virus/Species', command=Train_W2V_V_S,height = 2, width = 35)
Button_test_virus=Button(frame_Virus_Host,text='Virus to be tested', command=test_virus,height = 2, width = 25)
Button_set_protein_comp=Button(frame_Virus_Host,text='Protein Comparison', command=set_protein_comp,height = 2, width = 25)
Similarity_progressbar=ttk.Progressbar(frame_Virus_Host,orient="horizontal",length=300,mode="determinate")
#text_report = Text(frame_Smart_RNA_Control)
#frame_Smart_RNA_Control.grid(row=0, column=2)


Label_Corpus_Creation.grid(row =0, column =0)
Button_Species_PH_Axioms.grid(row =1, column =0)
Button_Virus_PH_Axioms.grid(row =1, column =1)
Button_Species_PH_Prot_asso.grid(row =2, column =0)
Button_Virus_PH_Prot_asso.grid(row =2, column =1)
Button_Species_Graph.grid(row =3, column =0)
Button_Virus_Graph.grid(row =3, column =1)
Button_V_S_Word2Vect.grid(row =4, column =0)

Label_PPI.grid(row =5, column =0)
Button_set_protein_comp.grid(row =6, column =0)
Button_test_virus.grid(row =6, column =1)
root.poll = True
root.mainloop()



