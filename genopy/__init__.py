# Copyright (c) 2019 Yann BOUYERON
#
#
# licensed under GNU GPL version 3 (or later)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 




import os

import matplotlib

if 'DISPLAY' not in os.environ:
	
	matplotlib.use('Agg') 

import numpy as np
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_nucleotide, generic_protein, SingleLetterAlphabet
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC, seq1, seq3
from Bio import SeqIO
from io import StringIO
from Bio import AlignIO
from Bio import Phylo
from Bio import Restriction as zym

from Bio.Phylo.TreeConstruction import *
from Bio.Phylo.Consensus import *

import fnmatch
import pandas as pd
import seaborn as sns
import unicodedata
import ipfshttpclient

from Bio.Emboss.Applications import NeedleCommandline
from Bio.Emboss.Applications import WaterCommandline
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

import os

__all__ = ['toSeq', 'matrix', 'matrix2ipfs', 'matrix2png', 'transcribe', 'translate', 'retro_transcribe', 'codegen', 'show', 'reademboss', 'needle', 'water', 'tree_nj', 'tree_upgma', 'parsimony_tree', 'showanag', 'showgenix', 'edi2fasta', 'search', 'creat', 'mkfas', 'mkfasx', 'clustal', 'phylo', 'draw_tree'] 

__version__ = "0.0.2"

__author__ = "Yann Bouyeron"

ipfs = ipfshttpclient.connect('/dns/ipfs.infura.io/tcp/5001/https')


wd = os.getcwd()



def recurliste(path):
	
	"""
	argument: path est une chaine de caracteres
	return: une liste recursive des dossiers et fichiers contenus dans path.
	
	c'est une alternative a glob.glob('**',recursive=True) qui ne fonctionne pas en python < 3.5"""

	r = []
	for root, dir, files in os.walk(path):
		if root != path:
			r.append(root.replace(path+'/', ''))
		for f in fnmatch.filter(files, "*"):
			x = root+'/'+f
			r.append(x)
	return r



	
		
	
	



###########################################   Alphabets et langages 	###########


def canbedna(seq):
	
	"""
	Teste si une séquence (str ou Seq ou SeqRecord) peut être un dna
		
		argument: seq (str, Seq, ou SeqRecord): la séquence à tester
		return: bool
	"""

	if type(seq) == SeqRecord:
		
		seq = seq.seq
		
	elif type(seq) != str and type(seq) != Seq:
		
		raise ValueError("L'argument seq doit être de type str, Seq ou SeqRecord")


	div = list(set(seq))
			
	for i in div:
		
		if i not in list(IUPAC.unambiguous_dna.letters) and i.upper() not in list(IUPAC.unambiguous_dna.letters):
			
			return False
			
	return True
	


							
def canberna(seq):
	
	"""
	Teste si une séquence (str ou Seq ou SeqRecord) peut être un rna
	
		argument: seq (str, Seq, ou SeqRecord): la séquence à tester
		return: bool
	"""
		
	if type(seq) == SeqRecord:
		
		seq = seq.seq
		
	elif type(seq) != str and type(seq) != Seq:
		
		raise ValueError("L'argument seq doit être de type str, Seq ou SeqRecord")
	
	div = list(set(seq))
	
	for i in div:
		
		if i not in list(IUPAC.unambiguous_rna.letters) and i.upper() not in list(IUPAC.unambiguous_rna.letters):
			
			return False
			
	return True


def canbeexprot(seq):
	
	"""
	Teste si une séquence (str ou Seq ou SeqRecord) peut être une extended proteine
	
		argument: seq (str, Seq, ou SeqRecord): la séquence à tester
		return: bool
	"""
	
	if type(seq) == SeqRecord:
		
		seq = seq.seq
		
	elif type(seq) != str and type(seq) != Seq:
		
		raise ValueError("L'argument seq doit être de type str, Seq ou SeqRecord")
		
	
	div = list(set(seq))
	
	for i in div:
		
		if i not in list(IUPAC.extended_protein.letters) + ['*']:
			
			return False
			
	return True


def canbeprot(seq):
	
	"""
	Teste si une séquence (str ou Seq ou SeqRecord) peut être une proteine
	
		argument: seq (str, Seq, ou SeqRecord): la séquence à tester
		return: bool
	"""
	
	if type(seq) == SeqRecord:
		
		seq = seq.seq
		
	elif type(seq) != str and type(seq) != Seq:
		
		raise ValueError("L'argument seq doit être de type str, Seq ou SeqRecord")
	
	div = list(set(seq))
	
	for i in div:
		
		if i not in list(IUPAC.protein.letters) + ['*']:
			
			return False
			
	return True






def testalpha(seq):
	
	"""
	Teste l'alphabet d'une séquence
	
	argument: Seq ou SeqRecord  , la séquence à tester
	
	return:
		generic_nucleotide si il y'a un doute entre adn ou arn
		IUPAC.unambigous_dna
		IUPAC.unambigous_rna
		IUPAC.protein
		IUPAC.extended_protein
		ou None si la séquence n'est pas reconnue
		
	"""
	
	if canbedna(seq) and canberna(seq):
	
		return generic_nucleotide
	
	elif canbedna(seq):
		
		return IUPAC.unambiguous_dna
		
	elif canberna(seq):
		
		return IUPAC.unambiguous_rna
		
	
	elif canbeprot(seq):
		
		return IUPAC.protein
		
	elif canbeexprot(seq):
		
		return IUPAC.extended_protein
		
	else:
		
		return None



def toSeq(seq):
	
	"""Teste si une séquence est de type str Seq ou SeqRecord. Return Seq object"""

	if type(seq) == SeqRecord:
		
		return seq.seq
		
	elif type(seq) != str and type(seq) != Seq:
		
		raise ValueError("L'argument seq doit être de type str, Seq ou SeqRecord")
	
	elif type(seq) == str:
		
		return Seq(seq, testalpha(seq))

	else: 
		
		return seq




def langage(seq):
	
	"""
	Determine le mode d'ecriture à 1 ou 3 lettres des séquences peptidiques
	
	arguments:
		
		seq: la séquence à tester (str) ou (objet Seq) ou (objet SeqRecord)
	
	return:
		1 (int) si langage à une lettre
		3 (int) si langage à trois lettres
		ou None si la séquence n'est pas reconnue comme une séquence proteique
	"""
	
	seq = toSeq(seq)
	
	seq = str(seq)

	alpha = testalpha(seq)

	if seq.isupper() and (alpha == IUPAC.protein or alpha == IUPAC.extended_protein):

		return 1

	elif seq.isupper() == False and seq.islower() == False and (testalpha(seq1(seq)) == IUPAC.protein or testalpha(seq1(seq)) == IUPAC.extended_protein):

		return 3

	else:

		return None








def matrix(align, model="identity", frame = False):
	
	"""
	Construction d'une matrice de similitudes à partir de l'alignement
	
		arguments: 
                    align: objet align retourné par clustal
                    model: voir biopython models
                    frame: bool
                
		return: matrice de distances (défaut: si frame=False), DataFrame (si frame=True)
	"""
	
	calc = DistanceCalculator(model)
	dm = calc.get_distance(align)

	if frame == True:
		
		df = pd.DataFrame(dm.matrix)
		df.columns = dm.names
		df.index = dm.names
		
		for i in range(len(dm.names)):
			
			for c in range(len(dm.names)):
				
				if df.isnull().iloc[c,i] == True:
					
					df.iloc[c,i] = df.iloc[i,c]
		
		return df
		
	elif frame == False:
		
		return dm



def matrix2ipfs(matrix):	
	
	"""transforme une matrice de distance en html publié sur IPFS, retourne le hash ipfs
	
	La page html est accessible à l'adresse https://ipfs.io/ipfs/<hash>
	
	Pour plus d'informations sur ipfs: https://gist.github.com/YannBouyeron/53e6d67782dcff5995754b0a7333fa0b
	
	"""
	
	if type(matrix) == DistanceMatrix:
		
		try:
			
			df = pd.DataFrame(matrix.matrix)
			df.columns = matrix.names
			df.index = matrix.names
	
		except:
			
			return False
			
	elif type(matrix) == pd.DataFrame:
		
		df = matrix
			

	html = df.to_html()
	
	html = html.replace("\n", "")
	
	r = ipfs.add_str(html)
	
	return r

	
	
def matrix2png(matrix, path = "", scale=0.6, cmap="YlGnBu"):
	
	"""
	Transforme une matrice de distance en image png
	
	arguments:
		matrix: pd.DataFrame ou distance matrice retournée par la fonction matrix
		path (otpionnel): path en .png de l'image
		scale (optionnel) [defaut=0.6] : taille du texte (légendes, titres, cellules)
		cmap (optionnel) [defaut="YlGnBu"] : palette de couleurs
		
	return: objet heatmap
	"""
	
	if type(matrix) == DistanceMatrix:
		
		try:
			
			df = pd.DataFrame(matrix.matrix)
			df.columns = matrix.names
			df.index = matrix.names
	
		except:
			
			return False
			
	elif type(matrix) == pd.DataFrame:
		
		df = matrix	
		
		
	sns.set(font_scale = scale)
	
	plt.close()
	
	m = sns.heatmap(df, annot=True, cmap=cmap)
	
	if 'DISPLAY' in os.environ:
		
		plt.show()
		
	if path != "":
	
		plt.savefig(path)
	
	return m
		
	

	

			
def transcribe(*seq, out = False):
	
	""" 
        Transcription d'un ou plusieurs dna SeqRecord. 
        
        arguments:
            *seq: liste de sequences Seq ou SeqRecord
            out: bool. Si out == True, chaque rna est sauvegardé dans un fasta portant son nom de séquence en .fas 
            
        return: rna SeqRecord ou une liste de rna SeqRecord
        
        """
	
	result = []
	
	for i in seq:
		
		if type(i) == SeqRecord:
			
			if testalpha(i) == IUPAC.unambiguous_dna or testalpha(i) == generic_nucleotide:
			
				dna = i.seq
				
				try:
					rna = dna.transcribe()
				
					sr = creat(i.name + '.rna', rna, 'transcription de [{0}]'.format(i.description), out=out)
				
					result.append(sr)
				
				except:
					pass
	
	if len(result) == 1:
		
		return result[0]
		
	else:
		
		return result
			

def translate(*seq, table_id = 1, to_stop = True, stop_symbol = " ", cds = True, out = False):	
	
	""" 
        Traduction d'un ou plusieurs dna ou rna SeqRecord. 
        
        arguments:
            *seq: liste de dna ou rna Seq ou SeqRecord
            table_id: indice de la table du code génétique à utiliser (par defaut: 1 c'est le code standard)
            to_stop: bool  True par défaut: la traduction s'arrete au codon stop
            stop_symbol: représentation de l'arret de la traduction
            cds: coding sequence: bool  True par défaut: la traduction commence au condon initiateur et se termine au codon stop
            out: bool False par défaut. Si out == True, chaque sequence peptidique est sauvegardée dans un fasta. 
            
        return: proteine SeqRecord ou une liste de proteines SeqRecord
        
        """
	
	result = []
	
	if table_id not in list(range(32))[1:]:
		
		table_id = 1
	
	for i in seq:
		
		if type(i) == SeqRecord:
			
			if testalpha(i) == IUPAC.unambiguous_rna or testalpha(i) == generic_nucleotide or testalpha(i) == IUPAC.unambiguous_dna:
			
				nucl = i.seq
				
				try:
					
					prot = nucl.translate(table = table_id, to_stop = to_stop, stop_symbol = stop_symbol, cds = cds)
					
					sr = creat(i.name + '.prot', prot, 'traduction de [{0}]'.format(i.description), out=out)
				
					result.append(sr)
				
				except:
					pass
	
	if len(result) == 1:
		
		return result[0]
		
	else:
		
		return result
		

def retro_transcribe(*seq, out = False):
	
	""" 
        Retro transcription d'un ou plusieurs rna SeqRecord. 
        
        arguments:
            *seq: liste de rna Seq ou SeqRecord
            out: bool False par défaut. Si out == True, chaque dna est sauvegardé dans un fasta. 
            
        return: dna SeqRecord ou une liste de dna SeqRecord
        
        """
	
	result = []
	
	for i in seq:
		
		if type(i) == SeqRecord:
			
			if testalpha(i) == IUPAC.unambiguous_rna or testalpha(i) == generic_nucleotide:
			
				rna = i.seq
				
				try:
					dna = rna.back_transcribe()
				
					sr = creat(i.name + '.dnac', dna, 'retro transcription de [{0}]'.format(i.description), out=out)
				
					result.append(sr)
				
				except:
					pass
	
	if len(result) == 1:
		
		return result[0]
		
	else:
		
		return result
			

def codegen(table_id = None, bydna = False):
	
	"""Affiche et retourne la table du code génétique correspondant à l'id entré en argument
	Si aucun id entré en argument: retourne menu de choix d'id
	"""
	
	d = CodonTable.unambiguous_dna_by_id
	r = CodonTable.unambiguous_rna_by_id
	
	print('\n')
	print('\n')
	
	while table_id not in list(range(32))[1:]:
		
		for i in d:
				
			print(d[i].id, ':', d[i].names[0])
			
		print('\n')
		print('\n')
		
		table_id = int(input("Entrez le numéro de la table à afficher: "))
			
	
	if bydna == True:
		
		print('\n')
		
		print('DNA table: \n')
		
		print(d[table_id])
		
		print('\n')
		
		return d[table_id]
	
	else:
		
		print('\n')
		
		print('RNA table: \n')
		
		print(r[table_id])	
		
		print('\n')
		
		return r[table_id]



################################### Représentation des séquences ####################################
			
def show(seq, start = 0, stop = None, width = None, peprep = 3):
	
	
	"""
        Affiche une représentation d'une séquence nucléotique ou peptidique avec règle graduée
	
        arguments:
            seq: Seq ou SeqRecord à représenter
            start et stop: limite de la représentation. Les arguments start et stop doivent être des entiers (int); 
            width: largeur de la représentation. width doit être un multiple de 10. 
            peprep: 1 ou 3: représentation des acides aminés. Par défaut, la représentation des seq peptidique (peprep) et de 3 lettres/aa
            
        """
	
	if testalpha(seq) == generic_nucleotide or testalpha(seq) == IUPAC.unambiguous_dna or testalpha(seq) == IUPAC.unambiguous_rna:
		
		return shown(seq, start=start, stop=stop, width=width)
	
	
	#elif testalpha(seq) == IUPAC.protein or testalpha(seq) == IUPAC.extended_protein:
	
	elif langage(seq) == 1 or langage(seq) == 3:
		
		if peprep == 3 and langage(seq) == 1:
			
			seq = toSeq(seq)
			seq = seq3(seq)
		
		return showp(seq, start=start, stop=stop, width=width)
	
	
def shown(seq, start = 0, stop = None, width = None):
	
	"""Présentation d'une séquence nucléotique découpée en tronçons avec règle graduée.
	Les arguments start et stop doivent être des entiers (int); width doit être un multiple de 10.
	
	"""
	
	seq = toSeq(seq)
	
	seq = str(seq)
	
	if width == None:
		width = 60

	if width/10 != int(width/10):
		raise ValueError("Le parametre width doit être un multiple de 10")

	if stop == None: 
		stop = len(seq)

	elif type(stop) == type(int()): 
		pass

	else: 
		raise ValueError("Le parametre stop doit être un entier")


	
	if type(start) != type(int()):
		
		raise ValueError("Le parametre start doit être un entier")

	d = start
	
	f = d + width
	
	

	txt = ''

	tir = "----:----|" * int(width/10)

	while d <= stop:
		
		if f > stop:
			
			tir = tir[:len(seq[d:stop])]
		
		
		
		num = "          "
		
		for k in np.arange(10, width + 10, 10):
				
			num = num + """{0}{1}""".format("         " if k == 10 else "        " if  d+k <= 100 else "       " if d+k <= 1000 else "      ", d + k if d+k <= stop else "")
		
		
		if f <= stop:
			
			txt = txt + "          " + seq[d:f] + "\n" + "          " + tir + "\n" + num + "\n\n"
		
		else:
			
			txt = txt + "          " + seq[d:stop] + "\n" + "          " + tir + "\n" + num + "\n\n"
		
		d, f = d + width, f + width

	print("\n")
	print(txt)
	print('\n')
	
	
	
	
def showp(seq, start = 0, stop = None, width = None):
	
	"""Présentation d'une séquence peptidique découpée en tronçons avec règle graduée.
	Les arguments start et stop doivent être des entiers (int); width doit être un multiple de 10.
	"""
	
	seq = toSeq(seq)
		
	seq = str(seq)
	
	if langage(seq) == 1:
	
		if width == None:
			width = 60

		if width/10 != int(width/10):
			raise ValueError("Le parametre width doit être un multiple de 10")
	
		if stop == None: 
			stop = len(seq)
	
		elif type(stop) == type(int()): 
			pass
	
		else: 
			raise ValueError("Le parametre stop doit être un entier")
	
	
		
		if type(start) != type(int()):
			raise ValueError("Le parametre start doit être un entier")
	
		
			
		d = start
		
		f = d + width
	
		txt = ''	
		
		
		
	
	
		tir = "----:----|" * int(width/10)
	
		while d <= stop:
			
			if f > stop:
				
				tir = tir[:len(seq[d:stop])]
			
			
			
			num = "          "
			
			for k in np.arange(10, width + 10, 10):
					
				num = num + """{0}{1}""".format("         " if k == 10 else "        " if  d+k <= 100 else "       " , d + k if d+k <= stop else "")
			
			
			if f <= stop:
				
				txt = txt + "          " + seq[d:f] + "\n" + "          " + tir + "\n" + num + "\n\n"
			
			else:
				
				txt = txt + "          " + seq[d:stop] + "\n" + "          " + tir + "\n" + num + "\n\n"
			
			d, f = d + width, f + width
	
	
	
		print("\n")
		print(txt)
		print('\n')
	
		
		
		
	elif langage(seq) == 3:
		
		if width == None:
			width = 30
		
		if width/10 != int(width/10):
			raise ValueError("Le le parametre width doit être un multiple de 10")
	
		if stop == None: 
			stop = int(len(seq)/3)
	
		elif type(stop) == type(int()): 
			pass
	
		else: 
			raise ValueError("Le parametre stop doit être un entier")
	
		if type(start) != type(int()):
			raise ValueError("Le parametre start doit être un entier")
	
		d = start
		
		f = d + width
	
		txt = ''	
		
		
		#decoupage de la séquence
		w = 0
		h = 3
		
		z = []
		
		while h <= len(seq):
			
			z.append(seq[w:h])
			
			w, h = h, h + 3
		
		
		tir = " -  -  -  -  :  -  -  -  -  | " * int(width/10)
	
		while d <= stop:
			
			if f > stop:
				
				tir = tir[:len(''.join(z[d:stop]))]
			
			
			
			num = "          "
			
			for k in np.arange(10, width + 10, 10):
					
				num = num + """{0}{1}""".format("                            " if  d+k <= 100 or k == 10 else "                           ", d + k if d+k <= stop else "")
			
			
			if f <= stop:
				
				txt = txt + "          " + ''.join(z[d:f]) + "\n" + "          " + tir + "\n" + num + "\n\n"
			
			else:
				
				txt = txt + "          " + ''.join(z[d:stop]) + "\n" + "          " + tir + "\n" + num + "\n\n"
			
			d, f = d + width, f + width
	
	
	
		print("\n")
		print(txt)
		print('\n')
		






############################### Représentation des alignements de séquences ############################


def reademboss(path = 'emb.aln'):
	
	"""Cette fonction affiche le contenu de l'alignement emboss contenu dans le fichier à l'adresse "path", et return un objet align"""
	
	try:
	
		align = AlignIO.read(path, "emboss")
		
	except:
		
		return None
		
	else:
		
		with open(path, 'rb') as nd:
			
			x = nd.read()
			x = x.decode()
		
		print("\n", x)
		
		return align
		


def showanag(align, start = 0, stop = None, width = 60):
	
	
	"""Présentation d'un alignement nucléotique découpé en tronçons avec règle graduée , format anagene-like.
	start et stop doivent être des entiers; width doit être un multiple de 10
	
	argument positionnel:
		align: l'alignement (str) retourné par la fonction AlignIO.read de biopython: exemple : align = AlignIO.read(path, "clustal")
	
	arguments facultatifs:
		start: (int) index du premier nucléotide représenté
		stop: (int) index du dernier nucléotide représenté
		width: (int) largeur des lignes; width doit être un multiple de 10
		
		Legendes: _ deletion, - similitude avec la référence, * similitude
	"""
	
	
	s = list(align)
	
	#sequence de reférence 
	sa = str(s[0].seq)
		
	#autres sequence de l'alignement
	ss = [str(i.seq) for i in s[1:]]


	

	#verification des options

	if width/10 != int(width/10):
		raise ValueError("Le parametre width doit être un multiple de 10")

	if stop == None: 
		stop = len(sa)

	elif type(stop) == type(int()) and stop <= len(sa): 
		pass

	else: 
		raise ValueError("Le parametre stop doit être un entier inférieur ou égal à la longueur de l'alignement")

	if type(start) != type(int()):
		raise ValueError("Le parametre start doit être un entier")

	
	#modification de sa: remplacement des - par des _
	san = ""
	
	for i in sa:
		if i == '-': san += "_"
		else: san += i

	#modif des sequences secondaires de ss, création des symboles *
	
	sss = []
	
	for sb in ss:
		
		sbn = ""
	
		for a , b in zip(sa,sb):
			
			if b == "-": 
				
				sbn += "_"
				
				
			elif b == a: 
				
				sbn += "-"
				
				
			elif b != a: 
				
				sbn += b
				
		sss.append(sbn)





	#récupération des noms des séquences
	
	#nom de la séquence de référence
	nsa = str(s[0].id)
	
	#liste des noms des autres séquences de l'alignement
	nnn = [str(i.id) for i in s[1:]]
	
	
	#recherche du nom de séquence le plus long et gestion des espaces entre nom et séquences
	lt = nnn
	lt.append(nsa)
	
	ltlen = []
	for i in lt:
		ltlen.append(len(i))
	lnm = max(ltlen) + 8     #c'est la longueur du nom le plus long plus 8 espaces
	
	
	
	#gestion des symboles * de similitudes
	symb = ""
	
	for i , j in enumerate(sa):
		
		count = 0
		
		for k in ss:
			
			if j != k[i]: count += 1
			
		if count == 0: symb += "*"
		else: symb += " "
	
	#gestion des tirets
	tir = "----:----|" * int(width/10)
	
	

	d = start
	
	f = d + width
	
	
	#initiation de la str devant contenir la repsrésentation de l'alignement
	txt = ""
	

	while d <= stop:
		
		if f > stop:
			
			tir = tir[:len(sa[d:stop])]
		
		
		
		num = " "*lnm
		
		for k in np.arange(10, width + 10, 10):
				
			num = num + """{0}{1}""".format("         " if k == 10 else "        " if  d+k <= 100 else "       " if d+k <= 1000 else "      " , d + k if d+k <= stop else "")
		
		
		if f <= stop:
			
			txt = txt + " "*lnm + symb[d:f] + "\n"
			txt = txt + nsa + " "*(lnm - len(nsa)) + san[d:f] + "\n"
			
			for nsb, sbn in zip(nnn,sss):
				txt = txt + nsb + " "*(lnm - len(nsb)) + sbn[d:f] + "\n"
				
			txt = txt + "\n"
			txt = txt + " "*lnm + tir + "\n" + num + "\n\n\n\n"
		
		else:
			
			txt = txt + " "*lnm + symb[d:stop] + "\n"
			txt = txt + nsa + " "*(lnm - len(nsa)) + san[d:stop] + "\n"
			
			for nsb, sbn in zip(nnn,sss):
				txt = txt + nsb + " "*(lnm - len(nsb)) + sbn[d:stop] + "\n"
				
			txt = txt + "\n"
			txt = txt + " "*lnm + tir + "\n" + num + "\n\n\n\n"
			
			
		
		d, f = d + width, f + width



	print("\n")
	print(txt)
	print('\n')
	
	





def showgenix(align, start = 0, stop = None, width = 60, aaa = True):
	
	
	"""Présentation d'un alignement nucléotique ou peptidique à une lettre découpé en tronçons avec règle graduée.
	start et stop doivent être des entiers; width doit être un multiple de 10
	
	argument positionnel:
		align: l'alignement (str) retourné par la fonction AlignIO.read de biopython: exemple : align = AlignIO.read(path, "clustal")
		
	arguments facultatifs:
		start: (int) index du premier nucléotide représenté
		stop: (int) index du dernier nucléotide représenté
		width: (int) largeur des lignes; width doit être un multiple de 10
		text: (str) déscription
		aaa: (bool) acides aminés à 3 lettres si True
	
	Legendes: - deletion, * similitude	
	"""
	
	s = list(align)
	
	#sequence de reférence 
	sa = str(s[0].seq)
		
	#autres sequence de l'alignement
	ss = [str(i.seq) for i in s[1:]]


	if aaa == True and (testalpha(sa) == IUPAC.protein or testalpha(sa) == IUPAC.extended_protein):
		
		return showgenixp3(align, start=start, stop=stop, stdout=stdout, text=text)

	

	#verification des options

	if width/10 != int(width/10):
		raise ValueError("Le parametre width doit être un multiple de 10")

	if stop == None: 
		stop = len(sa)

	elif type(stop) == type(int()) and stop <= len(sa): 
		pass

	else: 
		raise ValueError("Le parametre stop doit être un entier inférieur ou égal à la longueur de l'alignement")

	if type(start) != type(int()):
		raise ValueError("Le parametre start doit être un entier")


	#récupération des noms des séquences
	
	#nom de la séquence de référence
	nsa = str(s[0].id)
	
	#liste des noms des autres séquences de l'alignement
	nnn = [str(i.id) for i in s[1:]]
	
	
	#recherche du nom de séquence le plus long et gestion des espaces entre nom et séquences
	lt = nnn
	lt.append(nsa)
	
	ltlen = []
	for i in lt:
		ltlen.append(len(i))
	lnm = max(ltlen) + 8     #c'est la longueur du nom le plus long plus 8 espaces
	
	
	
	#gestion des symboles * de similitudes
	#symb = align._star_info ca ne marche pas toujours !!!
	symb = ""
	
	for i , j in enumerate(sa):
		
		count = 0
		
		for k in ss:
			
			if j != k[i]: count += 1
			
		if count == 0: symb += "*"
		else: symb += " "
	
	
	#gestion des tirets
	tir = "----:----|" * int(width/10)
	
	
	#initiation des intervals
	d = start
	
	f = d + width
	
	
	
	#initiation de la str devant contenir la repsrésentation de l'alignement
	txt = ""
	

	while d <= stop:
		
		if f > stop:
			
			tir = tir[:len(sa[d:stop])]
		
		
		
		num = " "*lnm
		
		for k in np.arange(10, width + 10, 10):
				
			num = num + """{0}{1}""".format("         " if k == 10 else "        " if  d+k <= 100 else "       " if d+k <= 1000 else "      ", d + k if d+k <= stop else "")
		
		
		if f <= stop:
			
			txt = txt + " "*lnm + symb[d:f] + "\n"
			txt = txt + nsa + " "*(lnm - len(nsa)) + sa[d:f] + "\n"
			
			for nsb, sb in zip(nnn,ss):
				txt = txt + nsb + " "*(lnm - len(nsb)) + sb[d:f] + "\n"
				
			txt = txt + "\n"
			txt = txt + " "*lnm + tir + "\n" + num + "\n\n\n\n"
		
		else:
			
			txt = txt + " "*lnm + symb[d:stop] + "\n"
			txt = txt + nsa + " "*(lnm - len(nsa)) + sa[d:stop] + "\n"
			
			for nsb, sb in zip(nnn,ss):
				txt = txt + nsb + " "*(lnm - len(nsb)) + sb[d:stop] + "\n"
				
			txt = txt + "\n"
			txt = txt + " "*lnm + tir + "\n" + num + "\n\n\n\n"
			
			
		
		d, f = d + width, f + width



	print("\n")
	print(txt)
	print('\n')
	
	




def treeliste(s):
	
	"""Conversion d'une str représentant une seq petpitidique à 3 lettres : MetValGlu en une liste d'acides aminés ['Met','Val','Glu']"""
	
	
	d = 0
	f = 3
	
	tl = []
	

	while f <= len(s):
		
		tl.append(s[d:f])
		d , f = f , f+3

	return tl




def showgenixp3(align, start = 0, stop = None, width = 30):
	
	
	"""Présentation d'un alignement peptidique avec langage 3 lettres decoupé en tronçons avec règle graduée.
	start et stop doivent être des entiers; width doit être un multiple de 10
	
	argument positionnel:
		align: l'alignement (str) retourné par la fonction AlignIO.read de biopython: exemple : align = AlignIO.read(path, "clustal")
		
	arguments facultatifs:
		start: (int) index du premier nucléotide représenté
		stop: (int) index du dernier nucléotide représenté
		width: (int) largeur des lignes; width doit être un multiple de 10
		
	Legendes: - deletion, * similitude
	"""
	
	
	s = list(align)
	
	#sequence de reférence 
	sa = str(s[0].seq)
		
	#autres sequences de l'alignement
	ss = [str(i.seq) for i in s[1:]]
	
	
	#verification des options

	if width/10 != int(width/10):
		raise ValueError("Le parametre width doit être un multiple de 10")

	if stop == None: 
		stop = len(sa)

	elif type(stop) == type(int()) and stop <= len(sa): 
		pass

	else: 
		raise ValueError("Le parametre stop doit être un entier inférieur ou égal à la longueur de l'alignement")

	if type(start) != type(int()):
		raise ValueError("Le parametre start doit être un entier")
	
		
	#passage au langage a 3 lettres
	sa = seq3(sa, custom_map={"*": ""}, undef_code=' - ')
	
	#transformation en list d'acides amines ['Met','Val','Glu']
	sa = treeliste(sa)
	
	#les deletions notées - dans l'alignement sont converties par seq3 en Xaa , il faut donc les corriger.... inutile avec undef_code=' - '
	#for i, j in enumerate(sa): 
		
		#if j == 'Xaa':
			
			#sa[i] = ' - '
	
	
	
	#creation d'une liste de liste des autres sequences avec langage a 3 lettres
	sss = []
	
	for i in ss:
		
		j = seq3(i, custom_map={"*": ""}, undef_code=' - ')
		
		ssl = treeliste(j)
		
		#for i, j in enumerate(ssl): 
		
			#if j == 'Xaa':
			
				#ssl[i] = ' - '
		
		sss.append(ssl)
		
	ss = sss
	
		
	#gestion des tirets
	tir = " -  -  -  -  :  -  -  -  -  | " * int(width/10)
	
	
	#récupération des noms des séquences
	
	#nom de la séquence de référence
	nsa = str(s[0].id)
	
	#liste des noms des autres séquences de l'alignement
	nnn = [str(i.id) for i in s[1:]]
	
	
	#recherche du nom de séquence le plus long et gestion des espaces entre nom et séquences
	lt = nnn
	lt.append(nsa)
	
	ltlen = []
	for i in lt:
		ltlen.append(len(i))
	lnm = max(ltlen) + 8     #c'est la longueur du nom le plus long plus 8 espaces
	
	
	
	#gestion des symboles * de similitudes
	symb = ""
	
	for i , j in enumerate(sa):
		
		count = 0
		
		for k in ss:
			
			if j != k[i]: 
				count += 1
			
		if count == 0: 
			symb += " * "
			
		else: 
			symb += "   "
	
	symb = treeliste(symb)
	
	
	
	#initiation des intervals
	d = start
	
	f = d + width
	

	#initiation de la str devant contenir la repsrésentation de l'alignement
	txt = ""
	
	while d <= stop:
		
		if f > stop:
			
			tir = tir[:len(''.join(sa[d:stop]))]
		
		
		num = " "*lnm
		
		for k in np.arange(10, width + 10, 10):
				
			num = num + """{0}{1}""".format("                            " if k == 10 else "                            " if  d+k <= 100 else "                           " , d + k if d+k <= stop else "")
		
		
		if f <= stop:
			
			txt = txt + " "*lnm + "".join(symb[d:f]) + "\n"
			txt = txt + nsa + " "*(lnm - len(nsa)) + "".join(sa[d:f]) + "\n"
			
			for nsb, sb in zip(nnn,ss):
				txt = txt + nsb + " "*(lnm - len(nsb)) + "".join(sb[d:f]) + "\n"
				
			txt = txt + "\n"
			txt = txt + " "*lnm + tir + "\n" + num + "\n\n\n\n"
		
		else:
			
			txt = txt + " "*lnm + "".join(symb[d:stop]) + "\n"
			txt = txt + nsa + " "*(lnm - len(nsa)) + "".join(sa[d:stop]) + "\n"
			
			for nsb, sb in zip(nnn,ss):
				txt = txt + nsb + " "*(lnm - len(nsb)) + "".join(sb[d:stop]) + "\n"
				
			txt = txt + "\n"
			txt = txt + " "*lnm + tir + "\n" + num + "\n\n\n\n"
			
					
		d, f = d + width, f + width

	print("\n")
	print(txt)
	print('\n')
	
	
	
	
####################################### Outils de conversions edi / fasta ###########################
	
def supraccent(x):
	
	y = unicodedata.normalize('NFKD', x).encode('ascii', 'ignore')
	
	return y.decode()
	
def edi2fasta(path_edi, path_fas):
	
	"""
	Transforme un fichier edi en fasta
	
	Arguments:
		
		path_edi: path du fichier en .edi
		path_fasta: path du fichier en .fas
		
		les extensions sont obligatoires !!!
		
	"""

	#verification des extensions
	if path_edi[len(path_edi)-4:] != '.edi':
		raise ValueError('path_edi recquiert .edi extension  !!')
		
		
	if path_fas[len(path_fas)-4:] != '.fas':
		#raise ValueError('path_fas recquiert .fas extension  !!')
		path_fas += ".fas" 
		
		
	#ouverture et lecture du fichier .edi
	
	with open(path_edi,'r', encoding='ISO-8859-1') as f:
		r = f.read()
		r = supraccent(r)
			
			
	

	#decoupage en liste de sequences (si plusieurs sequences dans le edi)
	l = r.split(';-')
	
	
	list_fasta = []
	
	for i in l[:len(l)-1]:
		
		#decoupage de la sequence
		a = i.split('\n;')
		
		try:
			name = a[1][1:] #recuperation du nom de la sequence et suppression de l'espace
			if '-' in name:
				name = name.replace('-', '_')
			
			
			if len(a) == 6:
				des = a[4][1:] #description de la sequence
			
			elif len(a) > 6:
				des = a[4:len(a)-1]
				des = "".join(des)
			else:
				des = ''
			
			seq = a[len(a)-1] #la seq est en dernier dans la liste
			name = name.replace(" ","_")
			
			
		except IndexError:
			return 'Le fichier .edi presente une structure anormale'
		

		else:
			#suppression des espaces et des \n de la sequence
			y = list(seq) #conversion str to list
			while ' ' in y:
				y.remove(' ')
			while '\n' in y:
				y.remove('\n')
				
			seq = "".join(y) #reconversion list to str
			
			fas = SeqRecord(Seq(seq), id = name , description=des)
			
			list_fasta.append(fas)
			
	
	#ecriture du fasta
	SeqIO.write(list_fasta, path_fas, "fasta")
			
	

			
						
			
############################################################################################################################






def search(*keys, path = wd):

	"""
	Recherche tous les fasta ayant dans leur nom l'un des mots clé recherché dans le repertoire et sous repertoires du path
	
	arguments:
		path: str : path du répertoire de recherche, défaut: répertoire courant
		keys: str : mots clés de recherche
	
	return:
		une liste des SeqRecord
	"""
	
	lso = []
	
	#si path ne se termine pas par un / on le rajoute de maniière à obtenir des path complets depuis la racine
	if path[len(path)-1] != '/':
		path = path + '/'
	
	#si la tuple de mot_cle est vide, on créé le mot_cle '' qui correspond à n'importe quel caractère
	if keys == ():
		keys = ('',)
	
	for i in recurliste(path):
	
		name = i.split('/')
		name = name[len(name)-1]
		
		if os.path.isfile(i) and '.fas' in name:
			for j in keys:
				if j in name:
					so = list(SeqIO.parse(i, "fasta"))
					lso = lso + so
					break  


	#elimination des doublons
	list_SeqRec_open = []
	
	list_name =[]
	
	for i in lso:
		
		if i.name not in list_name:
			list_SeqRec_open.append(i)
			list_name.append(i.name)

		
	return list_SeqRec_open


def creat(seq_name, seq, des = "", out = False):
	
	"""Création d'un SeqRecord et d'un fasta mono séquence
	
	arguments:
		seq_name: (str) nom de la séquence créée
		seq: (str) ou (Seq) ou (SeqRecord) la séquence
		des: (str) déscription du SeqRecord
		out (bool): si False (défaut) , le fasta n'est pas créé
		
	retrun: 
		le SeqRecord créé
	"""
	
	seq = toSeq(seq)
	
	#creation d'un objet SeqRecord
	
	try:
		
		f = SeqRecord(seq, id= seq_name, name = seq_name, description = des)
	
		x = f.format('fasta')
	
	except:
		
		return None
		
	
	if out == True:

		SeqIO.write(f, os.getcwd() + '/' + f.name + '.fas', "fasta")
		
	return f



def mkfas(*seq, path = wd):
	
	"""
	Cette fonction permet de créer des fasta mono séquence
	
	arguments: 
		*seq: str: liste de SeqRecord
		path: str: repertoire (absolu) dans lequel seront enregistrés les fasta créés. Le path par défaut est le repertoire courant

	return:
		cette fonction ne return rien, elle crée des fichiers fasta mono séquence
		les fichiers sont enregistrés dans le path indiqué en argument, ils portent leur nom de sequence en . fas
	"""
	
	if os.path.isdir(path):
		
		if path[len(path)-1] != '/':
			
			path = path + '/'
		
		for i in seq:
			
			try:
				SeqIO.write(i, path + i.name + '.fas', "fasta")
			
			except:
				pass


def mkfasx(out, *seq):
	
	"""
	Création de fasta multi séquence
	
	arguments: 

		*seq: liste de SeqRecord 
		out: str: path absolu du fasta mutli séquence qui sera créé
	
	return:
		
		cette fonction ne retourne rien, elle crée un fichier fasta multi séquences contenant chaque séquence passée en argument
	"""
	
	SeqIO.write(seq, out, "fasta")



def clustal(*seq, out = 'comp.aln', std = False):
	
	"""
	Alignement multiple ClustalW
	
	arguments:
		
		seq: liste des SeqRecord à aligner
		out: path du fichier contenant l'alignement créé
		std: bool, si True, return align et stdout , défaut = False
		
	Return: l'objet align
	"""
	
	mkfasx(out, *seq)

	cline = ClustalwCommandline("clustalw", infile=out, score='PERCENT')
	
	stdout, Stderr = cline()
	
	align = AlignIO.read(out, "clustal")
	
	if std == True:
		
		return align, stdout
	
	return align
	


	
	
def phylo(*seq, out = "comp.aln"):
	
	"""
	Méthode simplifiée de création d'un arbre phylogénétique à partir d'une liste de SeqRecord
	
		argument: *seq: liste de SeqRecord à aligner
		return: (align, stdout, tree)
	
	"""
	
	align, stdout = clustal(*seq, out = out, std = True)
	
	mpd =matrix(align, frame = True)
	
	print('\n')
	print(mpd)
	
	tree = Phylo.read("comp.dnd", "newick")
	
	print('\n')
	
	draw_tree(tree)
	
	return (align, stdout, tree)	
	

def tree_nj(dm):
	
	"""
	Création d'un arbre nj à partir de la matrice de distance retournée par la fonction matrix
	return: tree
	"""
	
	cst = DistanceTreeConstructor()
	
	tree = cst.nj(dm)

	draw_tree(tree)
	
	return tree
	

def tree_upgma(dm):
	
	"""
	Création d'un arbre upgma à partir de la matrice de distance retournée par la fonction matrix
	return: tree
	"""
	
	cst = DistanceTreeConstructor()
	
	tree = cst.upgma(dm)
	
	draw_tree(tree)
	
	return tree
	

def parsimony_tree(align, starting_tree):
	
	"""
	Création d'un arbre basé sur la méthode de parcimonie à partir de l'alignement et d'un arbre de départ (nj)

	arguments:
		align: objet align retourné par clustal()
		satrting_tree: objet tree retourné par tree_nj()
	
	return: tree
	"""
	
	scorer = ParsimonyScorer() #possibilité de passer une matrice de parsimony en argument ?!?
	
	searcher = NNITreeSearcher(scorer)
	
	cst = ParsimonyTreeConstructor(searcher, starting_tree)
	
	tree = cst.build_tree(align)
	
	draw_tree(tree)
	
	return tree
	
	

def draw_tree(tree, distance=False):
	
	if "DISPLAY" not in os.environ:
		
		print('\n')
		
		Phylo.draw_ascii(tree)
		print('\n')
		
	if distance == True:
		
		Phylo.draw(tree, branch_labels=lambda c: c.branch_length)
	
	elif distance == False:
		
		Phylo.draw(tree)


def needle(*seq, gapopen=10, gapextend=0.5, out ='emb.aln'):

	"""Alignement global par la méthode de Needleman
	
	arguments:
		
		seq: couple de 2 SeqRecord à aligner
		gapopen: pénalité de gap
		gapextend: pénalité d'expansion
		out: nom du fichier emboss créé
		
	return: un objet align
	"""
	
	
	mkfasx("seqa.fas", seq[0])
	
	mkfasx("seqb.fas", *seq[1:])
	
	needle_cline = NeedleCommandline(asequence='seqa.fas', bsequence='seqb.fas', gapopen=gapopen, gapextend=gapextend, outfile=out)

	stdout, stderr = needle_cline()

	os.remove('seqa.fas')
	os.remove('seqb.fas')

	if len(seq) < 3:
		align = AlignIO.read(out, "emboss")
		return align



def water(*seq, gapopen=10, gapextend=0.5, out ='emb.aln'):

	"""Alignement global par la méthode de Needleman
	
	arguments:
		
		seq: couple de 2 SeqRecord à aligner
		gapopen: pénalité de gap
		gapextend: pénalité d'expansion
		out: nom du fichier emboss créé
		
	return: un objet align
	"""
	
	
	mkfasx('seqa.fas', seq[0])
	
	mkfasx("seqb.fas", *seq[1:])
	
	water_cline = WaterCommandline(asequence='seqa.fas', bsequence='seqb.fas', gapopen=gapopen, gapextend=gapextend, outfile=out)

	stdout, stderr = water_cline()

	os.remove('seqa.fas')
	os.remove('seqb.fas')

	if len(seq) < 3:
		align = AlignIO.read(out, "emboss")
		return align




