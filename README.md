# Genopy


Analyses génétiques et phylogénétiques pour les SVT

<a name ="installgenoy"></a>

## Installation


### Installation des dépendances.

Genopy requiert python >= 3.4 et des dépendances python qui ne sont pas présentes dans la librairie standard (numpy, biopython, pandas, matplotlib, seaborn, ipfshttpclient) et des dépendances non python (Clustalw, Emboss, Rebase, Phylip).

#### Installation des dépendances Python.

	sudo pip install numpy biopython pandas matplotlib seaborn ipfshttpclient

Si vous avez plusieurs versions de python et que vous souhaitez installer Genopy pour la version de python 3.6 (par exemple):

	sudo python3.6 -m pip install numpy biopython pandas matplotlib seaborn ipfshttpclient

#### Installation de Clustalw

[http://www.clustal.org/clustal2/](http://www.clustal.org/clustal2/)

Sur Debian:

	sudo apt-get install clustalw

#### Installation de Emboss

[http://emboss.sourceforge.net/download/](http://emboss.sourceforge.net/download/)

Sur Debian:

	sudo apt-get install emboss


#### Installation rebase pour emboss. (Facultatif)

Rebase permet  d’utiliser certaines fonctions liées aux enzymes de restriction; il s'agit notamment des fonctions restrict et remap.

1: placez vous dans le répertoire de rebase

	cd /usr/share/EMBOSS/data/REBASE

2: télécharger les fichiers withrefm et proto depuis le serveur ftp de rebase

	sudo wget ftp://ftp.ebi.ac.uk/pub/databases/rebase/withrefm*

Puis 

	sudo wget ftp://ftp.ebi.ac.uk/pub/databases/rebase/proto*

3: décompresser 

	sudo uncompress withrefm.707.gz

Et 

	sudo uncompress proto.707.gz

Vous n'aurez pas forcément le version 707.... c'est à adapter

4: extraire rebase

	rebaseextract

Et voilà... c'est terminé.

#### Installation de Phylip

[http://evolution.genetics.washington.edu/phylip/install.html](http://evolution.genetics.washington.edu/phylip/install.html)


Sur Debian:

	sudo apt-get install phylip

### Installation de genopy

#### Avec pip:

    sudo pip install genopy
    
Ou, pour une version donnée de python (exemple 3.6):

    sudo python3.6 -m pip install genopy


#### Depuis github:

	git clone https://github.com/YannBouyeron/genopy
	
	cd genopy
	
	sudo python setup.py install


## Licence

Ce code est sous licence <a href="http://www.gnu.org/licenses">GPL3</a>


## Exemple


	>>> import genopy as gp
	>>>
	>>> list_seq = gp.search("primates")
	
	
	>>> list_seq
	[SeqRecord(seq=Seq('TTCTTTCATGGGGAAGCAGATTTGGGTGCCACCCAAGTATTGACTCACCCATCA...CAT', SingleLetterAlphabet()), id='Refhumaine.adn', name='Refhumaine.adn', description=' Refhumaine.adn', dbxrefs=[]), SeqRecord(seq=Seq('TTCTTTCATGGGGGATCAGATTTGGGTACCACCCCAGTACCGACCCATTTCCAG...GGA', SingleLetterAlphabet()), id='Orangoutan.adn', name='Orangoutan.adn', description='Orangoutan.adn  ORANG-OUTAN', dbxrefs=[]), SeqRecord(seq=Seq('TTCTTTCATGGGGAAGCAAATTTAAGTGCCACCCAAGTATTGGCTCATTCACTA...ACG', SingleLetterAlphabet()), id='Bonobo.adn', name='Bonobo.adn', description='Bonobo.adn  BONOBO', dbxrefs=[]), SeqRecord(seq=Seq('TTCTTTCATGGGGAGACAAATTTGGGTGCCACCCAAGTATTAGCTAACCCACCA...CCC', SingleLetterAlphabet()), id='Gorille.adn', name='Gorille.adn', description='Gorille.adn  GORILLE', dbxrefs=[]), SeqRecord(seq=Seq('TTCTTTCATGGGGAAGCAAATTTAAGTACCACCTAAGTATTGGCCTATTCATTA...CAC', SingleLetterAlphabet()), id='Pan_AJ586556.adn', name='Pan_AJ586556.adn', description='Pan_AJ586556.adn  AJ586556-PAN', dbxrefs=[]), SeqRecord(seq=Seq('TTCTTTCATGGGGAAGCAAATTTAAGTACCACCTAAGTACTGGCTCATTCATTA...CGG', SingleLetterAlphabet()), id='Pan_AJ586557.adn', name='Pan_AJ586557.adn', description='Pan_AJ586557.adn  AJ586557-PAN', dbxrefs=[]), SeqRecord(seq=Seq('TTCTATAGGGGGGAAAGAACTCAAAGAACAACCTAAGTACTAACTTAATCTCCC...AGT', SingleLetterAlphabet()), id='Colobe.adn', name='Colobe.adn', description='Colobe.adn  COLOBE', dbxrefs=[])]
	>>> 
	
	
	>>> list_seq[1]
	SeqRecord(seq=Seq('TTCTTTCATGGGGGATCAGATTTGGGTACCACCCCAGTACCGACCCATTTCCAG...GGA', SingleLetterAlphabet()), id='Orangoutan.adn', name='Orangoutan.adn', description='Orangoutan.adn  ORANG-OUTAN', dbxrefs=[])
	>>> 
	>>> for i in list_seq:
	...     print(i.name)
	... 
	Refhumaine.adn
	Orangoutan.adn
	Bonobo.adn
	Gorille.adn
	Pan_AJ586556.adn
	Pan_AJ586557.adn
	Colobe.adn
	>>> 
	
	
	
	>>> list_seq[1].name
	'Orangoutan.adn'
	>>> 
	
	>>> gp.show(list_seq[1])
	
	
          TTCTTTCATGGGGGATCAGATTTGGGTACCACCCCAGTACCGACCCATTTCCAGCGGCCT
          ----:----|----:----|----:----|----:----|----:----|----:----|
                   10        20        30        40        50        60
	
          ATGTATTTCGTACATTCCTGCCAGCCAACATGAATATCACCCAACACAACAATCGCTTAA
          ----:----|----:----|----:----|----:----|----:----|----:----|
                   70        80        90        100       110       120
	
          CCACCTATAACACATACAAAGCCCAATCCACACCCAACCTCCACCCCCCGCTTACAAGCA
          ----:----|----:----|----:----|----:----|----:----|----:----|
                   130       140       150       160       170       180
	
          AGTACCCCCCCATGCCCCCCCACCCAAACACATACATCGATTCCCCCACATAACCCCTTC
          ----:----|----:----|----:----|----:----|----:----|----:----|
                   190       200       210       220       230       240
	
          CCCCCCCGCATACCAACCAACCCAATCAAGCTTTAAAGTACATAGCACATAACACCCCTA
          ----:----|----:----|----:----|----:----|----:----|----:----|
                   250       260       270       280       290       300
	
          CCGTACATAGCACATTTCTACTAACTCCCTGCTTAACCCTACGGA
          ----:----|----:----|----:----|----:----|----:
                   310       320       330       340
	
	
	 
	 
	>>> c = gp.clustal(*list_seq[1:5])
	>>> 
	>>> gp.showgenix(c)
	
	
    	                    *************   ** ****  ** *****  ****    *  *      *      
	Bonobo.adn              TTCTTTCATGGGGAAGCAAATTTAAGTGCCACCCAAGTATTGGCTCATTCACTA-TAACC
	Pan_AJ586556.adn        TTCTTTCATGGGGAAGCAAATTTAAGTACCACCTAAGTATTGGCCTATTCATTA-CAACC
	Gorille.adn             TTCTTTCATGGGGAGACAAATTTGGGTGCCACCCAAGTATTAGCTAACCCACCAACAATT
	Orangoutan.adn          TTCTTTCATGGGGGATCAGATTTGGGTACCACCCCAGTACCGACCCATTT-CCAGCGGCC
	
    	                       
    	                    ----:----|----:----|----:----|----:----|----:----|----:----|
        	                         10        20        30        40        50        60
	
	
	
    	                       ****** **** **** **** ***** ********      *  ** * *  * * 
	Bonobo.adn              GCTATGTATTTCGTACATTACTGCCAGCCACCATGAATATTATATAGTACTATAATCACT
	Pan_AJ586556.adn        GCTATGTATTTCGTACATTACTGCCAGCCACCATGAATATTGTACAGTACTATAATCACT
	Gorille.adn             GTCATGTATATCGTGCATTACTGCTAGCCACCATGAATAATATACAGTACTATACACACC
	Orangoutan.adn          --TATGTATTTCGTACATTCCTGCCAGCCAACATGAATATCACCCAACACAACAATCGCT
	
        	                ----:----|----:----|----:----|----:----|----:----|----:----|
    	                             70        80        90        100       110       120
	
	
	
    	                     **  **** *** ****  **  *     *****     *       ***   ******
	Bonobo.adn              TAACCACCTATAACACATAAAAACCTACATCCACAT-TAAAACTTT-ACCCCATGCTTAC
	Pan_AJ586556.adn        CAACTACCTATAATACATTAAACCCACCCCCCACAT-TACAACCTCCACCCTATGCTTAC
	Gorille.adn             CAAACACCTGTAATACATACAAAACCCCCTCCACATCTACAGACCC--CCCCATGCTTAC
	Orangoutan.adn          TAACCACCTATAACACATACAAAGCCCAATCCACAC--CCAACCTCCACCCCCCGCTTAC
	
        	                ----:----|----:----|----:----|----:----|----:----|----:----|
            	                     130       140       150       160       170       180
	
	
	
                            ***** * **          ** *  **    * **** * *   *    **      **
    Bonobo.adn              AAGCACGAACAATAATC-AACCTCCAACTGTCAAACATAACACACAAC-TCCAAAGACAC
	Pan_AJ586556.adn        AAGCACGCACAACAGTC-AACCCCCAACTGTCACACATAAAATGCAAC-TCCAAAGACAC
	Gorille.adn             AAGCAAGAACAGTCTTCCAACCCCTAACTAGCACACATTAAATCCAACCTCCCCACTCAC
	Orangoutan.adn          AAGCAAGTACCCCCCCATGCCCCCCCACCCAAACACAT-ACATCGATTCCCCCACATAAC
	
    	                    ----:----|----:----|----:----|----:----|----:----|----:----|
    	                             190       200       210       220       230       240
	
	
	
    	                      ** * ***      *** **** ****        *    * ** ****** ***** 
	Bonobo.adn              TCCTCC-CCCACCCCGATATCAACAAACCTGACAATCCT-TGATAGTACATAGTACATAC
	Pan_AJ586556.adn        CCCTCC-CCCACCCCGATACCAACAAACCTATACCCCCT-TAACAGTACATAGTACATAC
	Gorille.adn             TGCTCCACCCAATGGAATACCAACCAACCTACCTCTTCCACAAAAGCACATAGTACATAA
	Orangoutan.adn          CCCTTC-CCCCCCCGCATACCAACCAACCCAATCAAGCT-TTAAAGTACATAGCACATAA
	
    	                        
    	                    ----:----|----:----|----:----|----:----|----:----|----:----|
        	                         250       260       270       280       290       300
	
	
	
    	                       *    * ** ************ *    ** **  * **   ***      
	Bonobo.adn              AGTCATACACCGTACATAGCACATTACAGTCAAATCCATTCTTGCCCCCACG--
	Pan_AJ586556.adn        AATCGTACATCGCACATAGCACATTACAGTCAAATCCATCCTTGTCCCCAC---
	Gorille.adn             GATCATTCACCGTACATAGCACATTACAGTTAAATC-ATCCTCGTCCC------
	Orangoutan.adn          CACC-CCTACCGTACATAGCACATTTCTACTAACTCCCTGCTTAACCCTACGGA
	
    	                    ----:----|----:----|----:----|----:----|----:----|----
        	                         310       320       330       340       350       
	
	
	
	
	
	
	>>> m = gp.matrix(c)
	
	>>> print(m)
	Bonobo.adn       0
	Pan_AJ586556.adn 0.13276836158192096    0
	Gorille.adn      0.25988700564971756    0.2542372881355932    0
	Orangoutan.adn   0.3220338983050848     0.3220338983050848    0.36723163841807904     0
       	             Bonobo.adn             Pan_AJ586556.adn      Gorille.adn             Orangoutan.adn
	
	
	>>> tree = gp.tree_nj(m)
	
	
	  _______________ Bonobo.adn
	 |
	_|______________ Pan_AJ586556.adn
	 |
	 |         ___________________________________________________ Orangoutan.adn
	 |________|
	          |___________________________________ Gorille.adn
	

## Tutoriel

<a name="creatseq"></a>

### Créer des séquences

Le module genopy dépend de biopython, on peut créer des sequences sous deux types d'objets différents: les objets Seq et SeqRecord. 

Créer un objet Seq

    >>> from genopy import *
       
    >>> s1 = Seq("ATTCGCTCGTAATAGATGGCTCTATATAGATAGGCGCATAGCATCGATGTAGCTTAGCCGT") 
    >>> s1
    Seq('ATTCGCTCGTAATAGATGGCTCTATATAGATAGGCGCATAGCATCGATGTAGCT...CGT')
    
La fonction toSeq de genopy permet de créer plus facilement des objets Seq. Elle ajoute automatiquement le type.
 
    >>> s2 = toSeq("ATTCGCTCGTAATAGATGGCTCTATATAGATAGGCGAATAGCATCGATGTAGCTTACCCGT")
    >>> s2
    Seq('ATTCGCTCGTAATAGATGGCTCTATATAGATAGGCGAATAGCATCGATGTAGCT...CGT', IUPACUnambiguousDNA()) 
    
Les objets SeqRecord contiennent davantage d'informations:

    >>> sr1 = SeqRecord(s1, id="sequence1", name="sequence1", description="une sequence")
    >>> sr1
    SeqRecord(seq=Seq('ATTCGCTCGTAATAGATGGCTCTATATAGATAGGCGCATAGCATCGATGTAGCT...CGT'), id='sequence1', name='sequence1',      description='une sequence', dbxrefs=[])
    
    >>> sr2 = SeqRecord(s2, id="sequence2", name="sequence2", description="une autre sequence")
    >>> sr2
    SeqRecord(seq=Seq('ATTCGCTCGTAATAGATGGCTCTATATAGATAGGCGAATAGCATCGATGTAGCT...CGT', IUPACUnambiguousDNA()), id='sequence2', name='sequence2', description='une sequence', dbxrefs=[])

<a name="creatfasta"></a>

### Sauvegarder un ou plusieurs SeqRecord dans un ou plusieurs fichiers fasta

Créer un fasta contenant un seul SeqRecord. Le fichier fasta serra enregistré dans le repertoire courant et portera le nom du SeqRecord avec l'extension .fas

    >>> mkfas(sr1)

    
Créer plusieurs fichiers fasta mono séquence:
    
    >>> mkfas(sr1, sr2)

Créer un fasta multisequence contenant plusieurs SeqRecord. Le fichier fasta serra enregistré dans le path indiqué en argument:

    >>> mkfasx("myfasta.fas", sr1, sr2)

##### Créer un SeqRecord et un fasta avec la fonction creat()

    >>> help(creat)
    Help on function creat in module genopy:
    
    creat(seq_name, seq, des='', out=False)
        Création d'un SeqRecord et d'un fasta mono séquence
    
    arguments:
            seq_name: (str) nom de la séquence créée
            seq: (str) ou (Seq) ou (SeqRecord) la séquence
            des: (str) déscription du SeqRecord
            out (bool): si False (défaut) , le fasta n'est pas créé
            
    retrun: 
            le SeqRecord créé



    >>> s = creat("un_adn", "ATCTCGTAGCTAGT", des="human adn", out=True)
    >>> s
    SeqRecord(seq=Seq('ATCTCGTAGCTAGT', IUPACUnambiguousDNA()), id='un_adn', name='un_adn', description='human adn', dbxrefs=[])



<a name="convedi"></a>

### Convertir des fichiers anagène avec l'extension .edi en fichiers fasta

On commence par importer les modules genopy et os:

    >>> from genopy import *
    >>> import os

On liste les fichiers du repertoire dans lequel on se trouve (ou un autre repertoire si on indique son path en argument de listdir):

    >>> ld = os.listdir()
    >>> 
    >>> ld
    ['__pycache__', 'genes-Opsines.edi', 'emb.aln', 'comp.aln', 'primates.fas', 'allelesFamillechoree.edi', 'genopy.py', 'tyrfamille4.fas', 'globines-beta-vertebres.fas', 'comp.dnd', 'myfasta.fas']

On remarque ici deux fichiers en .edi
On utilise une boucle for pour selectionner les fichiers .edi et les convertir en .fas grace à la fonction edi2fasta. Cette fonction attend 2 arguments: edi2fasta(file_name.edi, file_name.fas):

    >>> for i in ld:
    ...     if i[len(i)-4:] == ".edi":
    ...             edi2fasta(i, i[:len(i)-4]+".fas")
    ... 
    >>> 

On liste à nouveau les fichiers du repertoire:

    >>> os.listdir()
    ['__pycache__', 'genes-Opsines.edi', 'allelesFamillechoree.fas', 'emb.aln', 'comp.aln', 'primates.fas', 'allelesFamillechoree.edi', 'genopy.py', 'genes-Opsines.fas', 'tyrfamille4.fas', 'globines-beta-vertebres.fas', 'comp.dnd', 'myfasta.fas']
    

<a name="openseq"></a>

### Ouvrir des séquences

Les séquences doivent être hébergées localement.

On commence par rechercher les séquences avec la fonction 'search()' :

    >>> from genopy import *
    >>>
    >>> q = search("Opsine")
    >>> 
    >>> q
    [SeqRecord(seq=Seq('ATGGCCCAGCAGTGGAGCCTCCAAAGGCTCGCAGGCCGCCATCCGCAGGACAGC...TGA', SingleLetterAlphabet()), id='gene_opsine_rouge', name='gene_opsine_rouge', description='gene_opsine_rouge red opsin Homo sapiens opsin 1 (cone pigments), mRNA', dbxrefs=[]), SeqRecord(seq=Seq('ATGGCCCAGCAGTGGAGCCTCCAAAGGCTCGCAGGCCGCCATCCGCAGGACAGC...TGA', SingleLetterAlphabet()), id='gene_opsine_verte', name='gene_opsine_verte', description='gene_opsine_verte red opsin Homo sapiens opsin 1 (cone pigments), mRNA', dbxrefs=[]), SeqRecord(seq=Seq('ATGAGAAAAATGTCGGAGGAAGAGTTTTATCTGTTCAAAAATATCTCTTCAGTG...TGA', SingleLetterAlphabet()), id='gene_opsine_bleue', name='gene_opsine_bleue', description='gene_opsine_bleue Blue opsin Homo sapiens opsin 1 (cone pigments), mRNA', dbxrefs=[])]

On obtient une liste de SeqRecord correspondant à notre recherche.

    >>> len(q)
    3
    >>> 
    >>> for i, j in enumerate(q):
    ...     print(str(i) + " " + j.name)
    ... 
    0 gene_opsine_rouge
    1 gene_opsine_verte
    2 gene_opsine_bleue

On peut alors afficher la séquence du gène de l'opsine verte dont l'indice dans la liste est 1 avec la fonction 'show()':

    >>> show(q[1])
    
    
          ATGGCCCAGCAGTGGAGCCTCCAAAGGCTCGCAGGCCGCCATCCGCAGGACAGCTATGAG
          ----:----|----:----|----:----|----:----|----:----|----:----|
                   10        20        30        40        50        60
    
          GACAGCACCCAGTCCAGCATCTTCACCTACACCAACAGCAACTCCACCAGAGGCCCCTTC
          ----:----|----:----|----:----|----:----|----:----|----:----|
                   70        80        90        100       110       120
    
          GAAGGCCCGAATTACCACATCGCTCCCAGATGGGTGTACCACCTCACCAGTGTCTGGATG
          ----:----|----:----|----:----|----:----|----:----|----:----|
                   130       140       150       160       170       180
    
          ATCTTTGTGGTCATTGCATCCGTTTTCACAAATGGGCTTGTGCTGGCGGCCACCATGAAG
          ----:----|----:----|----:----|----:----|----:----|----:----|
                   190       200       210       220       230       240
    
          TTCAAGAAGCTGCGCCACCCGCTGAACTGGATCCTGGTGAACCTGGCGGTCGCTGACCTG
          ----:----|----:----|----:----|----:----|----:----|----:----|
                   250       260       270       280       290       300
    
          GCAGAGACCGTCATCGCCAGCACTATCAGCGTTGTGAACCAGGTCTATGGCTACTTCGTG
          ----:----|----:----|----:----|----:----|----:----|----:----|
                   310       320       330       340       350       360



(la séquence n'est pas représentée entierement dans ce tutoriel)



### Transcrire, Traduire, Rétro-transcrire

<a name="transc"></a>

    >>> from genopy import *
    >>>
    >>> q = search("Opsine")
    >>> dna = q[0]
    >>> 
    >>> dna
    SeqRecord(seq=Seq('ATGGCCCAGCAGTGGAGCCTCCAAAGGCTCGCAGGCCGCCATCCGCAGGACAGC...TGA', SingleLetterAlphabet()), id='gene_opsine_rouge', name='gene_opsine_rouge', description='gene_opsine_rouge red opsin Homo sapiens opsin 1 (cone pigments), mRNA', dbxrefs=[])

#### Transcrire

    >>> rna = transcribe(dna)
    >>> rna
    SeqRecord(seq=Seq('AUGGCCCAGCAGUGGAGCCUCCAAAGGCUCGCAGGCCGCCAUCCGCAGGACAGC...UGA', RNAAlphabet()), id='gene_opsine_rouge.rna', name='gene_opsine_rouge.rna', description='transcription de [gene_opsine_rouge red opsin Homo sapiens opsin 1 (cone pigments), mRNA]', dbxrefs=[])


<a name="transl"></a>

#### Traduire

    >>> help(translate)
    Help on function translate in module genopy:
    
    translate(*seq, table_id=1, to_stop=True, stop_symbol=' ', cds=True, out=False)
    Traduction d'un ou plusieurs dna ou rna SeqRecord. 
    
    arguments:
        *seq: liste de dna ou rna Seq ou SeqRecord
        table_id: indice de la table du code génétique à utiliser (par defaut: 1 c'est le code standard)
        to_stop: bool  True par défaut: la traduction s'arrete au codon stop
        stop_symbol: représentation de l'arret de la traduction
        cds: coding sequence: bool  True par défaut: la traduction commence au condon initiateur et se termine au codon stop
        out: bool False par défaut. Si out == True, chaque sequence peptidique est sauvegardée dans un fasta. 
        
    return: proteine SeqRecord ou une liste de proteines SeqRecord
 
> ORF: open reading frame ou phase ouverte de lecture, c'est la séquence du brin non transcrit (brin codant) de l'adn (ou de l'arn pré-messager) située entre le premier codon initiateur jusqu'au premier codon stop 
> 
> CDS: coding sequence ou séquence codante, c'est l'ORF sans les introns

Il est possible de choisir la table du code génétique utilisée (par défaut c'est le code standard)

L'argument to_stop = True permet d'arreter la traduction au codon stop

L'argument stop_symbol permet de spécifier la représentation de l'arret de la traduction

L'argument cds = True implique que la traduction commencera au premier codon initiateur rencontré et se terminera au premier codon stop rencontré sur le cadre de lecture. La majorité des séquences importées depuis anangene sont déjà des CDS, donc ca change rien; mais ce n'est pas le cas des fasta téléchargés depuis les banques de séquences. 

    >>> pep = translate(rna)
    >>> pep
    SeqRecord(seq=Seq('MAQQWSLQRLAGRHPQDSYEDSTQSSIFTYTNSNSTRGPFEGPNYHIAPRWVYH...SPA', ExtendedIUPACProtein()), id='gene_opsine_rouge.rna.prot', name='gene_opsine_rouge.rna.prot', description='traduction de [transcription de [gene_opsine_rouge red opsin Homo sapiens opsin 1 (cone pigments), mRNA]]', dbxrefs=[])

    >>> show(pep)
    
    
          MetAlaGlnGlnTrpSerLeuGlnArgLeuAlaGlyArgHisProGlnAspSerTyrGluAspSerThrGlnSerSerIlePheThrTyr
           -  -  -  -  :  -  -  -  -  |  -  -  -  -  :  -  -  -  -  |  -  -  -  -  :  -  -  -  -  | 
                                      10                            20                            30
    
          ThrAsnSerAsnSerThrArgGlyProPheGluGlyProAsnTyrHisIleAlaProArgTrpValTyrHisLeuThrSerValTrpMet
           -  -  -  -  :  -  -  -  -  |  -  -  -  -  :  -  -  -  -  |  -  -  -  -  :  -  -  -  -  | 
                                      40                            50                            60
    
          IlePheValValThrAlaSerValPheThrAsnGlyLeuValLeuAlaAlaThrMetLysPheLysLysLeuArgHisProLeuAsnTrp
           -  -  -  -  :  -  -  -  -  |  -  -  -  -  :  -  -  -  -  |  -  -  -  -  :  -  -  -  -  | 
                                      70                            80                            90
    
          IleLeuValAsnLeuAlaValAlaAspLeuAlaGluThrValIleAlaSerThrIleSerIleValAsnGlnValSerGlyTyrPheVal
           -  -  -  -  :  -  -  -  -  |  -  -  -  -  :  -  -  -  -  |  -  -  -  -  :  -  -  -  -  | 
                                      100                           110                           120
    
          LeuGlyHisProMetCysValLeuGluGlyTyrThrValSerLeuCysGlyIleThrGlyLeuTrpSerLeuAlaIleIleSerTrpGlu
           -  -  -  -  :  -  -  -  -  |  -  -  -  -  :  -  -  -  -  |  -  -  -  -  :  -  -  -  -  | 
                                      130                           140                           150

<a name="rt"></a>

#### Rétro transcrire

    >>> rna
    SeqRecord(seq=Seq('AUGGCCCAGCAGUGGAGCCUCCAAAGGCUCGCAGGCCGCCAUCCGCAGGACAGC...UGA', RNAAlphabet()), id='gene_opsine_rouge.rna', name='gene_opsine_rouge.rna', description='transcription de [gene_opsine_rouge red opsin Homo sapiens opsin 1 (cone pigments), mRNA]', dbxrefs=[])
    
    
    >>> dna_c = retro_transcribe(rna)
    >>> dna_c
    SeqRecord(seq=Seq('ATGGCCCAGCAGTGGAGCCTCCAAAGGCTCGCAGGCCGCCATCCGCAGGACAGC...TGA', DNAAlphabet()), id='gene_opsine_rouge.rna.dnac', name='gene_opsine_rouge.rna.dnac', description='retro transcription de [transcription de [gene_opsine_rouge red opsin Homo sapiens opsin 1 (cone pigments), mRNA]]', dbxrefs=[])
    
    >>> dna.seq == dna_c.seq
    True

<a name="tablecode"></a>

### Afficher les tables du code génétique

    >>> from genopy import *
    >>>
    >>> codegen()
    
    
    1 : Standard
    2 : Vertebrate Mitochondrial
    3 : Yeast Mitochondrial
    4 : Mold Mitochondrial
    5 : Invertebrate Mitochondrial
    6 : Ciliate Nuclear
    9 : Echinoderm Mitochondrial
    10 : Euplotid Nuclear
    11 : Bacterial
    12 : Alternative Yeast Nuclear
    13 : Ascidian Mitochondrial
    14 : Alternative Flatworm Mitochondrial
    15 : Blepharisma Macronuclear
    16 : Chlorophycean Mitochondrial
    21 : Trematode Mitochondrial
    22 : Scenedesmus obliquus Mitochondrial
    23 : Thraustochytrium Mitochondrial
    24 : Pterobranchia Mitochondrial
    25 : Candidate Division SR1
    26 : Pachysolen tannophilus Nuclear
    27 : Karyorelict Nuclear
    28 : Condylostoma Nuclear
    29 : Mesodinium Nuclear
    30 : Peritrich Nuclear
    31 : Blastocrithidia Nuclear




    Entrez le numéro de la table à afficher: 1
    
    
    RNA table: 
    
    Table 1 Standard, SGC0
    
      |  U      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    U | UUU F   | UCU S   | UAU Y   | UGU C   | U
    U | UUC F   | UCC S   | UAC Y   | UGC C   | C
    U | UUA L   | UCA S   | UAA Stop| UGA Stop| A
    U | UUG L(s)| UCG S   | UAG Stop| UGG W   | G
    --+---------+---------+---------+---------+--
    C | CUU L   | CCU P   | CAU H   | CGU R   | U
    C | CUC L   | CCC P   | CAC H   | CGC R   | C
    C | CUA L   | CCA P   | CAA Q   | CGA R   | A
    C | CUG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | AUU I   | ACU T   | AAU N   | AGU S   | U
    A | AUC I   | ACC T   | AAC N   | AGC S   | C
    A | AUA I   | ACA T   | AAA K   | AGA R   | A
    A | AUG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GUU V   | GCU A   | GAU D   | GGU G   | U
    G | GUC V   | GCC A   | GAC D   | GGC G   | C
    G | GUA V   | GCA A   | GAA E   | GGA G   | A
    G | GUG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--
    
    
    <Bio.Data.CodonTable.NCBICodonTableRNA object at 0x74a047d0>
    
    
    
    
    >>> cg = codegen(1, bydna=True)
    
    
    DNA table: 
    
    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--


    >>> cg.
    cg.back_table           cg.id                   cg.nucleotide_alphabet  cg.start_codons         
    cg.forward_table        cg.names                cg.protein_alphabet     cg.stop_codons          

    >>> cg.start_codons
    ['TTG', 'CTG', 'ATG']

    >>> cg.back_table
    {'K': 'AAG', 'N': 'AAT', 'T': 'ACT', 'R': 'CGT', 'S': 'TCT', 'I': 'ATT', 'M': 'ATG', 'Q': 'CAG', 'H': 'CAT', 'P': 'CCT', 'L': 'TTG', 'E': 'GAG', 'D': 'GAT', 'A': 'GCT', 'G': 'GGT', 'V': 'GTT', 'Y': 'TAT', 'C': 'TGT', 'W': 'TGG', 'F': 'TTT', None: 'TAA'}

    >>> cg.forward_table
    {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

    >>> cg.forward_table["AAC"]
    'N' 




<a name="needlewater"></a>

### Comparer deux séquences: alignement par Needle ou Water

On commence par rechercher des séquences:

    >>> from genopy import *
    >>>
    >>> q = search("primates")
    >>> 
    >>> len(q)
    7
    >>> 
    >>> q[0]
    SeqRecord(seq=Seq('TTCTTTCATGGGGAAGCAGATTTGGGTGCCACCCAAGTATTGACTCACCCATCA...CAT', SingleLetterAlphabet()), id='Refhumaine.adn', name='Refhumaine.adn', description=' Refhumaine.adn', dbxrefs=[])

On va ensuite comparer les séquences 2 à 2. 

Il est possible de faire un alignement global (needle) ou local (water)

    >>> help(needle)
    Help on function needle in module genopy:
    
    needle(*seq, gapopen=10, gapextend=0.5, out='emb.aln')
        Alignement global par la méthode de Needleman
    
    arguments:
            
            seq: couple de 2 SeqRecord à aligner
            gapopen: pénalité de gap
            gapextend: pénalité d'expansion
            out: nom du fichier emboss créé
            
    return: un objet align

L'alignement est enregistré par défaut dans le fichier emb.aln (qui serra écrasé au prochain alignement réalisé)

Dans les 2 cas on peut eventuellement modifier les pénalités pour les ouvertures ou les extensions de gap.


    >>> n = needle(q[0], q[1])
    >>> 
    >>> w = water(q[0], q[1])
    >>> 
    >>> n
    <<class 'Bio.Align.MultipleSeqAlignment'> instance (2 records of length 365, SingleLetterAlphabet()) at 70615df0>
    >>> w
    <<class 'Bio.Align.MultipleSeqAlignment'> instance (2 records of length 361, SingleLetterAlphabet()) at 7061def0>


On obtient des objets "align" que l'on peut alors exploiter de différentes façons:

Obtenir des informations de base sur l'alignement:

    >>> print(n)
    SingleLetterAlphabet() alignment with 2 rows and 365 columns
    TTCTTTCATGGGGAAGCAGATTTGGGTGCCACCCAAGTATTGAC...--- Refhumaine.adn
    TTCTTTCATGGGGGATCAGATTTGGGTACCACCCCAGT----AC...GGA Orangoutan.adn
    
    >>> n.
    n.add_sequence(          n.append(                n.extend(                n.get_alignment_length(  
    n.annotations            n.column_annotations     n.format(                n.sort(                  
    
    >>> n.get_alignment_length()
    365
    
    >>> len(n)
    2
    
    >>> n[0]
    SeqRecord(seq=Seq('TTCTTTCATGGGGAAGCAGATTTGGGTGCCACCCAAGTATTGACTCACCCATCA...---', SingleLetterAlphabet()), id='Refhumaine.adn', name='<unknown name>', description='Refhumaine.adn', dbxrefs=[])
 


Afficher l'alignement avec les fonctions shogenix(), showanag(), (ou showgenix3() pour les séquences peptidiques avec acides aminés représentés à trois lettres)

L'affichage via showgenix représente toutes les séquences. Les similitudes sont symbolisées par une * et les délétions par des - . 

    >>> showgenix(n)


                          ************* * *********** ****** ***    **  ******   ** * 
    Refhumaine.adn        TTCTTTCATGGGGAAGCAGATTTGGGTGCCACCCAAGTATTGACTCACCCATCAACAACC
    Orangoutan.adn        TTCTTTCATGGGGGATCAGATTTGGGTACCACCCCAGT----ACCGACCCATTTCCAGCG

                          ----:----|----:----|----:----|----:----|----:----|----:----|
                                   10        20        30        40        50        60



                          * ****************** ********** *********   **   ** **  *** 
    Refhumaine.adn        G-CTATGTATTTCGTACATTACTGCCAGCCACCATGAATATTGTACGGTAC-CATAAAT-
    Orangoutan.adn        GCCTATGTATTTCGTACATTCCTGCCAGCCAACATGAATAT--CACCCAACACAACAATC

                          ----:----|----:----|----:----|----:----|----:----|----:----|
                                   70        80        90        100       110       120


L'affichage via showanag utilise le même mode de symbolisation que le logiciel anagene.

    >>> showanag(n)


                          ************* * *********** ****** ***    **  ******   ** * 
    Refhumaine.adn        TTCTTTCATGGGGAAGCAGATTTGGGTGCCACCCAAGTATTGACTCACCCATCAACAACC
    Orangoutan.adn        -------------G-T-----------A------C---____--CG------TTC--G-G

                          ----:----|----:----|----:----|----:----|----:----|----:----|
                                   10        20        30        40        50        60



                          * ****************** ********** *********   **   ** **  *** 
    Refhumaine.adn        G_CTATGTATTTCGTACATTACTGCCAGCCACCATGAATATTGTACGGTAC_CATAAAT_
    Orangoutan.adn        -C------------------C----------A---------__C--CCA--A--AC---C

                          ----:----|----:----|----:----|----:----|----:----|----:----|
                                   70        80        90        100       110       120
                                                                  




Il est aussi possible de lire l'alignement enregistré dans le fichier emb.aln avec la fonction reademboss() 


   
    >>> n = needle(q[0], q[1])
    
    >>> reademboss("emb.aln")    # ou reademboss()
    
    ########################################
    # Program: needle
    # Rundate: Sun 14 Jul 2019 11:13:59
    # Commandline: needle
    #    -outfile emb.aln
    #    -asequence seqa.fas
    #    -bsequence seqb.fas
    #    -gapopen 10
    #    -gapextend 0.5
    # Align_format: srspair
    # Report_file: emb.aln
    ########################################

    #=======================================
    #
    # Aligned_sequences: 2
    # 1: Refhumaine.adn
    # 2: Orangoutan.adn
    # Matrix: EDNAFULL
    # Gap_penalty: 10.0
    # Extend_penalty: 0.5
    #
    # Length: 365
    # Identity:     255/365 (69.9%)
    # Similarity:   255/365 (69.9%)
    # Gaps:          40/365 (11.0%)
    # Score: 824.5
    # 
    #
    #=======================================

    Refhumaine.ad      1 TTCTTTCATGGGGAAGCAGATTTGGGTGCCACCCAAGTATTGACTCACCC     50
                         |||||||||||||.|.|||||||||||.||||||.|||    ||..||||
    Orangoutan.ad      1 TTCTTTCATGGGGGATCAGATTTGGGTACCACCCCAGT----ACCGACCC     46

    Refhumaine.ad     51 ATCAACAACCG-CTATGTATTTCGTACATTACTGCCAGCCACCATGAATA     99
                         ||...||.|.| ||||||||||||||||||.||||||||||.||||||||
    Orangoutan.ad     47 ATTTCCAGCGGCCTATGTATTTCGTACATTCCTGCCAGCCAACATGAATA     96


    #---------------------------------------
    #---------------------------------------

    <<class 'Bio.Align.MultipleSeqAlignment'> instance (2 records of length 365, SingleLetterAlphabet()) at 70603ab0>


(Les alignements ne sont pas représentés en entier dans ce tutoriel)

Il est possible d'afficher le % de similitudes ou de différences en utilisant la fonction reademboss (comme ci dessous) ou en affichant la matrice de distance avec la fonction matrix()

    >>> m = matrix(n)
    >>> print(m)
    Refhumaine.adn  0
    rangoutan.adn   0.3013698630136986      0
                    Refhumaine.adn          Orangoutan.adn

<a name="clustal"></a>

### Comparaison multiple de séquences avec Clustal

On commence par importer des séquences:

    >>> from genopy import *
    >>>
    >>> q = search("tyr")
    >>> 
    >>> len(q)
    25
    >>> 
    >>> q[0]
    SeqRecord(seq=Seq('ATGCTCCTGGCTGTTTTGTACTGCCTGCTGTGGAGTTTCCAGACCTCCGCTGGC...TAA', SingleLetterAlphabet()), id='Tyrcod2', name='Tyrcod2', description="Tyrcod2  Partie strictement codante d'un allele du gene de la tyrosinase (ref erence 2).", dbxrefs=[])

Puis on utilise la fonction clustal: exemple ici pour comparer les séquence d’index 5, 6, 7, 8

    >>> c = clustal(*q[5:9])  # ou clustal(q[5], q[6], q[7], q[8])
    >>> c
    <<class 'Bio.Align.MultipleSeqAlignment'> instance (4 records of length 1590, SingleLetterAlphabet()) at 704f06f0>
    

Affichage simple: 

    >>> print(c)
    SingleLetterAlphabet() alignment with 4 rows and 1590 columns
    ATGCTCCTGGCTGTTTTGTACTGCCTGCTGTGGAGTTTCCAGAC...TAA F4_I2all1
    ATGCTCCTGGCTGTTTTGTACTGCCTGCTGTGGAGTTTCCAGAC...TAA F4_I2all2
    ATGCTCCTGGCTGTTTTGTACTGCCTGCTGTGGAGTTTCCAGAC...TAA F4_I1all2
    ATGCTCCTGGCTGTTTTGTACTGCCTGCTGTGGAGTTTCCAGAC...TAA F4_I1all1
    >>> 

Méthodes applicables à l'objet align:

    >>> c.
    c.add_sequence(          c.append(                c.extend(                c.get_alignment_length(  
    c.annotations            c.column_annotations     c.format(                c.sort(                  

Affichage au format clustal:

    >>> print(c.format("clustal"))
    CLUSTAL 2.1 multiple sequence alignment

   
    F4_I2all1     ATGCTCCTGGCTGTTTTGTACTGCCTGCTGTGGAGTTTCCAGACCTCCGC
    F4_I2all2     ATGCTCCTGGCTGTTTTGTACTGCCTGCTGTGGAGTTTCCAGACCTCCGC
    F4_I1all2     ATGCTCCTGGCTGTTTTGTACTGCCTGCTGTGGAGTTTCCAGACCTCCGC
    F4_I1all1     ATGCTCCTGGCTGTTTTGTACTGCCTGCTGTGGAGTTTCCAGACCTCCGC
                  **************************************************


Affichage au format phylip:

    >>> print(c.format("phylip"))
    4 1590
    F4_I2all1  ATGCTCCTGG CTGTTTTGTA CTGCCTGCTG TGGAGTTTCC AGACCTCCGC
    F4_I2all2  ATGCTCCTGG CTGTTTTGTA CTGCCTGCTG TGGAGTTTCC AGACCTCCGC
    F4_I1all2  ATGCTCCTGG CTGTTTTGTA CTGCCTGCTG TGGAGTTTCC AGACCTCCGC
    F4_I1all1  ATGCTCCTGG CTGTTTTGTA CTGCCTGCTG TGGAGTTTCC AGACCTCCGC

               TGGCCATTTC CCTAGAGCCT GTGTCTCCTC TAAGAACCTG ATGGAGAAGG
               TGGCCATTTC CCTAGAGCCT GTGTCTCCTC TAAGAACCTG ATGGAGAAGG
               TGGCCATTTC CCTAGAGCCT GTGTCTCCTC TAAGAACCTG ATGGAGAAGG
               TGGCCATTTC CCTAGAGCCT GTGTCTCCTC TAAGAACCTG ATGGAGAAGG


Affichage au format genix:

    >>> showgenix(c)


                     ************************************************************
    F4_I2all1        ATGCTCCTGGCTGTTTTGTACTGCCTGCTGTGGAGTTTCCAGACCTCCGCTGGCCATTTC
    F4_I2all2        ATGCTCCTGGCTGTTTTGTACTGCCTGCTGTGGAGTTTCCAGACCTCCGCTGGCCATTTC
    F4_I1all2        ATGCTCCTGGCTGTTTTGTACTGCCTGCTGTGGAGTTTCCAGACCTCCGCTGGCCATTTC
    F4_I1all1        ATGCTCCTGGCTGTTTTGTACTGCCTGCTGTGGAGTTTCCAGACCTCCGCTGGCCATTTC

                     ----:----|----:----|----:----|----:----|----:----|----:----|
                              10        20        30        40        50        60

Il est aussi possible de faire un affichage de l'alignement au format anagene avec la fonction showanag().

<a name="matrix"></a>

#### Afficher une matrice de distances:

La fonction matrix attend un alignement (needle, water ou clustal) en argument obligatoire:

    >>> m = matrix(c)
    >>> print(m)
    F4_I2all1       0
    F4_I2all2       0.0006289308176100628   0
    F4_I1all2       0.0012578616352201255   0.0006289308176100628   0
    F4_I1all1       0.0012578616352201255   0.0006289308176100628   0.0012578616352201255   0
                    F4_I2all1               F4_I2all2               F4_I1all2               F4_I1all1

La matrice peut être convertie en table HTML et déployée sur [IPFS](https://gist.github.com/YannBouyeron/53e6d67782dcff5995754b0a7333fa0b) avec la fonction matrix2ipfs:

    >>> matrix2ipfs(m)
    'QmcMhGpqsiCoyoPmbUf4xuX7ayMNiVHEkU44t4XMyWm2AK'

On obtient le hash du fichier HTML, consultable sur votre noeud ipfs ou via ipfs.io: https://ipfs.io/ipfs/<hash>:

[https://ipfs.io/ipfs/QmcMhGpqsiCoyoPmbUf4xuX7ayMNiVHEkU44t4XMyWm2AK](https://ipfs.io/ipfs/QmcMhGpqsiCoyoPmbUf4xuX7ayMNiVHEkU44t4XMyWm2AK)

La matrice peut aussi être sauvegardée sous forme d'image png qui serra enregistrée localement:

    >>> matrix2png(m, path="mamatrice.png")
    <matplotlib.axes._subplots.AxesSubplot object at 0x704f0a50>


<p align="center">
  <img src="Images/62FFEC5E-3F75-4173-9CF2-222939C40CF8.png
">
</p>

<a name="phylomatrix"></a>

### Construire un arbre phylogénétique à partir d'une matrice de distances

    >>> from genopy import *
    >>>
    >>> q = search("tyr")
    >>> c = clustal(*q[5:10])
    >>> m = matrix(c)
    
    >>> print(m)
    F4_I2all2       0
    4_II1all1       0.0                     0
    F4_I2all1       0.0006289308176100628   0.0006289308176100628   0
    F4_I1all2       0.0006289308176100628   0.0006289308176100628   0.0012578616352201255   0
    F4_I1all1       0.0006289308176100628   0.0006289308176100628   0.0012578616352201255   0.0012578616352201255   0
                    F4_I2all2               F4_II1all1              F4_I2all1               F4_I1all2               F4_I1all1


Neighbor-joining:

    >>> help(tree_nj)
    
    Help on function tree_nj in module genopy:

    tree_nj(dm)
        Création d'un arbre nj à partir de la matrice de distance retournée par la fonction matrix
        return: tree

    
    
    >>> nj = tree_nj(m)


     , F4_I2all2
     |
     | F4_II1all1
     |
     |_________________________________________________________________ F4_I2all1
    _|
     |_________________________________________________________________ F4_I1all2
     |
     |_________________________________________________________________ F4_I1all1


    >>> nj
    Tree(rooted=True)


UPGMA:

    >>> help(tree_upgma)
    
    Help on function tree_upgma in module genopy:

    tree_upgma(dm)
        Création d'un arbre upgma à partir de la matrice de distance retournée par la fonction matrix
        return: tree


    >>> up = tree_upgma(m)


      _________________________________________________________________ F4_I2all1
    _|
     |         ________________________________________________________ F4_I1all2
     |________|
              |                   _____________________________________ F4_I1all1
              |__________________|
                                 |                                     , F4_II1all1
                                 |_____________________________________|
                                                                       | F4_I2all2


Parcimonie:

    >>> help(parsimony_tree)
    
    Help on function parsimony_tree in module genopy:

    parsimony_tree(align, starting_tree)
        Création d'un arbre basé sur la méthode de parcimonie à partir de l'alignement et d'un arbre de départ (nj)
    
    arguments:
            align: objet align retourné par clustal()
            satrting_tree: objet tree retourné par tree_nj()
    
    return: tree


    >>> pa = parsimony_tree(c, nj)


      _________________________________________________________________ F4_I2all1
     |
     , F4_II1all1
     |
     | F4_I2all2
    _|
     |_________________________________________________________________ F4_I1all1
     |
     |_________________________________________________________________ F4_I1all2



Pour afficher les distances sur les branches de l'arbre:

     >>> nj = tree_nj(m)
     
     >>> dnj = draw_tree(nj, distance=True)


<a name="phyloseq"></a>

### Construire un arbre phylogénétique à partir d'une liste de séquences

    >>> from genopy import *
    >>>
    >>> q = search("tyr")
    >>> 
    >>> len(q)
    25
    >>> 
    >>> help(phylo)
    
    Help on function phylo in module genopy:

    phylo(*seq, out='comp.aln')
        Méthode simplifiée de création d'un arbre phylogénétique à partir d'une liste de SeqRecord
    
            argument: *seq: liste de SeqRecord à aligner
            return: (align, stdout, tree)


    >>> p = phylo(*q[3:7])


               Tyralb_TS  F4_I1all2  Tyralb_A5  F4_I1all1
    Tyralb_TS   0.000000   0.000000   0.001887   0.001258
    F4_I1all2   0.000000   0.000000   0.001887   0.001258
    Tyralb_A5   0.001887   0.001887   0.000000   0.000629
    F4_I1all1   0.001258   0.001258   0.000629   0.000000




                                                   ______________________ Tyralb_A5
      ____________________________________________|
     |                                            | F4_I1all1
    _|
     | Tyralb_TS
     |
     | F4_I1all2

