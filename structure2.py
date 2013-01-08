'''Project steps:
#First of all we provide a human ensembl gene ID and some species we are interested in. After that the script will gather the ortholog gene ID from ensembl database for those species.

#Then for each gene it will download all the proteins coded by those genes. This information will be stored in a hierarchichal structure.

#From the protein length we have to select the most suitable combinations and align them. We will use MAFFT (or PRANK) for that purpose.

#We need a hierarchical list/dict class structure
#One specie have several genes, one gene several proteins and one protein corresponds to one transcript.
'''
import csv

class Protein:
    #The class protein contains information about the length of the protein and a list with it's exons (with it's own attributes)
    def __init__(self, name, len):
        self.name = name
        self.len = len

class Gene:
    #The class gene contains information about the gene and a dict with it's proteins (with it's own attributes)
    def __init__(self, name):
        self.name = name
        self.proteins = {}

    @property
    def nproteins(self):
        return len(self.proteins)	
		
class Species:
    '''This structure contains all the information needed for all genes.
    One specie have several genes, one gene several proteins'''

    def __init__(self, name):
        self.name = name #name of the GENE
        self.genes = {}

    def addProtein(self, gene, protname, len):
        #Converting a line from the input file into a protein and/or an exon
        if gene not in self.genes:
            self.genes[gene] = Gene(gene) 
        self.genes[gene].proteins[protname] = Protein(protname, len)

    @property
    def ngenes(self):
        return len(self.genes)

def import_gene(specie, filename):
	#Reads input file and stores information in one species structure
	specie = Species(specie)
	for line in csv.reader(open(filename), delimiter='\t'):
		try:
			specie.addProtein(line[0], line[1], int(line[2]))
		except:
			pass
	return specie

def homo_load(gene_list, folder):
	homo_load = {}
	#Load all the homologs.txt information into a dictionary structure called homologs
	for line in csv.reader(open(folder+'/homologs.txt'), delimiter='\t'):
		#exec 'line_sp'
		homo_load[line[0]]=line
		gene_list.append(line[0])
	return homo_load, gene_list
	return sp_copy