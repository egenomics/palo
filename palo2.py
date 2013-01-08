from structure2 import *
from heuristics3 import *
import sys, os

folder ='files/' #Define here where your input files are. Output will also be written here

#species list
homologs = {}
gene_list = []

#Loading the species information
print 'Loading species information. Be patient, it may take a while'
species = import_gene('species', folder+'/species.txt')

#loading homologs info from file and creating the gene_list
(homologs, gene_list) = homo_load(gene_list, folder)
no_processed=0

#Creating output file
output = open(folder+'palo_combs.txt', "w")

#Count total number of jobs
total = float(len(gene_list))
progress=0.0
numberofgenes = 0
file2 = open(folder+'/unprocessed_genes.txt', "a")
file2.close()


print 'Starting calculations'
for geneId in gene_list:
    numberofgenes = numberofgenes + 1
    proteins_palo = []
    if (numberofgenes % 50) == 0:
        progress = int(((numberofgenes/total)*100.0))
        print str(numberofgenes)+' genes processed. Progress is ' + str(progress) +'%'
    combinations=[]
    num_combinations = 1
    try:
        for gene_sp in homologs[geneId]: #Calculate number of possible combinations of a gene
            num_combinations = num_combinations * species.genes[gene_sp].nproteins
    except:
        num_combinations=-1
        file2 = open(folder+'/unprocessed_genes.txt', "a")
        print>>file2, geneId, num_combinations, 'transcript error'
        file2.close()
        no_processed = no_processed +1
    try:
        if (0<num_combinations<=100000000): #Limit genes with a very high number of combinations
            (proteins_palo) = calc_prots(geneId, species, homologs)
            print>>output, geneId+'\t'+('\t'.join(proteins_palo))
        elif num_combinations==-1:
            pass
        else:
            file2 = open(folder+'/unprocessed_genes.txt', "a")
            print>>file2, geneId, num_combinations, 'skipped'
            file2.close()
            no_processed = no_processed +1
    except:
        file2 = open(folder+'/unprocessed_genes.txt', "a")
        print>>file2, geneId, num_combinations, 'other error'
        file2.close()
        no_processed = no_processed +1

f = open(folder+'/stats.txt', 'w')
print>>f, 'Processed genes: '+ str(numberofgenes-no_processed)
print>>f, 'Unprocessed genes: '+ str(no_processed)
print>>f, 'Total genes requested: '+ str(numberofgenes)
f.close()

print 'Job done, check your directory for output, stats and errors.'