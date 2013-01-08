import numpy as np
def calc_prots(gene, species, homologs):
    #This function calculates the more similar in lengths protein combination among different species (palo)
    proteins_palo = 0
    candidates_len = []
    candidates_name = []
    for index, homolog in enumerate(homologs[gene]):
        #Initialize proteins_default with the same number of indexes as number of species
        candidates_len.append([])
        candidates_name.append([])
        for prot in species.genes[homologs[gene][index]].proteins:
            candidates_len[index].append(species.genes[homologs[gene][index]].proteins[prot].len)
            candidates_name[index].append(species.genes[homologs[gene][index]].proteins[prot].name)
    #Creating the possible combinations and store the lists of lengths in vector r
    items=cartesian(candidates_len)
    names=cartesian(candidates_name)
    best = np.inf
    proteins_palo = None
    for i,item in enumerate(items):
        ls = (np.subtract.outer(item,item) ** 2).sum()
        if ls < best:
            best = ls
            proteins_palo = i
    return names[proteins_palo]

def cartesian(arrays, out=None):
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype
    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)
    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out