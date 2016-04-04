# SJS
# Various functions, classes used in simulations

import re
import sys
import numpy as np
from scipy import linalg
from copy import deepcopy
from random import uniform, shuffle
amino_acids  = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
codon_dict = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
genetic_code = [["GCA", "GCC", "GCG", "GCT"], ["TGC","TGT"], ["GAC", "GAT"], ["GAA", "GAG"], ["TTC", "TTT"], ["GGA", "GGC", "GGG", "GGT"], ["CAC", "CAT"], ["ATA", "ATC", "ATT"], ["AAA", "AAG"], ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"], ["ATG"], ["AAC", "AAT"], ["CCA", "CCC", "CCG", "CCT"], ["CAA", "CAG"], ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"] , ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"], ["ACA", "ACC", "ACG", "ACT"], ["GTA", "GTC", "GTG", "GTT"], ["TGG"], ["TAC", "TAT"]]
ZERO = 1e-10



class dNdS_from_MutSel():
    '''
        Class to compute dN, dS, and/or dN/dS from mutation-selection model parameters, as described in Spielman and Wilke 2015.

        Positional arguments:
            codon_frequencies = either a list, numpy array, or dictionary of equilibrium codon frequencies. If list or numpy array, frequencies should be in alphabetical order, excluding stop codons. Relatively little sanity checking is done here, so provide something reasonable..
            mutation_dictionary = (optional) dictionary of mutation rates. This argument is analogous to the "mu" dictionary provided to pyvolve when simulating with custom mutation rates. If this argument is not provided, assumes equal mutation rates.
    '''
    def __init__(self, codon_frequencies, mutation_dictionary = None):
        
        self.ZERO = 1e-12
        self.codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
        self.codon_dict = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}

        # Frequency setup
        if type(codon_frequencies) is dict:
            self.codon_freqs_dict = codon_frequencies
        else:
            self.codon_freqs_dict = {}
            for c in range(len(self.codons)):
                self.codon_freqs_dict[self.codons[c]] = codon_frequencies[c]
        assert(abs(1. - np.sum(self.codon_freqs_dict.values())) < self.ZERO), "\n\nProvided codon frequencies must sum to 1."
        
        
        # Mutation setup        
        if mutation_dictionary is None:
            self.mu_dict = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}
        else:
            self.mu_dict = deepcopy(mutation_dictionary)
            for key in mutation_dictionary: 
                if key[1] + key[0] not in self.mu_dict:
                    self.mu_dict[key[1] + key[0]] = mutation_dictionary[key]



    def compute_dn(self):
        ''' 
            Compute dN from mutation-selection parameters: a dictionary of equilibrium codon frequencies and a dictionary of mutation rates.
        '''
        return self.compute_quantity("nonsyn")




    def compute_ds(self):
        ''' 
            Compute dS from mutation-selection parameters: a dictionary of equilibrium codon frequencies and a dictionary of mutation rates.
        '''
        return self.compute_quantity("syn")




    def compute_dnds(self):
        ''' 
            Compute dN/dS from mutation-selection parameters: a dictionary of equilibrium codon frequencies and a dictionary of mutation rates.
        '''
        dn = self.compute_dn()
        ds = self.compute_ds()
        
        if dn <= self.ZERO:
            return 0.
        elif ds <= self.ZERO:
            return np.inf
        else:
            return dn/ds


    def compute_quantity(self, type):
        '''
            Compute either dN or dS (type is 'nonsyn' or 'syn', respectively).
        '''
        numer = 0.
        denom = 0.

        for codon in self.codon_freqs_dict:
            rate = 0.
            sites = 0.
            rate, sites = self.calc_paths(codon, type)
            numer += rate
            denom += sites
        
        assert( denom > self.ZERO ), "\n\nProvided frequencies indicate no evolution is 'possible'."
        return numer/denom




     
    def calc_paths(self, source, type):
        ''' 
            Compute a term in the quantity numerator for a given source codon.
        '''
        rate = 0.
        sites = 0.
        source_freq = self.codon_freqs_dict[source]
        for target in self.codons:
            diff = self.get_nuc_diff(source, target) # only consider single nucleotide differences since are calculating instantaneous.
            if (type == 'nonsyn' and self.codon_dict[source] != self.codon_dict[target]) or (type == 'syn' and self.codon_dict[source] == self.codon_dict[target]):
                if len(diff) == 2:
                    rate  += self.calc_subst_prob( source_freq, self.codon_freqs_dict[target], self.mu_dict[diff], self.mu_dict[diff[1]+diff[0]] )
                    sites += self.mu_dict[diff]
        rate  *= source_freq
        sites *= source_freq
        return rate, sites



    def get_nuc_diff(self, source, target):
        '''
            Determine nucleotide difference between two codons.
        '''
        diff = ''
        for i in range(3):
            if source[i] != target[i]: 
                diff += source[i]+target[i]
        return diff

 
    
    def calc_subst_prob(self, pi, pj, mu_ij, mu_ji):
        '''
            Compute the substitution probability between two codons, as in Halpern and Bruno 1998.
        '''
        if abs(pi) <= self.ZERO or abs(pj) <= self.ZERO:
            fixation_rate = 0.
        else:
            p_mu = (mu_ji*pj)/(mu_ij*pi)
        
            # If p_mu == 1, L'Hopitals gives fixation rate of 1 (substitution probability is the forward mutation rate) 
            if abs(1. - p_mu) <= self.ZERO:
                fixation_rate = 1. 
            else:
                fixation_rate =  (np.log(p_mu))/(1. - (1./p_mu))
        return fixation_rate * mu_ij            




def codon_fitness_to_freqs(codon_fitness):
    ''' 
        Convert amino acid fitness values to stationary codon frequencies using Sella and Hirsh 2005 (Boltzmann).
        
    '''
 
    codon_freqs = np.zeros(61)
    count = 0
    for codon in codons:
        codon_freqs[count] = np.exp( codon_fitness[codon] )
        count += 1
    codon_freqs /= np.sum(codon_freqs)                   
    assert( abs(1. - np.sum(codon_freqs)) <= ZERO), "codon_freq doesn't sum to 1 in codon_fitness_to_freqs"
    return codon_freqs 



def add_bias_to_fitness(rawfitness, bias):
    '''
        Derive new fitness values which incorporate codon bias.
    '''
    new_fitness = np.zeros(61)
    
    for i in range(len(genetic_code)):
        # Determine the new preferred, non-preferred frequencies
        family = genetic_code[i]
        aa_fit = rawfitness[ codons.index(genetic_code[i][0]) ]
        k = len(family) - 1.
          
        nonpref = abs(aa_fit) * bias * -1 # Reduce fitness by 50-100%
        pref = deepcopy(aa_fit)
        
        # Assign randomly
        indices = [codons.index(x) for x in family]
        shuffle(indices)
        first = True
        for ind in indices:
            if first:
                new_fitness[ind] = pref
                first=False
            else:
                new_fitness[ind] = nonpref
    
    return new_fitness




def add_bias(rawfreqs, bias):
    ''' 
        Derive new codon frequencies which incorporate codon bias.
    '''

    new_freqs = np.zeros(61)

    for i in range(len(genetic_code)):
        # Determine the new preferred, non-preferred frequencies
        family = genetic_code[i]
        aa_freq = rawfreqs[ codons.index(genetic_code[i][0]) ]
        aa_freq_full = aa_freq * len(family)
        k = len(family) - 1.
        
        if k != 0:
            freq_pref = aa_freq_full * bias
            freq_nonpref = (aa_freq_full - freq_pref)/k
            assert( abs(aa_freq_full - (freq_pref + freq_nonpref*k)) <= ZERO), "New bias frequencies improperly calculated."
        else:
            freq_pref = aa_freq
            
        # Assign randomly
        indices = [codons.index(x) for x in family]
        shuffle(indices)
        first = True
        for ind in indices:
            if first:
                new_freqs[ind] = freq_pref
                first=False
            else:
                new_freqs[ind] = freq_nonpref

    assert(abs(1. - np.sum(new_freqs)) <= ZERO), "New bias frequencies don't sum to 1."
    new_freqs_dict = dict(zip(codons, new_freqs))
    
    return new_freqs, new_freqs_dict





def aa_to_codon_freqs(aa_freqs):
    
    codon_freqs = np.zeros(61)
    aa_dict = dict(zip(amino_acids, aa_freqs))
    for aa in aa_dict:
        syn_codons = genetic_code[ amino_acids.index(aa) ]
        cf = aa_dict[aa] / float(len(syn_codons))
        for syn in syn_codons:
            codon_freqs[ codons.index(syn) ] = cf
    
    assert(1. - np.sum(codon_freqs) <= ZERO), "bad codon freqs from aa freqs"
    return codon_freqs, dict(zip(codons, codon_freqs))


def get_eq_from_eig(m):   
    ''' get the equilibrium frequencies from the matrix. the eq freqs are the left eigenvector corresponding to eigenvalue of 0. 
        Code here is largely taken from Bloom. See here - https://github.com/jbloom/phyloExpCM/blob/master/src/submatrix.py, specifically in the fxn StationaryStates
    '''
    (w, v) = linalg.eig(m, left=True, right=False)
    max_i = 0
    max_w = w[max_i]
    for i in range(1, len(w)):
        if w[i] > max_w:
            max_w = w[i]
            max_i = i
    assert( abs(max_w) < ZERO ), "Maximum eigenvalue is not close to zero."
    max_v = v[:,max_i]
    max_v /= np.sum(max_v)
    eq_freqs = max_v.real # these are the stationary frequencies
    
    # SOME SANITY CHECKS
    assert np.allclose(np.zeros(61), np.dot(eq_freqs, m)) # should be true since eigenvalue of zero
    pi_inv = np.diag(1.0 / eq_freqs)
    s = np.dot(m, pi_inv)
    assert np.allclose(m, np.dot(s, np.diag(eq_freqs)), atol=ZERO, rtol=1e-5), "exchangeability and equilibrium does not recover matrix"
    
    # And for some impressive overkill, double check pi_i*q_ij = pi_j*q_ji
    for i in range(61):
        pi_i = eq_freqs[i]
        for j in range(61):
            pi_j = eq_freqs[j]
            forward  = pi_i * m[i][j] 
            backward = pi_j * m[j][i]
            assert(abs(forward - backward) < ZERO), "Detailed balance violated."    
    return eq_freqs




def build_mutsel_matrix(mu_dict, codon_fitness):
    ''' 
        Build the mutsel model matrix according to,  q_ij = mu_ij * Sij / (1-e^(-Sij)).
        The resulting *left* eigenvector of this matrix gives the state frequencies.
    '''
    matrix = np.zeros([61,61])
    for s in range(61):
        for t in range(61):
        
            nucdiff = get_nuc_diff(codons[s], codons[t])
            if len(nucdiff) != 2:
                continue
                
            # Non-diagonal
            sij = codon_fitness[t] - codon_fitness[s]  
            if abs(sij) < ZERO:
                rate = mu_dict[nucdiff]
            else:
                rate = mu_dict[nucdiff] *  (sij)/(1 - np.exp(-1.*sij))
            matrix[s][t] = rate
                
        # Fill in the diagonal position so the row sums to 0, but ensure it doesn't became -0
        matrix[s][s]= -1. * np.sum( matrix[s] )
        if matrix[s][s] == -0.:
            matrix[s][s] = 0.
        assert ( abs(np.sum(matrix[s])) < ZERO ), "Row in instantaneous matrix does not sum to 0."
    return matrix

    



def derive_dnds(codon_freqs_dict, mu_dict = None):
    ''' By default, calculate dS. If no bias and symmetric mutation rates, it will be 1 anyways at virtually no computational cost... '''

    if mu_dict is None:
        mu_dict = {'AC':1.,  'CA':1.,  'AG':1.,  'GA':1.,  'AT':1.,  'TA':1.,  'CG':1.,  'GC':1.,  'CT':1.,  'TC':1.,  'GT':1.,  'TG':1.}

    numer_dn = 0.; denom_dn = 0.;
    numer_ds = 0.; denom_ds = 0.;

    for codon in codon_freqs_dict:
        rate = 0.
        sites = 0.
        rate, sites = calc_paths(codon, codon_freqs_dict, mu_dict, 'nonsyn')
        numer_dn += rate
        denom_dn += sites

        rate = 0.
        sites = 0.
        rate, sites = calc_paths(codon, codon_freqs_dict, mu_dict, 'syn')
        numer_ds += rate
        denom_ds += sites

    assert( denom_dn != 0. and denom_ds != 0.), "dN/dS calc indicates no evolution, maybe????"
    return (numer_dn/denom_dn)/(numer_ds/denom_ds)
    


     
def calc_paths(source, cfreqs, mu_dict, type):
    ''' type is either syn or nonsyn '''
    rate = 0.
    sites = 0.
    source_freq = cfreqs[source]
    for target in codons:
        diff = get_nuc_diff(source, target) # only consider single nucleotide differences since are calculating instantaneous.
        if (type == 'nonsyn' and codon_dict[source] != codon_dict[target]) or (type == 'syn' and codon_dict[source] == codon_dict[target]):
            if len(diff) == 2:
                rate  += calc_subst_prob( source_freq, cfreqs[target], mu_dict[diff], mu_dict[diff[1]+diff[0]] )
                sites += mu_dict[diff]
    rate  *= source_freq
    sites *= source_freq
    return rate, sites



def get_nuc_diff(source, target):
    diff = ''
    for i in range(3):
        if source[i] != target[i]: 
            diff += source[i]+target[i]
    return diff

 
    
def calc_subst_prob(pi, pj, mu_ij, mu_ji):
    if abs(pi) <= ZERO or abs(pj) <= ZERO:
        fixation_rate = 0.
    else:
        p_mu = (mu_ji*pj)/(mu_ij*pi)
        
        # If p_mu == 1, L'Hopitals gives fixation rate of 1 (substitution probability is the forward mutation rate) 
        if abs(1. - p_mu) <= ZERO:
            fixation_rate = 1. 
        else:
            fixation_rate =  (np.log(p_mu))/(1. - (1./p_mu))
    return fixation_rate * mu_ij            










