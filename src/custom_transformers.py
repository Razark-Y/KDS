from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.decomposition import PCA
from itertools import product
import pandas as pd

codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S', 'CCT': 'P', 'CCC': 'P',
    'CCA': 'P', 'CCG': 'P', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
    'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*'
}

def find_all_orfs_fast(dna, codon_table, min_aa_length=30):
    dna = dna.upper()
    stop_codons = {'TAA', 'TAG', 'TGA'}
    orfs = []
    for frame in range(3):
        i = frame
        while i <= len(dna) - 3:
            codon = dna[i:i+3]
            if codon == 'ATG':
                protein = []
                j = i
                while j <= len(dna) - 3:
                    next_codon = dna[j:j+3]
                    if next_codon in stop_codons:
                        break
                    aa = codon_table.get(next_codon)
                    if aa is None:
                        break
                    protein.append(aa)
                    j += 3
                if len(protein) >= min_aa_length:
                    orfs.append(''.join(protein))
                i = j
            else:
                i += 3
    return orfs
def gc_content(genome):
    genome = genome.upper()
    g = genome.count("G")
    c = genome.count("C")
    a = genome.count("A")
    t = genome.count("T")
    total = a + t + g + c
    return 100 * (g + c) / total if total > 0 else 0

def count_differences(seq1, seq2):
    return sum(1 for a, b in zip(seq1.upper(), seq2.upper()) if a != b)
class NRemover(BaseEstimator, TransformerMixin):
    def fit(self, X, y=None):
        return self
    def transform(self, X):
        X = X.copy()
        X['FASTA'] = X['FASTA'].str.upper().str.replace('N', '', regex=False)
        X['length'] = X['FASTA'].str.len()
        return X

class TwoMerPCA(BaseEstimator, TransformerMixin):
    def __init__(self, n_components=10):
        self.n_components = n_components
        self.kmers = [''.join(p) for p in product('ACGT', repeat=2)]
        self.pca = PCA(n_components=n_components)

    def fit(self, X, y=None):
        kmer_df = X['FASTA'].apply(self._extract_kmer_counts)
        self.pca.fit(kmer_df.values)
        return self

    def transform(self, X):
        X = X.copy()
        kmer_df = X['FASTA'].apply(self._extract_kmer_counts)
        reduced = self.pca.transform(kmer_df.values)
        for i in range(self.n_components):
            X[f'2mer_pca_{i+1}'] = reduced[:, i]
        return X

    def _extract_kmer_counts(self, seq):
        seq = seq.upper()
        counts = dict.fromkeys(self.kmers, 0)
        for i in range(len(seq) - 1):  
            kmer = seq[i:i+2]
            if set(kmer).issubset({'A', 'C', 'G', 'T'}):  
                counts[kmer] += 1
        return pd.Series(counts)

class ThreeMerPCA(BaseEstimator, TransformerMixin):
    def __init__(self, n_components=10):
        self.n_components = n_components
        self.kmers = [''.join(p) for p in product('ACGT', repeat=3)]
        self.pca = PCA(n_components=n_components)
    def fit(self, X, y=None):
        kmer_df = X['FASTA'].apply(self._extract_kmer_counts)
        self.pca.fit(kmer_df.values)
        return self
    def transform(self, X):
        X = X.copy()
        kmer_df = X['FASTA'].apply(self._extract_kmer_counts)
        reduced = self.pca.transform(kmer_df.values)
        for i in range(self.n_components):
            X[f'3mer_pca_{i+1}'] = reduced[:, i]
        return X

    def _extract_kmer_counts(self, seq):
        seq = seq.upper()
        counts = dict.fromkeys(self.kmers, 0)
        for i in range(len(seq) - 2):
            kmer = seq[i:i+3]
            if kmer in counts:
                counts[kmer] += 1
        return pd.Series(counts)

class FeatureDropper(BaseEstimator, TransformerMixin):
    def __init__(self, columns=None):
        self.columns = columns if columns is not None else []

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        return X.drop(columns=self.columns, errors='ignore')

class FeatureCreator(BaseEstimator, TransformerMixin):
    def __init__(self):
        with open("src/Sars.txt", 'r') as f:  
            self.sars_sequence = f.read().strip()

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        X = X.copy()
        base_counts = X['FASTA'].apply(lambda seq: pd.Series({
            'A': seq.upper().count('A'),
            'C': seq.upper().count('C'),
            'G': seq.upper().count('G'),
            'T': seq.upper().count('T')
        }))
        X = pd.concat([X, base_counts], axis=1)
        X['gc_content'] = X['FASTA'].apply(gc_content)
        X['mutation'] = X['FASTA'].apply(
            lambda seq: count_differences(seq, self.sars_sequence)
        )
        return X
