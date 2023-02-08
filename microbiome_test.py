# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 17:32:48 2018

@author: jdkan
"""
from audioop import avg
from typing import *
import numpy
import time
import alignment
import copy
import random
import matplotlib.pyplot as plt

def Load16SFastA(path, fraction = 1.0):
    # from a file, read in all sequences and store them in a dictionary
    # sequences can be randomly ignored for testing by adjusting the fraction
    random.seed(11)
    
    infile = open(path, 'r')
    sequences_16s = {}
    c = 0
    my_seq = ""
    for line in infile:
        if ">" in line:
            my_id = line[1:-1]
            if random.random() < fraction:
                sequences_16s[my_id] = ""
             
        else:
            if my_id in sequences_16s:
                sequences_16s[my_id] += line[:-1]
    
       
    return sequences_16s


def ConvertLibaryToKmerSets(library, K=2) -> Dict[str, Set[str]]:
    import hashlib
    # use the hashlib function to generate unique values to store for each k-mer
    
    #m = hashlib.md5()
    #m.update('TEST'.encode('utf-8'))
    #m.digest() -> returns the hashed value of the updated input
    new_lib = {}
    c = 0
    for k in library.keys():
        new_lib[k] = { library[k][i:(i+K)] for i in range(len(library[k]) - K + 1) }
        
    return new_lib

def JaccardIndex(s1: Set[str], s2: Set[str]) -> float:
    numerator = float(len(s1.intersection(s2)))
    denominator = float(len(s1.union(s2)))
    return numerator/denominator

def KmerMatch(sequence_kmer_set: Set[str], library_kmer_set: Dict[str, Set[str]]) -> Tuple[float, str]:
    best_score = 0.0
    best_match = None
    
    #add your code here to find the best kmer match
    for k in library_kmer_set.keys():
        score = JaccardIndex(sequence_kmer_set, library_kmer_set[k])
        if score > best_score:
            best_score = score
            best_match = k
    
    return best_score, best_match


def AlignmentMatch(sequence, library):
    best_score = -10000000000
    best_match = None

    for k in library.keys():
        score, _, _ = alignment.local_align(sequence, library[k])
        if score > best_score:
            best_score = score
            best_match = k
    
    # or length of sequence
    return best_score / (10 * len(sequence)), best_match

def ThresholdOption(item: Tuple[float, str], threshold: float) -> Optional[str]:
    # print(item, threshold)
    if item[0] > threshold:
        return item[1]
    else:
        return None

ALIGNMENT_THRESHOLD = 0.7

def MeasureAccuracy(alignments: Dict[str, Tuple[float, str]], sequences_16s: Dict[str, str], sample_sequences: Dict[str, str], K: int, threshold: float) -> float:
    kmer_16s_sequences = ConvertLibaryToKmerSets(sequences_16s, K=K)
    kmer_sample_sequences = ConvertLibaryToKmerSets(sample_sequences, K=K)

    reference = { name: ThresholdOption(alignments[name], ALIGNMENT_THRESHOLD) for name in alignments.keys() }
    kmers = { name: ThresholdOption(KmerMatch(kmer_sample_sequences[name], kmer_16s_sequences), threshold) for name in kmer_sample_sequences.keys() }

    # if all(kmers[name] is None for name in kmers.keys()):
    #     print('allnone :)')

    # for name in reference.keys():
    #     print(reference[name] == kmers[name])
    #     print(reference[name])
    #     print(kmers[name])

    # print(reference)
    # print(kmers)

    return sum(1 if reference[name] == kmers[name] else 0 for name in reference.keys()) / len(reference)

def GetPhylum(seq_id: str) -> str:
    return seq_id.split(';')[1]

def GetSampleLocation(sample_name: str) -> str:
    return sample_name.split('_')[0]

if __name__ == "__main__":
    # stuff only to run when not called via 'import' here

    location_dict = {}
    location_dict['C1'] = ("Control", 1)
    location_dict['C2'] = ("Control", 2)
    location_dict['C3'] = ("Control", 3)

    location_dict['S8'] = ("Point-Mon", 1)
    location_dict['S9'] = ("Point-Mon", 2)
    location_dict['S11'] = ("Point-Mon", 3)  ###########

    location_dict['S2'] = ("Point-Allegheny", 1)
    location_dict['S3'] = ("Point-Allegheny", 2)
    location_dict['S4'] = ("Point-Allegheny", 3)

    location_dict['S1'] = ("Sharpsburg", 1)
    location_dict['S15'] = ("Sharpsburg", 2)
    location_dict['S16'] = ("Sharpsburg", 3)
    location_dict['S17'] = ("Sharpsburg", 4)

    location_dict['S12'] = ("Braddock", 1)
    location_dict['S13'] = ("Braddock", 2)
    location_dict['S14'] = ("Braddock", 3)

    location_dict['S5'] = ("Neville Island", 1)
    location_dict['S6'] = ("Neville Island", 2)
    location_dict['S7'] = ("Neville Island", 3)

    random.seed(2261)

    # # SETTINGS FOR ALIGNMENTS
    # fn = "bacterial_16s_genes.fa"
    # sequences_16s = Load16SFastA(fn, fraction = 0.001)
    # fn = "Fall2018CleanReads.fa"
    # sample_sequences = Load16SFastA(fn, fraction = 0.0001)

    # SETTINGS FOR KMERS
    fn = "bacterial_16s_genes.fa"
    sequences_16s = Load16SFastA(fn, fraction = 0.01)
    fn = "Fall2018CleanReads.fa"
    sample_sequences = Load16SFastA(fn, fraction = 0.001)

    print ("Loaded %d 16s sequences." % len(sequences_16s))
    print ("Loaded %d sample sequences." % len(sample_sequences))

    filename = 'alignments.txt'

    # TASK 3
    # # WRITES ALIGNMENT FILE
    # with open(filename, 'w') as f:
    #     for name in sample_sequences.keys():
    #         align_score, align_match = AlignmentMatch(sample_sequences[name], sequences_16s)
    #         f.write(f"{name}, {align_score}, {align_match}\n")
    #         f.flush()
    #     f.close()

    # READS ALIGNMENT FILE
    alignments = {}
    with open(filename, 'r') as f:
        for line in f:
            name, align_score, align_match = line.replace('\n', '').split(', ')
            alignments[name] = (float(align_score), align_match)

    # # GENERATES ACCURACY PLOTS
    # for k in (1, 3, 5, 7, 9, 11):
    #     x = numpy.arange(0.0, 1.0, 0.05)
    #     y = [MeasureAccuracy(alignments, sequences_16s, sample_sequences, k, threshold) for threshold in x]
    #     plt.plot(x, y, label=f'K={k}')
    #     plt.title('Kmer Accuracy Curve')
    #     plt.xlabel('Match Threshold')
    #     plt.ylabel('Accuracy')
    #     plt.legend(loc='best')
    #     plt.savefig(f'plots3/accuracy_k={k}')
    #     plt.cla()

    # # GENERATES TOTAL ACCURACY PLOT
    # for k in (1, 3, 5, 7, 9, 11):
    #     x = numpy.arange(0.0, 1.0, 0.05)
    #     y = [MeasureAccuracy(alignments, sequences_16s, sample_sequences, k, threshold) for threshold in x]
    #     plt.plot(x, y, label=f'K={k}')
    #     plt.title('Kmer Accuracy Curve')
    #     plt.xlabel('Match Threshold')
    #     plt.ylabel('Accuracy')
    #     plt.legend(loc='best')
    # plt.savefig(f'accuracy_all')
    # plt.cla()

    # # TASK 4
    CHOSEN_K=5
    CHOSEN_THRESHOLD=0.35

    kmer_16s_sequences = ConvertLibaryToKmerSets(sequences_16s, K=CHOSEN_K)
    kmer_sample_sequences = ConvertLibaryToKmerSets(sample_sequences, K=CHOSEN_K)

    # location_to_phylum = {}
    # for sample in sample_sequences.keys():
    #     location = GetSampleLocation(sample)
    #     if location not in location_to_phylum:
    #         location_to_phylum[location] = dict()
    #     match = ThresholdOption(KmerMatch(kmer_sample_sequences[sample], kmer_16s_sequences), CHOSEN_THRESHOLD)
    #     if match is not None:
    #         phylum = GetPhylum(match)
    #         if phylum not in location_to_phylum[location]:
    #             location_to_phylum[location][phylum] = 0
    #         location_to_phylum[location][phylum] += 1

    # # print(location_to_phylum)

    # norm_sample_to_phylum = {
    #     sample: {
    #         phylum: location_to_phylum[sample][phylum] / (3 * sum(location_to_phylum[sample].values()))
    #         for phylum in location_to_phylum[sample]
    #     }
    #     for sample in location_to_phylum
    # }
    
    # LOCATIONS = { tup[0] for tup in location_dict.values() }
    # SAMPLE_MAP = { loc: [sample for sample in location_dict.keys() if location_dict[sample][0] == loc] for loc in LOCATIONS }

    # # print(SAMPLE_MAP)

    # loc_phylum_average = {
    #     loc: {
    #         phylum: sum(norm_sample_to_phylum[sample].get(phylum, 0) for sample in SAMPLE_MAP[loc])
    #         for phylum in set().union(*[norm_sample_to_phylum[sample].keys() for sample in SAMPLE_MAP[loc]])
    #     }
    #     for loc in LOCATIONS
    # }

    # # GENERATE PIE CHART 

    # for loc in LOCATIONS:
    #     labels = [phyl[3:] for phyl in loc_phylum_average[loc].keys()]
    #     sizes = [loc_phylum_average[loc][label] for label in loc_phylum_average[loc].keys()]

    #     # print(loc_phylum_average)

    #     fig1, ax1 = plt.subplots()
    #     plt.pie(sizes, labels=labels, autopct='%1.1f%%', radius=2)
    #     plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

    #     plt.title(f'Bacteria Types in {loc}')
    #     plt.savefig(f'pieplots/{loc}')
    #     plt.cla()

    # TASK 5 (9 units -> no function)

    # FIND AN UNMATCHED SEQUENCE

    for sample in sample_sequences.keys():
        match = ThresholdOption(KmerMatch(kmer_sample_sequences[sample], kmer_16s_sequences), CHOSEN_THRESHOLD)
        if match is None:
            print(sample)
            print(sample_sequences[sample])
