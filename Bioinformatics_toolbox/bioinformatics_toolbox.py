LICENSE = 'Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.'

print("--"*20, "\r\n"*2, "This python script contains functions useful for bioinfromatics, created by completing end of chapter 8 exercises in the SysMIC python course", "\r\n"*2, "Last modified: 11/04/20", "\r\n"*2, "Copyright (C) 2020 - Maxim Karpov", "\r\n"*2, "License:", LICENSE, "\r\n"*2, "--"*20, "\r\n")


#TASK #1 - function to read a FASTA txt file and produce a single, continuous string that doesn't contain spaces or line breaks

def read_fasta(path):
    '''
    Reads a FASTA file and produces a single,
    continuous string that doesn't contain 
    spaces or line breaks.

    :param path - Path to the FASTA file.
    :type path - str
    :return - Continuous FASTA sequence.
    :type return - str
    '''
    with open(path, mode='r') as fasta_file:
        fasta_data = fasta_file.read()
    
    fasta_list = fasta_data.splitlines()

    sequence_list = fasta_list[1:]
    sequence_joined = ''.join(sequence_list)

    results = sequence_joined
    return results


#TASK 2 - converting DNA sequence to cDNA sequence

def complimentary_dna(path):
    '''
    Produces a complimentary DNA (cDNA) sequence 
    from a FASTA file.
    
    :param path - Path to the FASTA file.
    :type path - str
    :return - Continuous cDNA sequence.
    :type return - str
    '''
    
    dna2cdna = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

    dna_seq = read_fasta(path)
    cdna_seq = str()
    
    for nucleotide in dna_seq:
        cdna_seq += dna2cdna[nucleotide]
    
    return cdna_seq


#TASK 3 - converting cDNA sequence to mRNA 

def transcribe(path):
    '''
    Transcribes a FASTA file sequence to
    an mRNA sequence.

    :param path - Path to the FASTA file.
    :type path - str
    :return - Continuous mRNA sequence.
    :type return - str
    '''

    cdna2mrna = {'T':'A','A':'U','C':'G','G':'C'}

    cdna_seq = complimentary_dna(path)
    mrna_seq = str()

    for nucleotide in cdna_seq:
        mrna_seq += cdna2mrna[nucleotide]

    return mrna_seq


#TASK 4 - converting an mRNA sequence to a peptide sequence at 3 different reading frames to choose from

def translate(path, reading_frame):
    '''
    When given a path to a FASTA file and a request
    for translation of a desired reading frame,
    working in conjunction with previous functions
    (read_fasta, complimentary_dna and transcribe),
    the function completes user's request
    for FASTA sequence translation to protein sequence.
    
    Default argument = print all reading frames.

    :param path - Path to the FASTA file.
    :type path - str
    :param reading_frame - Desired reading frame.
    :type reading_frame - int
    :return - Continuous protein sequence in str or dict formats.
    :type return - dict (default arguement); str (when reading_frame = 1 or 2 or 3)
    '''

    codon2protein = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "CUU": "L",
    "CUC": "L", "CUA": "L", "CUG": "L", "AUU": "I", "AUC": "I",
    "AUA": "I", "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S",
    "AGC": "S", "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCU": "A",
    "GCC": "A", "GCA": "A", "GCG": "A", "UAU": "Y", "UAC": "Y",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAU": "N",
    "AAC": "N", "AAA": "K", "AAG": "K", "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E", "UGU": "C", "UGC": "C", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R",
    "AGG": "R", "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "AUG": "<Met>", "UAA": "<STOP>", "UAG": "<STOP>", "UGA": "<STOP>"
    }

    codon = str()

    if not reading_frame == 1 and not reading_frame == 2 and not reading_frame == 3 :
        reading_frame = (1,2,3)
        prot_seqs = {'Frame 1':'','Frame 2':'','Frame 3':''}
        for frame in reading_frame:
            for nucleotide in range(frame-1, len(transcribe(path))-((len(transcribe(path))-(frame-1))%3), 3):
                codon = transcribe(path)[nucleotide:nucleotide+3]
                prot_seqs['Frame ' + str(frame)] += codon2protein[codon]
        return prot_seqs

    else:
        prot_seq = str()
        for nucleotide in range(reading_frame-1, len(transcribe(path))-((len(transcribe(path))-(reading_frame-1))%3), 3):
            codon = transcribe(path)[nucleotide:nucleotide+3]
            prot_seq += codon2protein[codon]
        return_seq = ('Reading frame ' + str(reading_frame) + ' = ' + str(prot_seq))
        return return_seq


#TASK 5 - a function to translate and print all reading frames of an mRNA sequence

def translate_all_frames(path):
    '''
    When given a path to a FASTA file, working in 
    conjunction with previous functions (read_fasta,
    complimentary_dna and transcribe), the function 
    translates FASTA sequence to protein sequence at
    every possible reading frame.

    :param path - Path to the FASTA file.
    :type path - str.
    :return - Continuous protein sequence at each reading frame.
    :type return - dict.
    '''

    codon2protein = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "CUU": "L",
    "CUC": "L", "CUA": "L", "CUG": "L", "AUU": "I", "AUC": "I",
    "AUA": "I", "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S",
    "AGC": "S", "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCU": "A",
    "GCC": "A", "GCA": "A", "GCG": "A", "UAU": "Y", "UAC": "Y",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAU": "N",
    "AAC": "N", "AAA": "K", "AAG": "K", "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E", "UGU": "C", "UGC": "C", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R",
    "AGG": "R", "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "AUG": "<Met>", "UAA": "<STOP>", "UAG": "<STOP>", "UGA": "<STOP>"
    }

    codon = str()
    reading_frame = (1,2,3)
    prot_seqs = {'Frame 1':'','Frame 2':'','Frame 3':''}

    for frame in reading_frame:
        for nucleotide in range(frame-1, len(transcribe(path))-((len(transcribe(path))-(frame-1))%3), 3):
            codon = transcribe(path)[nucleotide:nucleotide+3]
            prot_seqs['Frame ' + str(frame)] += codon2protein[codon]
    return prot_seqs


#TASK 6 - producing exon sequences at each reading frame

def find_exon(path, frame):
    '''
    Takes an input for a path to a FASTA file
    and for desired reading frame perspective 
    to find all exons in such reading frame.

    :param path - Path to the FASTA file.
    :type path - str
    :param reading_frame - Desired reading frame.
    :type reading_frame - int
    :return - List of all exons in desired reading frame.
    :type return - list
    '''
    
    exons = list()

    sequence = translate_all_frames(path)['Frame ' + str(frame)]
    sequence = str(sequence.replace('<Met>','M'))

    while sequence.count('<STOP>')>0 and sequence.count('M')>0:

            exons.append(sequence[sequence.find('M'):sequence.find('<STOP>')])
            sequence = sequence[sequence.find('<STOP>')+6:]
            sequence = sequence[sequence.find('M'):]
    
    for member in exons:
        if len(member)<2:
            exons.remove(member)

    return exons


def prep_exon_dict_list(path, frame):
    '''
    Given a path to a FASTA file and a desired 
    reading frame number, prepares a list of 
    dictionaries containing sequences of exons
    and their properties: exon number; total 
    exons; length of exon; sequence.

    :param path - Path to the FASTA file.
    :type path - str
    :param reading_frame - Desired reading frame.
    :type reading_frame - int
    :return - List of dictionaries of exon sequences and their properties.
    :type return - list
    '''
    
    exons = find_exon(path,frame)

    frame_exon_dict = dict()
    exon_dict_list = list()
    
    total_exons = len(exons)
    exon_number = int()
    length = int()
    sequence = str()

    for member in exons:
        exon_number = exons.index(member)+1
        length = len(member)
        sequence = member
        frame_exon_dict = {'exon number':exon_number, 
        'total exons':total_exons, 'length':length, 
        'sequence':sequence
        }
        exon_dict_list.append(frame_exon_dict)
    
    return exon_dict_list
        

def extract_exons(path):
    '''
    Given a path to a FASTA file, returns
    a list of dictionaries containing all exons
    and their properties for all reading frames
    of a DNA sequence.

    :param path - Path to the FASTA file.
    :type path - str
    :return - Dictionary of all exons and their properties in all reading frames.
    :type return - dict
    '''

    frame = [1, 2, 3]

    all_exons = dict()

    for number in frame:
        all_exons['frame ' + str(number)] = prep_exon_dict_list(path, number)
    
    return all_exons
