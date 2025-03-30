#!/usr/bin/env python

# Author: <YOU> <optional@email.address>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.'''

__version__ = "0.1"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set('ATCGNatcgn')
RNA_bases = set('AUCGNaucgn')

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter)-33

def qual_score(phred_score: str) -> float:
    '''Write your own doc string'''
    total=0 #total starts at 0 so we have a number to add the first ph_score to. 
    for char in phred_score: #for each character in the string phred_score
        ph_score = convert_phred(char) #converts ascii to numerical phred score       
        total = total + ph_score #same as total+=ph_score     
    return total / len(phred_score) #returns total sum of all phred scores divided by the length of the string

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)
    # pass

def gc_content(DNA):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used a DNA sequence?"    
    DNA = DNA.upper()
    return (DNA.count("G")+DNA.count("C"))/len(DNA)

#this needs validate base seq in order to work 

def calc_median(lst: list) -> float:
    N = 0
    if len(lst) % 2 != 0: 
        pos = len(lst) // 2
        # for position, number in enumerate lst: 
        median = lst[pos]
        # print(median)
    else: 
        upper = (len(lst) // 2)
        lower = (len(lst) // 2) - 1
        median = (lst[upper] + lst[lower]) / 2
        # median = 
    return(median)
    pass

def oneline_fasta(file):
    # file = args.f
    curr_seq = ""

    with open (file, "r") as FH1, open(file + "one_line.fa", "w") as FH2:
        #open two files, fasta as reading and write to out.fasta
        while True:
            line = FH1.readline()
            # print(line)
            if line == "":
                FH2.write(curr_seq)
                break
            #bc read line finishes w empty srting we use that to break the while true loop
            #BUT we need to make sure it writes out the last line of the file first hence: .write(curr_seq)
            if line[0] == ">": #header
                if len(curr_seq) > 0: #if there is something assigned to curr_seq
                    FH2.write(curr_seq + '\n') #WRITE IT and add a new line (to start the next header)
                FH2.write(line) #write the seq
                curr_seq = "" #reseat curr_seq to empty after we've writen it
            # if line == "":
            #     break
            else:
                curr_seq = curr_seq + line.strip('\n') #if its not a header, strip new lines and assign it to curr_seq
    return()
    pass

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")


#######################
## Deduper functions ##
#######################

### check strandedness function 
def reverse_strand(sam_line: str) -> bool:
    '''This function will access the bitflag of a sam header and return 
    True if reverse, False if foward strand'''
    spline = sam_line.split()
    flag = int(spline[1])
    if flag & 16 == 16:
        return True
    else: 
        return False 


def get_5_start_pos(sam_line: str) -> int:
    '''This function will take access the cigar string and position of a sam header 
    and return the adjusted 5' start position '''
    spline = sam_line.split()
    pos = int(spline[3])
    cigar = spline[5]
    clip_num = str("0")
    rev_strand:bool = reverse_strand(sam_line)
    #if the read is on the + strand, just adjust for left soft clipping
    if rev_strand == False:
        cigar_hit = re.findall(r'(\d+)([A-Z]{1})', cigar)
        # print(matches)
        pos_adj = 0 
        for i, hit in enumerate(cigar_hit):
            #if the first index position of the cigar string is an S 
            if i == 0 and hit[1] == "S":
                #set position adjust = to the integer of soft clipping
                pos_adj += int(hit[0])
            #subtract the soft clipping from the given position to get the true 5' start position 
            new_pos = pos - pos_adj
    # if the match is on the reverse strand 
    if rev_strand == True: 
        pos_adj = 0 
        #create tuple holding the letter and corresponding number 
        cigar_hit = re.findall(r'(\d+)([A-Z]{1})', cigar)
        for hit in cigar_hit:
            #add match number to position adjust 
            if hit[1] == "M":
                pos_adj += int(hit[0])
            #add deletion number to  position adjust
            if hit[1] == "D":
                pos_adj += int(hit[0])
            # adjust for N for deletions 
            if hit[1] == "N":
                pos_adj += int(hit[0])
        for i, hit in enumerate(cigar_hit):
            #if it's 3' clipping, skip
            if i == 0:
                pass
            #if it's 5' clipping, add to the position adjust 
            elif i != 0:
                if hit[1] == "S":
                    pos_adj += int(hit[0])
        #
        new_pos = pos + pos_adj
    return(new_pos)


def get_line_info(sam_line: str) -> tuple:
    '''This function will take in a sam file line and return a touple containing
    [chrom, true 5' start position, strand, UMI]'''
    line_info = ()
    start_pos = get_5_start_pos(sam_line)
    strand = reverse_strand(sam_line)
    spline = sam_line.split()
    UMI = spline[0].split(":")
    UMI = UMI[-1]
    chrom = spline[2]
    line_info = (chrom, start_pos, strand, UMI)
    return(line_info)


###########################
## Demultiplex functions ##
###########################

#function to append the index to the header 
def append_header(record: str, Index_1: str, Index_2: str) -> str:
    '''This function will take in the header of the read list [0], the seq of index 1 and the rev 
    compliment of index 2 [1] and output the header w indexes appended'''
    Index_2 = rev_comp(Index_2)
    return(f"{record} {Index_1}-{Index_2}")

##########################
## Motif-Mark functions ##
##########################

#def find and replace the regex expressions of the
def reg_ex_replace(motif) -> str: 
    '''Takes in motif and outputs seq containing possible variations in regex'''
    motif = motif.upper()
    for char in motif: 
        if char in IUPAC_DICT:
            motif = motif.replace(char, IUPAC_DICT[char])
    return(motif)

#function to find the exon start and length
def find_exon(sequence: str) -> tuple:
    '''Takes a DNA seq and returns the start position and total length'''
    for i in range(len(sequence)):
        if sequence[i].isupper():
            start = int(i)
            break
    for i in range(start, len(sequence)):
        if sequence[i].islower():
            end = int(i)
            break
    length = int(end) - int(start)
    return (start, length)
