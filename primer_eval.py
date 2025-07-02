#PURPOSE: helper functions for primerfinder.py

#calculate the GC content of a sequence
def calculate_gc_content(sequence):
    if not sequence:
        return 0.0
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100

#see if contains dimers
def contains_dimer(sequence, dimer_length=2):
    for i in range(len(sequence) - dimer_length + 1):
        dimer = sequence[i:i + dimer_length]
        if sequence.count(str(dimer)) > 1:
            return True
    return False

#contains hairpin structures
def contains_hairpin(sequence, hairpin_length=3):
    """Check if the sequence contains a hairpin structure of a specified length."""
    for i in range(len(sequence) - hairpin_length + 1):
        subsequence = sequence[i:i + hairpin_length]
        if( str(subsequence.seq) == str(subsequence.seq)[::-1]): # Check if the subsequence is palindromic
            return True
    return False

#contains self-complementarity
def contains_self_complementarity(sequence, complement_length=4):
    """Check if the sequence contains self-complementarity of a specified length."""
    for i in range(len(sequence) - complement_length + 1):
        subsequence = sequence[i:i + complement_length]
        complement = subsequence[::-1].translate(str.maketrans("ATCG", "TAGC"))
        if subsequence == complement:
            return True
    return False

#melting temperature
    #melting temperature from https://www.rosalind.bio/en/knowledge/what-formula-is-used-to-calculate-tm
def calculate_melting_temperature(sequence):
    if not sequence:
        return 0.0
    a_count = sequence.count("A")
    t_count = sequence.count("T")
    g_count = sequence.count("G")
    c_count = sequence.count("C")
    
    tm =  64.9 + 41 * (g_count + c_count - 16.4) / (a_count + t_count + g_count + c_count)
    return tm

#contains repeats
def contains_repeats(sequence, repeat_length=3):
    """Check if the sequence contains repeats of a specified length."""
    for i in range(len(sequence) - repeat_length + 1):
        subsequence = sequence[i:i + repeat_length]
        if sequence.count(subsequence) > 1:
            return True
    return False

#starts/ends with G/C sequence
def starts_ends_with_gc(sequence):
    """Check if the sequence starts or ends with G or C."""
    if sequence.startswith("G") or sequence.startswith("C") or sequence.endswith("G") or sequence.endswith("C"):
        return True
    return False

#determines if a sequence is suitable for primer design
###NOT FINDING SUITABLE PRIMERS
def is_suitable_for_primer(sequence, min_length=18, max_length=25):
    #GC content between 40-60%
    if (
        calculate_gc_content(sequence) < 40 or calculate_gc_content(sequence) > 60  # GC content NOT between 40-60%
        or contains_dimer(sequence)  # contains dimers
        or contains_hairpin(sequence)  # contains hairpins
        or contains_self_complementarity(sequence)  # contains self-complementarity
        or contains_repeats(sequence)  # contains repeats
        or not starts_ends_with_gc(sequence)  # starts/ends with G/C
        or calculate_melting_temperature(sequence) < 50  # melting temperature < 50C
        or calculate_melting_temperature(sequence) > 68  # melting temperature > 65C
    ):
        return False
    else: 
        return True
    
#returns the best suitable primer sequence
    #upstrm_seq is the candidate sequence to find a suitable primer from
    #min_length and max_length are the minimum and maximum lengths of the primer
    ## If no suitable primer is found, returns None
    ## returns the longest suitable primer found within the specified length range and its index in the sequence
def find_suitable_primer(upstrm_seq, min_length=18, max_length=25):
    curr_best_primer = None
    curr_best_length = 0

    #exception handling for invalid input
    if min_length > max_length or len(upstrm_seq) < min_length:
        return None  # Invalid input

    #loop though the sequence to find suitable primers
    #i is the starting index, j is the length of the candidate primer
    for i in range(len(upstrm_seq) - min_length + 1):
        for j in range(min_length, max_length + 1):
            if i + j <= len(upstrm_seq):
                candidate_primer = upstrm_seq[i:i + j]
                if is_suitable_for_primer(candidate_primer, min_length, max_length):
                    if len(candidate_primer) > curr_best_length:
                        curr_best_primer = candidate_primer
                        curr_best_length = len(candidate_primer)
                        if curr_best_length == max_length:
                            return curr_best_primer, i  # Early exit
    if curr_best_primer is None:
        return None
    return curr_best_primer, i
    
