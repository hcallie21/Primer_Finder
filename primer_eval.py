#PURPOSE: helper functions for primerfinder.py
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#calculate the GC content of a sequence
def calculate_gc_content(sequence):
    if not sequence:
        return 0.0
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100


def normalize_sequence(seq):
    """Convert BioPython SeqRecord or Seq to a plain string."""
    if isinstance(seq, SeqRecord):
        return str(seq.seq)
    return str(seq)

def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    seq = ''.join(seq) if not isinstance(seq, str) else seq
    complement = str.maketrans("ATCGatcg", "TAGCtagc")
    return seq.translate(complement)[::-1]

def has_self_dimer(primer, check_length=4):
    """
    Check if the reverse complement of the 3' end exists elsewhere in the primer,
    which suggests possible self-dimer formation.
    """
    if len(primer) < check_length * 2:
        return False  # Not enough length to compare meaningfully

    tail = primer[-check_length:]
    rc_tail = reverse_complement(tail)

    # Look for the RC of the tail elsewhere in the sequence
    return rc_tail in primer[:-check_length]



# #contains hairpin structures
# def contains_hairpin(sequence, hairpin_length=3):
#     """Check if the sequence contains a hairpin structure of a specified length."""
#     for i in range(len(sequence) - hairpin_length + 1):
#         subsequence = sequence[i:i + hairpin_length]
#         if( str(subsequence) == str(subsequence)[::-1]): # Check if the subsequence is palindromic
#             return True
#     return False

# #contains self-complementarity
# def contains_self_complementarity(sequence, complement_length=4):
#     """Check if the sequence contains self-complementarity of a specified length."""
#     for i in range(len(sequence) - complement_length + 1):
#         subsequence = sequence[i:i + complement_length]
#         complement = subsequence[::-1].translate(str.maketrans("ATCG", "TAGC"))
#         if subsequence == complement:
#             return True
#     return False

#melting temperature
    #melting temperature from https://www.rosalind.bio/en/knowledge/what-formula-is-used-to-calculate-tm
def calculate_melting_temperature(sequence):
    if not sequence:
        return 0.0
    a_count = sequence.count("A")
    t_count = sequence.count("T")
    g_count = sequence.count("G")
    c_count = sequence.count("C")
    
    #formula for sequences greater than 13 nucleotides in length
    tm =  64.9 + 41 * (g_count + c_count - 16.4) / (a_count + t_count + g_count + c_count)
    return tm

#contains repeats
def contains_repeats(sequence, repeat_length=2):
    sequence = normalize_sequence(sequence)
    """Check if the sequence contains repeated subsequences of a given length."""
    seen = set()
    for i in range(len(sequence) - repeat_length + 1):
        subsequence = sequence[i:i + repeat_length]
        if subsequence in seen:
            return True
        seen.add(subsequence)
    return False


#starts/ends with G/C sequence
def starts_ends_with_gc(sequence):
    """Check if the first or last 4 characters contain G or C."""
    if not sequence:
        return False  # Handle empty string

    upstrm = sequence[:4]     # First 4 characters
    dwnstrm = sequence[-4:]   # Last 4 characters

    return any(base in {'G', 'C'} for base in upstrm) or any(base in {'G', 'C'} for base in dwnstrm)

#determines if a sequence is suitable for primer design
###NOT FINDING SUITABLE PRIMERS
def is_suitable_for_primer(sequence):
    """Check if a sequence meets primer suitability criteria."""

    sequence = str(sequence)  # normalize to string in case it's Seq/SeqRecord

    gc = calculate_gc_content(sequence)
    tm = calculate_melting_temperature(sequence)

    if (
        gc < 40 or gc > 60
        or has_self_dimer(sequence)
        or contains_repeats(sequence)
        or not starts_ends_with_gc(sequence)
        or tm < 45 or tm > 70
    ):
        return False
    return True

    
#returns the best suitable primer sequence
    #upstrm_seq is the candidate sequence to find a suitable primer from
    #min_length and max_length are the minimum and maximum lengths of the primer
    ## If no suitable primer is found, returns None
    ## returns the longest suitable primer found within the specified length range and its index in the sequence
def find_suitable_primer(upstrm_seq, min_length=18, max_length=25):
    print("Calling find_suitable_primer")
    curr_best_primer = None
    curr_best_length = 0

    if min_length > max_length or len(upstrm_seq) < min_length:
        print("Invalid input to suitable primer")
        return None

    for i in range(len(upstrm_seq) - min_length + 1):
        for j in range(min_length, max_length + 1):
            if i + j <= len(upstrm_seq):
                candidate_primer = upstrm_seq[i:i + j]
                if is_suitable_for_primer(candidate_primer):
                    print(f"Found suitable primer: {candidate_primer} at index {i}")
                    if len(candidate_primer) > curr_best_length:
                        curr_best_primer = candidate_primer
                        curr_best_length = len(candidate_primer)
                        if curr_best_length == max_length:
                            return curr_best_primer, i
    return (curr_best_primer, i) if curr_best_primer else None
