from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

def circular_permuted(x):
    """
    Args:
        x (iterator)
    Returns:
        list: All cicular permutations of x
    """
    return([x[i:] + x[:i] for i in range(len(x))])

def self_and_rev_complement(in_dna):
    """
    Args:
        in_dna (string): DNA sequence
    Returns:
        list of strings: The original and reverse complement of in_dna
    """
    all_possible = [in_dna]

    # Get reverse complement
    dna = Seq(in_dna)#, generic_dna)
    rev_complement = str(dna.reverse_complement())
    all_possible.append(rev_complement)
    return(all_possible)

def normalise_str(in_dna):
    """Find all possible eqivalent STR sequences.
    And return the first alphabetically.
    For example, TA = AT. But would return AT.
    """
    all_possible = []
    # Circularly permute original sequence and reverse complement
    for seq in self_and_rev_complement(in_dna):
        for permuted_seq in circular_permuted(seq): # Switch to faster permutation (6)
            all_possible.append(permuted_seq)

    # Sort and take the first
    all_possible.sort()
    return(all_possible[0])

