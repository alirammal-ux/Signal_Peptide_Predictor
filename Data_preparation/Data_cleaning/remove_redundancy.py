import sys 

def extract_ids_rep(clustered_fasta):
    '''This function extracts from a file the ids of the sequence and returns a set with them'''
    reps_ids=set()   # set: membership test in remove_redundancy is O(1) instead of scanning a list
    with open(clustered_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                reps_ids.add(line.split()[0][1:])
    return reps_ids

def remove_redundancy(original_tsv, representative, new_tsv):
    '''This function reads entries from original_tsv and if it is a representative it is written in the new tsv'''
    with open(original_tsv,'r') as tsv1:
        with open(new_tsv, 'w') as tsv2:
            for line in tsv1:
                if line.split('\t')[1] in representative:
                    print(line.rstrip(), file=tsv2)


                
if __name__=="__main__":
    if len(sys.argv)!=4:
        print('usage: python3 remove_redundancy.py <cluster_rep_seq.fasta> <original.tsv> <output.clean>')
        sys.exit(1)
    cluster_reps_file=sys.argv[1]
    original_tsv=sys.argv[2]
    new_tsv=sys.argv[3]

    rep_ids=extract_ids_rep(cluster_reps_file)

    remove_redundancy(original_tsv, rep_ids, new_tsv)
