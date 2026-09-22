import sys 

def extract_ids_rep(clustered_fasta):
    '''This function extracts from a file the ids of the sequence and returns a list with them'''
    list_reps_ids=[]
    with open(clustered_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                list_reps_ids.append(line.split()[0][1:])
    return list_reps_ids

def remove_redundancy(original_tsv, representative, new_tsv):
    '''This function reads entries from original_tsv and if it is a representative it is written in the new tsv'''
    with open(original_tsv,'r') as tsv1:
        with open(new_tsv, 'w') as tsv2:
            for line in tsv1:
                if line.split('\t')[1] in representative:
                    print(line.rstrip(), file=tsv2)


                
if __name__=="__main__":
    cluster_reps_file=sys.argv[1]
    original_tsv=sys.argv[2]
    new_tsv=sys.argv[3]

    list_ids=extract_ids_rep(cluster_reps_file)
       
    remove_redundancy(original_tsv, list_ids, new_tsv)
