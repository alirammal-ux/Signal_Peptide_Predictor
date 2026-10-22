import requests
from requests.adapters import HTTPAdapter, Retry
import json
import re
from sys import argv

##############################
# FUNCTIONS for POSITIVE SET
##############################

def valid_signal(entry):
  '''This function returns the signal peptide feature that starts at the beginning of the chain and is
   longer or equal to 14 amino acids, or None if there is none or if it is not cleaved'''
  for i in entry['features']:
    if i['type']=='Signal':
      if i['location']['start']['value']==1:
        if i['description']=='Not Cleaved': return None
        elif i['location']['end']['value']!='?' and i['location']['end']['value']!=None and i['location']['end']['value']>13:
          return i

def filter_function_positive(entry):
  '''This function returns true if the entry has a signal peptide that starts at the
   beginning of the chain and is longer or equal to 14 amino acids, excluding not cleaved ones'''
  return valid_signal(entry) is not None

def extract_fields_positive(entry):
    UPaccession=entry['primaryAccession']
    org_name=entry['organism']['scientificName']
    seq_len=entry['sequence']['length']
    sequence = entry.get("sequence", {}).get("value", "")
    class_= 1
    early_trans = False


    # kingdom: here we classify the organism into one of the three kingdoms: Metazoa, Fungi, or Plants.
    #  If the organism does not belong to any of these kingdoms, it is classified as "Other".
    #  We do this expecting fot the other kingdoms to be way less than the three main ones. --> to be checked during data statistics.
    lineage = entry.get("organism", {}).get("lineage", [])
    if "Metazoa" in lineage:
        kingdom = "Metazoa"
    elif "Fungi" in lineage:
        kingdom = "Fungi"
    elif "Viridiplantae" in lineage:
        kingdom = "Plants"
    else:
        kingdom = "Other"

    # cleavage site = end of the same signal peptide accepted by filter_function_positive
    cleavage_pos = valid_signal(entry)['location']['end']['value']

    return [class_, UPaccession, org_name, kingdom, seq_len, cleavage_pos, early_trans, sequence]



############################
# FUNCTIONS for NEGATIVE SET
############################

def extract_fields_negative(entry):

    accession = entry.get("primaryAccession", "")
    organism_name = entry.get("organism", {}).get("scientificName", "")
    length = entry.get("sequence", {}).get("length", 0)
    sequence = entry.get("sequence", {}).get("value", "")
    class_ = 0
    cleavage_pos = 0


    #kindom
    lineage = entry.get("organism", {}).get("lineage", [])
    if "Metazoa" in lineage:
        kingdom = "Metazoa"
    elif "Fungi" in lineage:
        kingdom = "Fungi"
    elif "Viridiplantae" in lineage:
        kingdom = "Plants"
    else:
        kingdom = "Other"

    early_trans = False

    for f in entry.get("features", []):
        if f.get("type") == "Transmembrane" and f.get("description") == "Helical":
            try:
                start_val = f["location"]["start"]["value"]
                if start_val < 90:
                    early_trans = True
                    break
            except KeyError:

                continue

    return (class_, accession, organism_name, kingdom, length, cleavage_pos, early_trans, sequence)



######################
# COMMON PART
######################

retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

# Function that extract from the header the cursor/token for the next batch to be retrived and generate the url for the next request
def get_next_link(headers):
    if "Link" in headers:
        # The regular expression is used to extract the next link for pagination
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

# This function actually retrieve the next data batch in the search.
# The function act as an iterator, yielding the next result batch at every call
# The function terminates after the last batch has been returned. In this case,
# the next link will be None
def get_batch(batch_url):
    while batch_url:
        # Run the API call
        response = session.get(batch_url)
        # Will raise an error if an error status code is obtained
        response.raise_for_status()
        # Get the total number of entries in the search
        total = response.headers["x-total-results"]
        # Yield the response and the total number of entries
        yield response, total #--> GENERATOR: obtain the data of current batch and yield them
        # Get the link to the API call for the next data batch
        batch_url = get_next_link(response.headers)

def get_dataset(search_url, extract_function, output_file_name, *filter_function):
    # If we specify a function as argument of another function means that we can call that function by using the name of the argument
    filtered_json = []
    n_total, n_filtered = 0, 0
    pos = False
    if filter_function:
        pos = True
        filter_function = filter_function[0]

    # Run the API call in batches
    for batch, total in get_batch(search_url):
        # parse the JSON body of the response
        batch_json = json.loads(batch.text)
        # filter the entries
        for entry in batch_json["results"]:
            n_total += 1
            # POSITIVE DATASET
            if pos:
              # Check if the entry passes the filter (if not ignore it)
              if filter_function(entry):
                  n_filtered += 1
                  filtered_json.append(extract_function(entry))
              '''
              ## JUST FOR TEST ##
              else:
                print(extract_function(entry))
              '''
            # NEGATIVE DATASET
            else:
              n_filtered += 1
              filtered_json.append(extract_function(entry))
              
    print(n_total, n_filtered)

    tsv = f"{output_file_name}.tsv"
    fasta = f"{output_file_name}.fasta"

    with open(tsv, "w") as ofs:
        for entry in filtered_json:
            # Print the fields in TSV format
            print(*entry[:7], sep="\t", file=ofs)

    with open(fasta, "w") as ofs:
        for entry in filtered_json:
            # Print the fields in FASTA format
            print(f">{entry[1]}", file=ofs)
            print(entry[-1], file=ofs)




if __name__ == "__main__":
    batch_size = 500


    if len(argv)==3 and argv[1] in ('-pos','-neg'):  
        if argv[1]=='-pos':
            url_positive = "https://rest.uniprot.org/uniprotkb/search?compressed=false&format=json&query=%28%28ft_signal_exp%3A*%29+AND+%28reviewed%3Atrue%29+AND+%28fragment%3Afalse%29+AND+%28taxonomy_id%3A2759%29+AND+%28length%3A%5B40+TO+*%5D%29+AND+%28existence%3A1%29%29&size=500"
            output_file_positive = argv[2]
            get_dataset(url_positive, extract_fields_positive, output_file_positive, filter_function_positive)
        else:
            url_negative = "https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28taxonomy_id%3A2759%29+AND+%28existence%3A1%29+AND+%28reviewed%3Atrue%29+AND+%28fragment%3Afalse%29+AND+%28length%3A%5B40+TO+*%5D%29+AND+%28%28cc_scl_term_exp%3ASL-0091%29+OR+%28cc_scl_term_exp%3ASL-0191%29+OR+%28cc_scl_term_exp%3ASL-0173%29+OR+%28cc_scl_term_exp%3ASL-0039%29+OR+%28cc_scl_term_exp%3ASL-0209%29+OR+%28cc_scl_term_exp%3ASL-0204%29%29+NOT+%28ft_signal%3A*%29+NOT+%28%28cc_scl_term%3ASL-0095%29+OR+%28cc_scl_term%3ASL-0132%29+OR+%28cc_scl_term%3ASL-0112%29+OR+%28cc_scl_term%3ASL-0243%29%29%29"
            output_file_negative = argv[2]
            get_dataset(url_negative, extract_fields_negative, output_file_negative)

    elif len(argv)==4 and argv[1]=='-all':
        url_negative = "https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28taxonomy_id%3A2759%29+AND+%28existence%3A1%29+AND+%28reviewed%3Atrue%29+AND+%28fragment%3Afalse%29+AND+%28length%3A%5B40+TO+*%5D%29+AND+%28%28cc_scl_term_exp%3ASL-0091%29+OR+%28cc_scl_term_exp%3ASL-0191%29+OR+%28cc_scl_term_exp%3ASL-0173%29+OR+%28cc_scl_term_exp%3ASL-0039%29+OR+%28cc_scl_term_exp%3ASL-0209%29+OR+%28cc_scl_term_exp%3ASL-0204%29%29+NOT+%28ft_signal%3A*%29+NOT+%28%28cc_scl_term%3ASL-0095%29+OR+%28cc_scl_term%3ASL-0132%29+OR+%28cc_scl_term%3ASL-0112%29+OR+%28cc_scl_term%3ASL-0243%29%29%29"
        output_file_negative = argv[2]
        get_dataset(url_negative, extract_fields_negative, output_file_negative)

        url_positive = "https://rest.uniprot.org/uniprotkb/search?compressed=false&format=json&query=%28%28ft_signal_exp%3A*%29+AND+%28reviewed%3Atrue%29+AND+%28fragment%3Afalse%29+AND+%28taxonomy_id%3A2759%29+AND+%28length%3A%5B40+TO+*%5D%29+AND+%28existence%3A1%29%29&size=500"
        output_file_positive = argv[3]
        get_dataset(url_positive, extract_fields_positive, output_file_positive, filter_function_positive)

    else:
        print('usage: python3 data_collection.py -pos <output prefix>\n'
              '       python3 data_collection.py -neg <output prefix>\n'
              '       python3 data_collection.py -all <output prefix negative> <output prefix positive>\n'
              'each prefix produces <prefix>.tsv and <prefix>.fasta')
        exit()