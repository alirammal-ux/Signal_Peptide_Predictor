# OBSOLETE: older version kept for reference, use Data_collection/data_collection.py instead.
# This script writes a 5-column TSV (no class column), while remove_redundancy.py and datasplitting.py
# expect the 7-column TSV written by data_collection.py.

import requests
from requests.adapters import HTTPAdapter, Retry
import json
import re
import pprint
from pathlib import Path

pp = pprint.PrettyPrinter()

#allow retries for UniProt API calls in case of errors 500, 502, 503, 504

retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
# force the retrival to recall the api if errors 500, 502, 503, 504 happens
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

# NB: the last batch has no cursor so the last url will be empty

def valid_signal(entry):
  '''This function returns the signal peptide feature that starts at the beginning of the chain and is
   longer or equal to 14 amino acids, or None if there is none or if it is not cleaved'''
  for i in entry['features']:
    if i['type']=='Signal':
      if i['location']['start']['value']==1:
        if i['description']=='Not Cleaved': return None
        elif i['location']['end']['value']!='?' and i['location']['end']['value']!=None and i['location']['end']['value']>13:
          return i

def filter_function(entry):
  '''This function returns true if the entry has a signal peptide that starts at the
   beginning of the chain and is longer or equal to 14 amino acids, excluding not cleaved ones'''
  return valid_signal(entry) is not None

def extract_fields(entry):
    UPaccession=entry['primaryAccession']
    org_name=entry['organism']['scientificName']
    seq_len=entry['sequence']['length']
    # cleavage site = end of the same signal peptide accepted by filter_function
    cleavage_pos = valid_signal(entry)['location']['end']['value']

    lineage = entry.get("organism", {}).get("lineage", [])
    if "Metazoa" in lineage:
        kingdom = "Metazoa"
    elif "Fungi" in lineage:
        kingdom = "Fungi"
    elif "Viridiplantae" in lineage:
        kingdom = "Plants"
    else:
        kingdom = "Other"

    return [UPaccession, org_name, kingdom, seq_len, cleavage_pos]

def extract_fasta(entry):
    UPaccession=entry['primaryAccession']
    fasta_seq=entry['sequence']['value']
    return [f'>{UPaccession}\n{fasta_seq}']

def get_dataset(search_url, filter_function, extract_function,extract_fasta, output_file_name):
    # If we specify a function as argument of another function means that we can call that function by using the name of the argument
    filtered_json = []
    fasta_file=[]
    n_total, n_filtered = 0, 0
    # Run the API call in batches
    for batch, total in get_batch(search_url):
        # parse the JSON body of the response
        batch_json = json.loads(batch.text)
        # filter the entries
        for entry in batch_json["results"]:
            n_total += 1
            # Check if the entry passes the filter (if not ignore it)
            if filter_function(entry):
                n_filtered += 1
                filtered_json.append(extract_function(entry))
                fasta_file.append(extract_fasta(entry))
            # ## JUST FOR TEST ##
            # else:
            #   print(extract_function(entry))

    print(n_total, n_filtered)
    output_file_name_tsv = output_file_name+".tsv"
    output_file_name_fasta= output_file_name+".fasta"
    with open(output_file_name_tsv, "w") as ofs:
        for entry in filtered_json:
            # Extract the fields of interest
            #fields = extract_function(entry)
            # Print the fields in TSV format
            print(*entry, sep="\t", file=ofs)
        ofs.close
    with open(output_file_name_fasta, "w") as ofs:
            for entry in fasta_file:
                print(*entry, sep="\n", file=ofs)
            ofs.close

if __name__=="__main__":
    batch_size = 500
    url = "https://rest.uniprot.org/uniprotkb/search?compressed=false&format=json&query=%28%28ft_signal_exp%3A*%29+AND+%28reviewed%3Atrue%29+AND+%28fragment%3Afalse%29+AND+%28taxonomy_id%3A2759%29+AND+%28length%3A%5B40+TO+*%5D%29+AND+%28existence%3A1%29%29&size=500"
    # We set the name of the output file, we want TSV output
    output_file = str(Path(__file__).resolve().parent / "positive_set")
    # compressed=false --> in order to be able to manage the response
    # format=json --> (NOT fasta)   
    get_dataset(url, filter_function, extract_fields, extract_fasta, output_file)

  