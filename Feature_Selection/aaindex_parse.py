def parse_aaindex_file(filepath, record_id):
    """Extracts the 20 numerical values for a given ID from an aaindex file."""
    with open(filepath, 'r') as file:
        found_record = False
        
        for line in file:
            # 1. Find the specific record
            if line.startswith(f"H {record_id}"):
                found_record = True
                
            # 2. Once inside the record, look for the 'I' line
            if found_record and line.startswith("I "):
                # The standard amino acid order in AAindex1
                aa_order = "ARNDCQEGHILKMFPSTWYV"
                
                # The values are on the two lines immediately after the 'I' line
                val_line_1 = file.readline().strip().split()
                val_line_2 = file.readline().strip().split()
                
                values_str = val_line_1 + val_line_2
                
                # Convert 'NA' to 0.0, and parse the rest as floats
                values = [float(v) if v != 'NA' else 0.0 for v in values_str]
                
                return dict(zip(aa_order, values))
                
            # 3. Stop searching if we hit the end of the current record
            if found_record and line.startswith("//"):
                break
                
    return None



if __name__=='__main__':
# Example usage:
    from pathlib import Path
    result = parse_aaindex_file(Path(__file__).resolve().parent / "aaindex1", "BIGC670101")
    print(result)

