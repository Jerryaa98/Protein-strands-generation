import csv
if __name__ == "__main__":
    
    database_file = "/root/Biology_project/data/OMBB_data.csv"
    pdb_dir = '/root/Biology_project/pdb_files'

    with open(database_file, 'r') as file1:
        reader = csv.reader(file1)
        next(reader)  # Skip the header
        for row1 in reader:
            seq_id = row1[0]
            seq_len = int(row1[-1])
            seq = row1[2]
            assert len(seq) == seq_len, 'Erorr'
            pdb_file = f'{pdb_dir}/{seq_id}.pdb'
            CA_count = 0
            with open(pdb_file, 'r') as file2:
                for row2 in file2:
                    if row2.startswith("ATOM") and "CA" in row2: 
                        CA_count += 1
            if seq_len != CA_count:
                print(f'Id : {seq_id} | SEQ LEN : {seq_len} | CA Count : {CA_count}')
           

