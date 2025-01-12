# Import modules from different libraries:

import os
import re
import csv
from glob import glob
import pandas as pd
from Bio import SeqIO
from Bio import Entrez
from collections import defaultdict

# Define directories (input and output files for the different functions in this script)

GFFs = '/home/alumno/Compartida24/TFM/000_TFM/GFFs' # Input directory (it includes GFF files with the genomic coordinates of 18S and 28S)
ITS_coordinates = '/home/alumno/Compartida24/TFM/000_TFM/ITS_coordinates.csv' # Output file with coordinates of ITS
sequences_dir = '/home/alumno/Compartida24/TFM/000_TFM/Sequences_ITS'  # Directory to save the sequences retrieved from the National Center for Biotechnology Information
output_FASTAfile = '/home/alumno/Compartida24/TFM/000_TFM/Sequences_ITS/Allseqs.fasta' # Fasta file including all ITS sequences retrieved
taxonomy_file = '/home/alumno/Compartida24/TFM/000_TFM/fungi_taxid_taxonomy.txt'  # Text file with taxonomy data
taxonomy_fasta = '/home/alumno/Compartida24/TFM/000_TFM/Sequences_ITS/Allseqs_taxonomy.fasta'  # Annotated fasta file
output_table = '/home/alumno/Compartida24/TFM/000_TFM/Sequences_ITS/Allseqs_taxonomy_Shared_seqs.txt'  # Table file with info on the shared sequences among different species
Uniqueseqs_taxonomy = '/home/alumno/Compartida24/TFM/000_TFM/Sequences_ITS/UniqueSeqs_taxonomy.fasta' # Fasta file with non identical sequences (obtained from the fasta with all seqs)
Subseq_table = '/home/alumno/Compartida24/TFM/000_TFM/Sequences_ITS/Subsequences.txt'  # Table file with info on the subsequences
Nosubseq_fasta = '/home/alumno/Compartida24/TFM/000_TFM/Sequences_ITS/UniqueSeqs_taxonomy_FINAL.fasta'  # Output FASTA file with no subsequences (FINAL FASTA FILE)



###############################################################################################################################################################################################################################
#################################################################################################### EXTRACT ITS COORDINATES ################################################################################################## 
###############################################################################################################################################################################################################################

def extract_ITS(GFFs, ITS_coordinates): # Create the function "extract_ITS" to extract the coordinates of ITS from each GFF file
    
    gff_files = glob(os.path.join(GFFs, '*.gff')) # Find all GFFs files in my input directory
    
    # Create an empty list to store the results of the function
    results_positive = []  # List for the results obtained for positive strands
    results_negative = []  # List for negative strands
    
    # For each file in my input directory
    for file in gff_files:
        print(f"\nProcessing genome {'_'.join(file.split('/')[-1].split('_')[:2])}") # Print "Processing genome" followed by the name of the genome

        # Create a dataframe (df) by reading the GFF files. They are text files in CSV delimited by tab (\t)      
        df = pd.read_csv(file, sep='\t', header=None,
                         names=['contig_name', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute', '']) 
        df = df.dropna(axis=1, how='all') # Delete the columns with only NaN values
        df = df.sort_values(by=['strand', 'start'], ascending=[True, True]) # Order the dataframe by the columns 'strand' and 'start'

        
        df['attribute'] = df['attribute'].astype(str).fillna('') # Transform the column 'attribute' to string   

        # Create a separate dataframe to store the data for each gene 18S and 28S
        ssu_df = df[df['attribute'].str.contains('18s_rRNA', na=False, case=False)] 
        lsu_df = df[df['attribute'].str.contains('28s_rRNA', na=False, case=False)]

        # Check if there are records for 18S (SSU) and 28S (LSU)
        if ssu_df.empty or lsu_df.empty:
            print("There are no records of SSU nor LSU in this file.")
            continue

        # Print the number of records for 18S and 28S found in each GFF
        #print(f"Number of 18S sequences found: {len(ssu_df)}")
        #print(f"Number of 28S sequences found: {len(lsu_df)}")

        # Process pairs of strands
        for i in range(len(df) - 1): # for each row in the range of the length of my dataframe, except the last one...
            current_row = df.iloc[i] # "current_row" is the current row
            next_row = df.iloc[i + 1] # "next_row" is the current row plus 1

            if current_row['strand'] == next_row['strand']: #  Check if both strands are positive or negative (the strands MUST have the same orientation to find the coordinates of ITS in a correct way)
                if current_row['strand'] == '+' and ('18s_rRNA' in current_row['attribute'] and '28s_rRNA' in next_row['attribute']): # IF THE STRANDS ARE POSITIVE (18S is located before 28S)...
                    start_18S = current_row['start']
                    end_18S = current_row['end']
                    start_28S = next_row['start']
                    end_28S = next_row['end']

                    
                    if 100 <= (start_28S - end_18S) <= 250: # If the distance between 18S and 28S is >= 100 and <= 250 bp...
                        results_positive.append({ # Add info to the list "results_positive"
                            'contig_name': current_row['contig_name'],
                            'start_18S': start_18S,
                            'end_18S': end_18S,
                            'start_28S': start_28S,
                            'end_28S': end_28S,
                            'ITS_start': end_18S + 1,
                            'ITS_end': start_28S - 1,
                            'strand': current_row['strand'],
                            'score': current_row['score'],
                            'genome_name': '_'.join(file.split('/')[-1].split('_')[:2]),
                            'ITS_length': (start_28S - 1) - (end_18S + 1)
                        })
                        print(f"Pair 18S-28S found in the positive strand: {results_positive[-1]}")
                        i += 2  # The counter i is increase by 2, so the program does not go the immediately next row, but the following one

                
                elif current_row['strand'] == '-' and ('28s_rRNA' in current_row['attribute'] and '18s_rRNA' in next_row['attribute']): # IF THE STRANDS ARE NEGATIVE (28S is located before 18S)...
                    start_18S = next_row['start']
                    end_18S = next_row['end']
                    start_28S = current_row['start']
                    end_28S = current_row['end']

                    if 100 <= (start_18S - end_28S) <= 250: # If the distance between the start of 18S and the end of 28S is >= 100 and <= 250 bp...
                        results_negative.append({ # Add the following info to the list "results_negative"
                            'contig_name': current_row['contig_name'],
                            'start_18S': start_18S,
                            'end_18S': end_18S,
                            'start_28S': start_28S,
                            'end_28S': end_28S,
                            'ITS_start': end_28S + 1,
                            'ITS_end': start_18S - 1,
                            'strand': current_row['strand'],
                            'score': current_row['score'],
                            'genome_name': '_'.join(file.split('/')[-1].split('_')[:2]),
                            'ITS_length': (start_18S - 1) - (end_28S + 1)
                        })
                        print(f"Pair 28S-18S found in the negative strand: {results_negative[-1]}")
                        i += 2  # The counter i is increase by 2, so the program does not go the immediately next row, but the following one
                else:
                    print(f"The 18S or 28S sequence found has no partner to be combined with")


    if results_positive or results_negative: # If both lists are not empty, combine them into a single dataframe
        results_combined = results_positive + results_negative
        results_df = pd.DataFrame(results_combined)
        results_df.to_csv(ITS_coordinates, index=False)
        print("\nITS coordinates extraction done. Data saved in file", ITS_coordinates)
    else:
        print("\nWarning message: no valid pairs 18S/28S (+strand) or 28S/18S (-strand) were found in the gff files processed.") # If one of the list is empty, give me a warning message

# Call the funtion "extract_ITS". We obtain a csv file (ITS_coordinates.csv) with the coordinates of ITS. The next step is to get the sequences corresponding to those coordinates using the module Entrez.
extract_ITS(GFFs, ITS_coordinates)


###############################################################################################################################################################################################################################
#################################################################################################### GETTING ITS SEQUENCES USING Entrez ####################################################################################### 
###############################################################################################################################################################################################################################
Entrez.email = "kina@usal.es" # Configure Entrez (replace the email account!!!)

def obtain_seq(accession, start, end, strand): # Define the function "obtain_seq", which receives four parameters (accession, start position of the sequence to be extracted, end position and strand)
    #print(f"Getting ITS from contig {accession} (positions {start}-{end}; {strand} strand)")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") # db="nucleotide" indicates that the data will be obtained from the nucleotide database; Use rettype="fasta" to obtain the sequence in fasta file and retmode="text" to get the result as text
    record = SeqIO.read(handle, "fasta") # Read the seq in fasta format and give it back as an object (record)
    handle.close() # Close the object handle after reading the seq
    sequence = record.seq[start-1:end]  # Extract the subsequence limited by the positions start and end from the seq downloaded by Entrez 
    if strand == '-':  # If the strand is negative, we extract the reverse complementary seqeunce
        sequence = sequence.reverse_complement()
    return sequence
  
def process_csv(ITS_coordinates, sequences_dir): # Define the function "process_csv" (remember! the input file for this function is the csv file with the coordinates "ITS_coordinates.csv"S)
    os.makedirs(sequences_dir, exist_ok=True)    # Create a new output directory
    with open(ITS_coordinates, newline='') as csvfile: # Read the csv
        reader = csv.DictReader(csvfile)
        for row in reader:
            accession = row['contig_name']
            start = int(float(row['ITS_start']))
            end = int(float(row['ITS_end']))     
            strand = row['strand']
            genome_name = '_'.join(row['genome_name'].split('/')[-1].split('_')[:2]) # Extract the value of the name of the genome from "genome_name" using a tab as delimiter and keeping only the last part. Then divide it using "_" as delimiter and join together the two first fragment
            #genomename = row['genome_name'] # This is basically the same, but in several steps.
            #genoma = genomename.split('/')[-1]
            #genoma2 = genoma.split('_')[:2]
            #genome_name = '_'.join(genoma2)

            sequence = obtain_seq(accession, start, end, strand=strand) # Get the ITS sequence
            if sequence: # if a sequence exists...
                filename = f"{sequences_dir}/{genome_name}_{accession}_{start}_{end}.fasta" # Save it as fasta
                with open(filename, "w") as output_FASTAfile:
                    output_FASTAfile.write(f">{genome_name}|{accession}_{start}_{end}\n{sequence}\n")
                #print(f"The sequence has been saved as {filename.split('/')[1]} in the directory {sequences_dir}\n")
            else:
                print(f"Error, no ITS seq was found for the accession {accession}")

# Call the function "process_csv", with "ITS_coordinates" as input file and "sequences_dir" as output directory
process_csv(ITS_coordinates, sequences_dir)

###############################################################################################################################################################################################################################
########################################################################### COMBINE ALL THE ITS SEQUENCES INTO A SINGLE FASTA FILE ############################################################################################ 
###############################################################################################################################################################################################################################

with open(output_FASTAfile, 'a') as outfile: # Open the fasta file with the mode 'a' (append) to avoid overwriting the file 
    for filename in os.listdir(sequences_dir): # For each file in my directory of sequences...
        with open(os.path.join(sequences_dir, filename), 'r') as infile:
            outfile.write(infile.read()) # Write the content in the output file
print(f"All ITS sequences have been combined into a fasta file ({output_FASTAfile})")


###############################################################################################################################################################################################################################
#################################################### ADD TAXONOMIC INFORMATION TO ALL SEQUENCES (SO IT CAN BE CHECKED IF DIFFERENT SPECIES SHARE THE SAME SEQ) ################################################################
###############################################################################################################################################################################################################################


sequences = []  # Create an empty list to store the information of the sequences (name, sequence and contig/positions)

with open(output_FASTAfile, 'r') as file: # Open and read the fasta with all seqs and associate names and sequences
    nombre = None
    contig_pos = None
    sequence = ""

    for line in file: # For each line in my file...
        if line.startswith(">"):  # if the line starts with ">" ...
            if nombre:  # If there is a name of a sequence...
                sequences.append({"nombre": nombre, "contig_pos": contig_pos, "sequence": sequence})  # Add the information (genome_ID, contig positios and sequence to the list "sequence"
            nombre = line.strip().split('|')[0][1:]  # Extract the name of the genome (part of the line between ">" and "|")
            contig_pos = line.strip().split('|')[1].split(',')[0] # Extract contig and positions
            sequence = ""  # Restart the sequence for each new contig
        else:
            sequence += line.strip()  # If the lines do not start with ">", they are not sequences and can be concatenated
    # Add the last sequence
    if nombre:
        sequences.append({"nombre": nombre, "contig_pos": contig_pos, "sequence": sequence})

# Create an empty dictionary to store the taxonomical information from the file "fungi_taxid_taxonomy"
tax_info_dict = {}

# Open the taxonomy file
with open(taxonomy_file, 'r') as file:
    for line in file:
        cols_taxonomy_file = line.strip().split('\t')  # The file "fungi_taxid_taxonomy.txt" has 2 columns separated by tab
        genome_code = '_'.join(cols_taxonomy_file[0].split('_')[:2])  # Extract genome code (first two elements)
        tax_info = cols_taxonomy_file[1]  # Extract taxonomic info
        tax_info_dict[genome_code] = tax_info  # Add genome code into the dictionary "tax_info_dict"

# Create another dictionary to group sequences by contig and positions
contig_dict = {}


for item in sequences:   # For each processed sequence...
    nombre = item["nombre"]
    contig_pos = item["contig_pos"]
    sequence = item["sequence"]
    nombre = re.sub(r"sp\.\s*'([^']+)'\b", r" \1", nombre) # Replace names like "Genus sp. 'specie'" by "Genus species"
    
    # if the contig already exist in the dictionary, combine the sequence (avoiding duplications)
    if contig_pos not in contig_dict:
        contig_dict[contig_pos] = {"sequences": [sequence], "tax_info": tax_info_dict.get(nombre, "Unknown")}
    else:
        contig_dict[contig_pos]["sequences"].append(sequence)

# Create the output file (a fasta file with associated taxonomy information)
with open(taxonomy_fasta, 'w') as outfile:
    for contig_pos, data in contig_dict.items():
        combined_sequence = ''.join(data["sequences"])# Combine all sequence of the same "contig_pos" into a single sequence              
        header = f">{contig_pos}\t{data['tax_info']}" # Change the header so it includes taxonomy info in the correct format
        header = re.sub(r"sp\.\s*'([^']+)'", r"\1", header)# Replace "Genus sp. 'specie' by "Genus species" (needed because in the taxonomy file there are some weird species names like, for example, Acidomyces sp. 'richmondensis')
        outfile.write(f"{header}\n{combined_sequence}\n") #Write the headear and the sequence

print("Fasta with taxonomic information saved at", taxonomy_fasta)



##############################################################################################################################################################################################################################
########################################################################### SEARCH FOR SHARED SEQUENCES AMONG DIFFERENT SPECIES ##############################################################################################
##############################################################################################################################################################################################################################

def process_fasta(taxonomy_fasta): # Define the function to process the fasta file with taxonomy data already added
    seqs = defaultdict(list)  # Create a dictinary
    seq_counter = defaultdict(int)  # Create a dictinary to count how many times each sequence appears in the fasta (number of repeats)
    
    with open(taxonomy_fasta, 'r') as f:
        header = None
        sequence = None
        species = set()  # Create a set to avoid duplicates
        contig = None  # Variable to store the contig
        
        for line in f:
            line = line.strip()
            
            if line.startswith(">"):  # If the line starts with ">" ...
                if header and sequence:  # If there is a header and a sequence...
                    seqs[sequence].append((contig, species)) # Add info
                    seq_counter[sequence] += 1  # Sum 1 to the counter of this sequence
                
                # Process the new line
                header = line 
                contig = re.findall(r">([A-Za-z0-9\.\-]+)", line) # Extract the ID with no position range (find ">" followed by any letter, number, dots or  hyphens"
                contig = '|'.join(contig)  # Combine several identifiers with "|"
                
                species = set(re.findall(r"S__([a-zA-Z0-9_ ]+)", line))  # Extract the species name
                sequence = ""  # Restart sequence
            else:
                sequence += line  # Continue adding sequence
        
        # Add the last sequence
        if header and sequence:
            seqs[sequence].append((contig, species))
            seq_counter[sequence] += 1  # Add 1 to the counter
    
    return seqs, seq_counter


def filter_seqs(seqs): # Function to filter the sequences associated with more than one species
    result = []
    for sequence, list_info in seqs.items():
        different_species = set() # Combine all the species found having the same sequence into a set 
        contigs = set()
        for contig, species in list_info:
            contigs.add(contig)
            different_species.update(species)
        
        if len(different_species) > 1: # If the lenght of "different_species" is higher than 1, i.e., if there is more than one species associated with the same sequence...
            result.append((sequence, contigs, different_species)) # Add sequence, contigs and species 
            print("The species", different_species, "share the same sequence\n")
    return result


def save_result(result, output_table, seq_counter): #Save the sequences shared by several species 
    with open(output_table, 'w') as f:
        f.write("Species\tContig(s)\tShared sequence\tNo. of repeats\n") # Write the header of the file (terms separated by tab)
         
        for sequence, contigs, species in result:
            ordered_contigs = sorted(contigs)
            repeats = seq_counter[sequence]
            f.write(f"{'|'.join(species)}\t{'|'.join(ordered_contigs)}\t{sequence}\t{repeats}\n") # Write the data (species, contigs, sequence adn number of repeats in the fasta)
            

# Call the functions
seqs, seq_counter = process_fasta(taxonomy_fasta) # Process the taxonomy_fasta file
result = filter_seqs(seqs) # Result of the funtion that filters those seqs with more than one species associated
save_result(result, output_table, seq_counter) # Save the result as a table file
print(f"Shared sequences by different species have been saved in: {output_table}")



###############################################################################################################################################################################################################################
######################################################### GETTING UNIQUE SEQUENCES FROM THE FASTA FILE CONTAINING ALL SEQUENCES ALREADY ANNOTATED WITH TAXONOMY DATA  #########################################################
###############################################################################################################################################################################################################################

unique_sequences = {}  # Create a dictionary to store all unique sequences

with open(taxonomy_fasta, 'r') as infile:  # Read the fasta with all sequences
    header = None
    sequence = []

    for line in infile:
        line = line.strip()  # Delete unnecesary breaks

        if line.startswith(">"): # if my line starts with ">"...
            if header and sequence: # If there is a sequence, save it
                seq_str = ''.join(sequence)  # Save the seq as a string
                # Extract ID of the genome and taxonomic info
                header_info = header.split("\t", 1)  # Divide the info in the header (contig and positions are separated from the taxonomical info by one tab)
                genome_code = header_info[0]  # Get contig and positions
                taxoinfo = header_info[1] if len(header_info) > 1 else ""  # If header_info has more than one element, get only the second one (header_info[1])

                if seq_str not in unique_sequences:  # If the sequence does not exit yet, add it
                    unique_sequences[seq_str] = [taxoinfo, genome_code, seq_str]

            # Process the new header
            header = line.strip()
            sequence = []  # Restart "sequence" for the next header
        else:
            # Add the seq
            sequence.append(line.strip())

    # Process the last seq after the for loop
    if header and sequence:
        seq_str = ''.join(sequence)
        header_info = header.split("\t", 1)  # Divide the header into contig info and taxonomy info
        genome_code = header_info[0]  # Get the contig or genome code (first part of of the header)
        taxoinfo = header_info[1] if len(header_info) > 1 else ""  # Taxonomy info
        if seq_str not in unique_sequences:  # If the sequence does not exit, add it
            unique_sequences[seq_str] = [taxoinfo, genome_code, seq_str]

# Write the output file with unique sequences
with open(Uniqueseqs_taxonomy, 'w') as outfile:
    for taxoinfo, genome_code, seq in unique_sequences.values():
        outfile.write(f"{genome_code}\t{taxoinfo}\n")
        outfile.write(f"{seq}\n")
print(f"Unique sequences have been saved in {Uniqueseqs_taxonomy}")


###############################################################################################################################################################################################################################
######################################################### LAST FILTERING STEP: DOUBLE CHECKING THAT THERE ARE NO IDENTICAL SEQS AND FINDING AND DELETING SUBSEQUENCES  ########################################################
###############################################################################################################################################################################################################################


def process_fasta(Uniqueseqs_taxonomy): 
    seqs = defaultdict(list)  # Create a dictionary to store sequences
    seq_counter = defaultdict(int)  # Create a dictionary to count how many times each sequence appears
    headers = []  # List to store the order of sequence headers
    
    with open(Uniqueseqs_taxonomy, 'r') as f:
        header = None
        sequence = None
        species = set()  # Set to store species (to avoid duplicates)
        contig = None  # Variable to store the contig identifier (singular, since each sequence is related to one contig)
        
        for line in f:
            line = line.strip()  # Remove leading/trailing spaces
            
            if line.startswith(">"):  # If the line starts with ">", it's a header line
                if header and sequence:  # If we have a previous header and sequence...
                    seqs[sequence].append((contig, species, header))  # Store the sequence information along with the header
                    seq_counter[sequence] += 1  # Increment the counter for this sequence
                
                header = line  # Store the current header
                contig = re.findall(r">(\S+)", line)  # Extract the full contig ID
                contig = contig[0] if contig else None  # Ensure contig is a single string, not a list
                
                species = set(re.findall(r"S__([a-zA-Z0-9_ ]+)", line))  # Extract species names
                sequence = ""  # Reset the sequence for the new record
                headers.append(header)  # Keep track of the header order
            else:
                sequence += line  # Add the current line to the sequence
        
        # After the loop, save the last sequence
        if header and sequence:
            seqs[sequence].append((contig, species, header))
            seq_counter[sequence] += 1  # Increment the counter for this sequence
    
    return seqs, seq_counter, headers


def lastfilter_seqs(seqs): 
    last_result = []  # List to store the filtered results
    identical_sequences_found = False  # Flag to check if we find identical sequences
    subsequences = defaultdict(list)  # Dictionary to store subsequences (smaller sequences within larger ones)
    
    for sequence, list_info in seqs.items():
        different_species = set()  # Set to store all species found for the current sequence
        contigs = set()  # Set to store contig identifiers
        
        # Go through each contig and species associated with the current sequence
        for contig, species, header in list_info:
            contigs.add(contig)
            different_species.update(species)  # Add species to the set
        
        if len(different_species) > 1:  # If more than one species is associated with the same sequence...
            last_result.append((sequence, contigs, different_species))  # Add to the result list
            identical_sequences_found = True  # We found identical sequences
    
    if not identical_sequences_found:  # If no identical sequences were found, look for subsequences
        print("No identical sequences found, searching for subsequences...\n")
        
        # Now check for subsequences (sequences that are part of a longer sequence)
        for sequence in seqs:
            for other_sequence in seqs:
                if other_sequence != sequence and sequence in other_sequence:  # If the sequence is part of another
                    subsequences[sequence].append(other_sequence)  # Add the subsequence
        
        # Include the subsequences in the result
        for subsequence, subseq_list in subsequences.items():
            for contig, species, header in seqs[subsequence]:
                last_result.append((subsequence, contig, species))  # Add subsequence to the result list
    
    # Remove subsequences from the original sequences
    all_sequences = set(seqs.keys())  # Get all the sequences in the original FASTA
    subsequences_found = set(subsequences.keys())  # Get the subsequences found
    
    remaining_sequences = all_sequences - subsequences_found  # Remove the subsequences from the original sequences
    #print(f"Remaining sequences after removing subsequences: {len(remaining_sequences)}\n")

    return last_result, remaining_sequences


def save_last_result(last_result, Subseq_table, seq_counter): 
    # Function to save the filtered results into a file
    with open(Subseq_table, 'w') as f:
        f.write("Species\tContig and positions\tShared sequence/Subsequence\tNo. of repeats\n")  # Write header to the output file
         
        # Write each entry in the result to the file
        for sequence, contigs, species in last_result:
            repeats = seq_counter[sequence]  # Get the number of repeats for the sequence
            f.write(f"{'|'.join(species)}\t{(contigs)}\t{sequence}\t{repeats}\n")  # Write data in the table format


def save_fasta(seqs, remaining_sequences, Nosubseq_fasta): 
    # Function to save the remaining sequences (without subsequences) into a new FASTA file
    with open(Nosubseq_fasta, 'w') as f:
        for sequence in remaining_sequences:
            for contig, species, header in seqs[sequence]:
                f.write(f"{header}\n")  # Write the original header
                f.write(f"{sequence}\n")


# Call the functions
seqs, seq_counter, headers = process_fasta(Uniqueseqs_taxonomy)  # Process the FASTA file with taxonomy data
last_result, remaining_sequences = lastfilter_seqs(seqs)  # Filter sequences that are associated with more than one species
save_last_result(last_result, Subseq_table, seq_counter)  # Save the result to an output table file
save_fasta(seqs, remaining_sequences, Nosubseq_fasta)  # Save the remaining sequences to a new FASTA file
print(f"Subsequences have been saved in: {Subseq_table}\n")
print(f"Filtered FASTA without subsequences has been saved in: {Nosubseq_fasta}\n")