# allowed to use biopython and re libraries
from Bio import SeqIO
import re

# the rest of the libraries are built in with python
import random
import statistics
from sys import argv, exit
import os



def read_promoters_file(promoters_file):
    """
    This function reads the list of promoter motifs from the given file and returns them as a list.
    """
    # create an empty list to store the motifs
    promoters = []
    # open the promoters txt file and store its content in the list
    with open(promoters_file, "r") as file:
        for line in file:
            lines = line.strip()
            promoters.append(lines)
    return promoters

def read_zea_mays_genes_file(zea_mays_genes_file):
    """
    This function reads the Zea mays gene IDs from the given file and returns them in a list.
    """
    # create an empty list to store the zea mays genes
    zea_mays_genes = []
    # open the zea mays gene txt file and store its content in the list
    with open(zea_mays_genes_file, "r") as file:
        for line in file:
            lines = line.strip()
            zea_mays_genes.append(lines)
    return zea_mays_genes

def read_gff_file(gff_file):
    """
    This function reads the GFF3 file to get the TSS for each gene.
    It returns a dictionary that has the location and strand information for each gene.
    """

    # open the zea mays whole genome file and store relevant information in a list
    Zea_mays_dict = {}
    with open(gff_file, "r") as file:
        for line in file:
            # skip the lines with #
            if line.startswith("#"):
                continue
            line_list = line.split("\t")
            # only store the information if it is in chromosomes 1 to 10 and is a gene
            if line_list[0] in map(str, range(1, 11)) and line_list[2] == "gene":
                gene_id = line_list[8].split(";")
                gene_id_string = gene_id[0].replace("ID=gene:", "")
                chrom = line_list[0]
                start = int(line_list[3])
                end = int(line_list[4])
                strand = line_list[6]

                # Some genes have multiple TSS, as mentioned in the hints section of the assignment
                # we need to only select the most upstream one
                
                # Check if the gene is already in Zea_mays_dict
                if gene_id_string in Zea_mays_dict:
                    # For positive strand, keep the smallest start position
                    if strand == "+" and start < int(Zea_mays_dict[gene_id]["Start Position"]):
                        Zea_mays_dict[gene_id] = {"Chromosome": chrom, "Start Position": start, "End Position": end, "Strand": strand}
                    # For negative strand, keep the largest start position, since it is a reverse complement strand
                    elif strand == "-" and start > int(Zea_mays_dict[gene_id]["Start Position"]):
                        Zea_mays_dict[gene_id] = {"Chromosome": chrom, "Start Position": start, "End Position": end, "Strand": strand}
                else:
                    # Add the gene to the dictionary for the first time
                    Zea_mays_dict[gene_id_string] = {"Chromosome": chrom, "Start Position": start, "End Position": end, "Strand": strand}
    return Zea_mays_dict

def load_fasta_files(fasta_dir):
    """
    This function loads the FASTA files for chromosomes 1 to 10.
    It returns a dictionary where keys are chromosome identifiers and values are SeqIO index objects.
    """
    # store all the fasta file sequences as values for each chromosome IDs as keys
    fasta_files = {}
    for i in range(1, 11):
        fasta_file = os.path.join(fasta_dir, f"Zea_mays.B73_RefGen_v4.dna.chromosome.{i}.fa")
        fasta_files[str(i)] = SeqIO.index(fasta_file, "fasta")
    return fasta_files

def get_promoter_sequences(Zea_mays_dict, fasta_files, promoter_length=500):
    """
    This function goes through each gene in the dictionary and gets the upstream 500bp promoter sequence.
    It modifies Zea_mays_dict by adding a 'Sequence' key for each gene.
    """

    # store the 500nt sequences for each gene. 
    # initially, I was going to do this only for the genes in the zea_mays_genes.txt file
    # but we need to obtain the sequences of randomly selected genes as well
    # so we might as well store the 500nt sequences of all genes so it's easier to retrieve it later on

    for key, value in Zea_mays_dict.items():
        chrom_id = value["Chromosome"]
        strand = value["Strand"]
        start = int(value["Start Position"])

        # read individual fasta file based on the chromosome position of the gene
        fasta_seq = fasta_files[chrom_id]
        chrom_seq = fasta_seq[chrom_id].seq

        # depending on positive or negative strand, store the 500nt sequence before or after the start position of the gene
        if strand == "+":
            startSeqPos = max(0, start - promoter_length) - 1
            endSeqPos = start - 1
            sequence = chrom_seq[startSeqPos:endSeqPos]
        elif strand == "-":
            startSeqNeg = start
            endSeqNeg = start + promoter_length
            # convert the sequence to the reverse compliment since it is in the negative strand
            sequence = chrom_seq[startSeqNeg:endSeqNeg].reverse_complement()
        
        # Check if sequence contains 'N', indicating unknown regions
        if "N" in sequence:
            # Find the position of the first N and shorten the sequence to exclude it and all sequence after it
            first_N = sequence.index("N")
            sequence = sequence[:first_N]
        
        # add the 500nt long sequence to the dictionary as a new key, value pair in the main dictionary, assigned to the key of the gene it is currently on
        Zea_mays_dict[key]["Sequence"] = str(sequence)

def convert_motifs_to_regex(promoters):
    """
    This function converts each motif into a regular expression pattern that can work with the re library.
    For example, [AG]CCGAC becomes (A|G)CCGAC so it can be  used with Python's re module.
    """

    # convert motifs into regex patterns for use in re library
    motif_patterns = {}
    for motif in promoters:
        changed_pattern = motif.replace("[", "(").replace("]", ")")
        compiled_pattern = re.compile(changed_pattern, re.IGNORECASE)
        # creates a dictionary with motifs as keys, and values as patterns that are compatible with the re library
        # e.g Key: [AG]CCGAC , Value: (A|G)CCGAC
        motif_patterns[motif] = compiled_pattern
    return motif_patterns

def count_motifs_in_genes(zea_mays_genes, Zea_mays_dict, motif_patterns):
    """
    This function counts how many times each motif show up in the promoter regions of the selected genes.
    It returns a dictionary that keeps track of the total count of each motif (motif_counts_from_file) and a per gene count of each motif (gene_motif_counts)
    """

    # set initial count of 0 to all the motifs from the file
    motif_counts_from_file = {}
    for motif in motif_patterns.keys():
        motif_counts_from_file[motif] = 0
        
    # Added to store counts for selected genes
    gene_motif_counts = {}

    # count the number of times the motif appears in the gene promoters from the text file provided
    for key in zea_mays_genes:
        if key in Zea_mays_dict:
            sequence = Zea_mays_dict[key]["Sequence"]
            gene_motif_counts[key] = {}
            for motif, pattern in motif_patterns.items():
                # use the re library pattern to find matches
                matches = pattern.findall(sequence)
                motif_counts_from_file[motif] += len(matches)
                gene_motif_counts[key][motif] = len(matches)
    return motif_counts_from_file, gene_motif_counts

def count_motifs_in_random_genes(Zea_mays_dict, motif_patterns, num_genes_to_sample, iterations=5):
    """
    This function picks random sets of genes (equal in number to the selected gene list) multiple times (5 as per the assignment instructions).
    It counts how many times the motif shows up in these random sets, so that we can compare how the selected genes differ from random ones.
    It returns the total count of each motif across all genes in each iteration of the run (random_sample_counts) and a per gene count of each motif (random_gene_motif_counts)
    """

    # do the same analysis for randomly selected genes 5 times
    # create an empty list to store total counts
    all_genes = list(Zea_mays_dict.keys())
    random_sample_counts = []
    random_gene_motif_counts = []
    
    #run the analysis 5 times
    for i in range(iterations):
        random_genes = random.sample(all_genes, num_genes_to_sample)
        # set initial motif count to 0
        motif_counts_random = {}
        for motif in motif_patterns.keys():
            motif_counts_random[motif] = 0
            
        run_gene_motif_counts = {}
        
        # calculate the number of times motifs show up in the sequence
        for gene in random_genes:
            sequence = Zea_mays_dict[gene]["Sequence"]
            run_gene_motif_counts[gene] = {}
            for motif, pattern in motif_patterns.items():
                # use the re library pattern to find matches
                matches = pattern.findall(sequence)
                motif_counts_random[motif] += len(matches)
                run_gene_motif_counts[gene][motif] = len(matches)
        random_sample_counts.append(motif_counts_random)
        random_gene_motif_counts.append(run_gene_motif_counts)

    return random_sample_counts, random_gene_motif_counts

def write_motif_counts_to_file(promoters, motif_counts_from_file, random_sample_counts, output_file="motif_counts.txt"):
    """
    This function writes motif counts to a file with the following columns:
    1) Motif
    2) Total counts of that motif across all selected genes combined
    3-7) The total counts of that motif for each of the random sampling runs
    """

    # store the counts in a newly created file
    with open(output_file, "w") as out:
        # We know there are 5 random runs, so we'll have 7 columns total
        header = ["Motif", "Selected", "Random Sample Run 1", "Random Sample Run 2", "Random Sample Run 3", "Random Sample Run 4", "Random Sample Run 5"]
        out.write("\t".join(header) + "\n")

        for motif in promoters:
            # Get the total selected count for this motif
            selected_total = motif_counts_from_file[motif]
            # Get the counts for each random run
            random_counts = []
            for run in random_sample_counts:
                random_counts.append(run[motif])
                
            # Combine all counts into one line
            counts = [motif, selected_total] + random_counts
            counts_str = "\t".join(map(str, counts))
            out.write(f"{counts_str}\n")


def calculate_t_test(selected_counts, random_counts):
    """
    This function helps calculate a basic t-test for comparing two sets of motif counts:
    one from the selected genes and one from the random sets.
    This uses the statistics library which is included with python. 
    This function would not be needed if we could use libraries such as scipy.
    Library manual used: https://docs.python.org/3/library/statistics.html
    """

    # Calculate t-test using statistics library
    mean_selected = statistics.mean(selected_counts)
    mean_random = statistics.mean(random_counts)
    var_selected = statistics.variance(selected_counts)
    var_random = statistics.variance(random_counts)
    t_stat = (mean_selected - mean_random) / (((var_selected / len(selected_counts)) + (var_random / len(random_counts))) ** 0.5)
    return t_stat

def get_over_under_represented_motifs(promoters, zea_mays_genes, gene_motif_counts, random_gene_motif_counts, motif_counts_from_file, output_file="over_under_represented_motifs.txt", threshold=2):
    """
    This function checks which motifs are significantly more or less frequent in the selected genes compared to random sets.
    It uses a simple t-test and writes over- or under-represented motifs (based on a threshold) to a file.
    """

    # Write over- or under-represented motifs in a separate file
    num_genes_to_sample = len(zea_mays_genes)
    with open(output_file, "w") as out:
        out.write("Motif\tAverage of selected genes in file provided\tAverage of Random sample runs\tT-Stat\tStatus\n")
        selected_avg_counts = {}
        for motif in promoters:
            selected_avg_counts[motif] = motif_counts_from_file[motif] / num_genes_to_sample
        
        # skip averages with 0, since there is nothing to compare
        for motif in promoters:
            selected_avg = selected_avg_counts[motif]
            if selected_avg == 0:
                continue
            
            # get the counts for the selected genes
            selected_counts = []
            for key in zea_mays_genes:
                selected_counts.append(gene_motif_counts[key][motif])
            
            # get the counts for the random runs
            random_counts = []
            for run_counts in random_gene_motif_counts:
                for counts in run_counts.values():
                    random_counts.append(counts[motif])
                
            # perform t-test using the function created earlier
            t_stat = calculate_t_test(selected_counts, random_counts)
            mean_random = statistics.mean(random_counts)

            # using a threshold of 2, determine if motif is over or under represented
            if t_stat > 2:
                status = "Over-represented"
            elif t_stat < -2:
                status = "Under-represented"
            else:
                continue

            out.write(f"{motif}\t{selected_avg:.2f}\t{mean_random:.2f}\t{t_stat:.2f}\t{status}\n")


if __name__ == "__main__":
    
    # As per the hints and warnings section in the assignment instructions, I have not hard-coded any files or directories.
    # So we need to get user input for the files
    
    # Check if the correct number of arguments is provided
    if len(argv) < 5:
        print("Usage: python3 script.py <promoters_file, e.g. promoters.txt> <genes_file, e.g zea_mays_genes.txt> <gff_file> <fasta_dir (name of folder containing the 1-10 fasta files>")
        exit()
        
    # Below are the steps done to complete the assignment. 
    # Each step relates to the steps indicated in the assignment as best as it can, so its hopefully easy for reader to follow along.
    # Each of the steps below uses the functions made above.    
    
    # Step 1: Get the file names from unix arguments
    promoters_file = argv[1]
    zea_mays_genes_file = argv[2]
    gff_file = argv[3]
    fasta_dir = argv[4]
    

    # Step 2. Read GFF to get TSS
    Zea_mays_dict = read_gff_file(gff_file)

    # Step 2 continued: Load gene IDs
    zea_mays_genes = read_zea_mays_genes_file(zea_mays_genes_file)

    # Step 3: Load FASTA files for each chromosome
    fasta_files = load_fasta_files(fasta_dir)

    # Step 3 continued: Get the 500nt upstream sequences
    get_promoter_sequences(Zea_mays_dict, fasta_files, promoter_length=500)

    # Step 4: Read promoter motifs file and convert motifs to regex to be used with re library
    promoters = read_promoters_file(promoters_file)
    motif_patterns = convert_motifs_to_regex(promoters)

    # Step 4 continued: Count motifs in the selected genes
    motif_counts_from_file, gene_motif_counts = count_motifs_in_genes(zea_mays_genes, Zea_mays_dict, motif_patterns)

    # Step 5: Do the random sampling and counting
    random_sample_counts, random_gene_motif_counts = count_motifs_in_random_genes(Zea_mays_dict, motif_patterns, len(zea_mays_genes), iterations=5)

    # Step 5 continued: Write motif counts to file
    write_motif_counts_to_file(promoters, motif_counts_from_file, random_sample_counts, output_file="motif_counts.txt")

    # Step 6: Find the over or under represented motifs and write to a separate file
    get_over_under_represented_motifs(promoters, zea_mays_genes, gene_motif_counts, random_gene_motif_counts, motif_counts_from_file, output_file="over_under_represented_motifs.txt", threshold=2)
