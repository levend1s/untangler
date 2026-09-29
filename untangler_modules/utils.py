import os
import math
import gffpandas.gffpandas as gffpandas

def process_input_files(INPUT):
    input_files = {}

    with open(INPUT, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("#") or line.strip() == "":
                continue
            else:
                line = line.split()

                input_files[line[0]] = {
                    'group': line[1], 
                    'type': line[2], 
                    'path': line[3]
                }

    return input_files

def process_annotation_file(annotation_file_path=None):
    print("LOG - loading gff file...")
    # load annotation file and find indexes for all parent children
    ANNOTATION_FILE = gffpandas.read_gff3(annotation_file_path)

    gff_df = ANNOTATION_FILE.attributes_to_columns()
    gff_df['strand'] = gff_df['strand'].astype('category')
    gff_df['seq_id'] = gff_df['seq_id'].astype('category')
    gff_df['ID'] = gff_df['ID'].astype('category')
    gff_df['type'] = gff_df['type'].astype('category')
    gff_df['start'] = gff_df['start'].astype('int')
    gff_df['end'] = gff_df['end'].astype('int')

    return gff_df

def process_genome_file(genome_file_path=None, contig_lengths = None):
    print("LOG - loading fasta file...")

    FASTA_DICT = {}
    FASTA_LINE_LENGTH = 60

    # read two lines and set the line length as the second line
    if os.path.isfile(genome_file_path):
        with open(genome_file_path) as f:
            line = f.readline()
            line = f.readline()
            FASTA_LINE_LENGTH = len(line.strip())
            print("LOG - fasta line len = {}...".format(FASTA_LINE_LENGTH))
            print(contig_lengths)



    if os.path.isfile(genome_file_path):
        with open(genome_file_path) as f:
            line = f.readline()
            current_contig = None

            while line:
                # is line a header
                # PLASMODIUM SPECIFIC
                if line.startswith('>'):
                    if not contig_lengths:
                        header_attrs = line.split('|')
                        for h in header_attrs:
                            h = h.strip()
                            if h.startswith('>'):
                                contig = h[1:]
                                FASTA_DICT[contig] = {}
                                current_contig = contig
                            else:
                                k, v = h.split('=')
                                FASTA_DICT[current_contig][k] = v
                        FASTA_DICT[current_contig]['length'] = int(FASTA_DICT[current_contig]['length'])
                    else:
                        header_attrs = line.split(' ')
                        current_contig = header_attrs[0][1:]
                        FASTA_DICT[current_contig] = {}
                        FASTA_DICT[current_contig]['length'] = contig_lengths[current_contig]
                    
                    num_lines_seq_fasta = math.ceil(FASTA_DICT[current_contig]['length'] / FASTA_LINE_LENGTH) 
                    FASTA_DICT[current_contig]['sequence'] = [None] * num_lines_seq_fasta
                    i = 0
                else:
                    if line.startswith('>'):
                        current_contig = None

                    elif current_contig:
                        # FASTA_DICT[current_contig]['sequence'] += line.strip()
                        FASTA_DICT[current_contig]['sequence'][i] = line.strip()
                        i += 1

                line = f.readline()
    
    # verify each chromosome is the correct length
    for k, v, in FASTA_DICT.items():
        FASTA_DICT[k]['sequence'] = "".join(FASTA_DICT[k]['sequence'])
        print("{} - loaded seq: {}, header length: {}".format(k, len(v['sequence']), v['length']))

    return FASTA_DICT


