chrom_lengths = {'chr1':249250621, 'chr2':243199373, 'chr3':198022430, 'chr4':191154276,
                 'chr5':180915260, 'chr6':171115067, 'chr7':159138663, 'chr8':146364022,
                 'chr9':141213431, 'chr10':135534747, 'chr11':135006516, 'chr12':133851895,
                 'chr13':115169878, 'chr14':107349540, 'chr15':102531392, 'chr16':90354753,
                 'chr17':81195210, 'chr18':78077248, 'chr19':59128983, 'chr20':63025520,
                 'chr21':48129895, 'chr22':51304566, 'chrX':155270560, 'chrY':59373566}

alignment_channel_model = "UTR"
query_channel_model = "CDS"
molecule_list_input_path = "C:\\cygwin64\\home\\Hila\\bin\\mol_gene_exp\\merged.conf.12.min.alignment.60.positions.list.intersect.gene.locations.3.sorted.txt"
gff_output_path = "C:\\cygwin64\\home\\Hila\\bin\\mol_gene_exp\\merged.conf.12.min.alignment.60.positions.list.intersect.gene.locations.3.gff3"

prev_chrom = ""

with open(molecule_list_input_path) as input_file:
    with open(gff_output_path, 'w') as output_file:
        output_file.write("##gff-version 3\n")
        for line in input_file:
            #info: molID, refID, molStart, molEnd, alignmentChannelPosList, queryChannelPosList
            info = line.split()
            if info[1] == "23":
                ref = "chrX"
            elif info[1] == "24":
                ref = "chrY"
            else:
                ref = "chr"+info[1]

            if ref != prev_chrom:
                chrom_len = chrom_lengths[ref]+1
                seq_region_line = "##sequence-region "+ref+" 1 " + str(chrom_len)+"\n"
                output_file.write(seq_region_line)
                prev_chrom = ref

            molID = "molecule"+info[0]
            mol_line = "\t".join([ref, ".", "transcript", info[2], info[3], ".", ".", ".", "ID="+molID])
            output_file.write(mol_line+"\n")

            alignment_positions = info[4].split(",")
            for alignment_pos in alignment_positions:
                posStart = int(alignment_pos) - 10
                posEnd = int(alignment_pos) + 10
                pos_line = "\t".join([ref, ".", alignment_channel_model, str(posStart), str(posEnd), ".", ".", "0", "Parent="+molID])
                output_file.write(pos_line+"\n")

            if len(info) > 5:
                query_positions = info[5].split(",")    
                for query_pos in query_positions:
                    posStart = int(query_pos) - 10
                    posEnd = int(query_pos) + 10
                    pos_line = "\t".join([ref, ".", query_channel_model, str(posStart), str(posEnd), ".", ".", ".", "Parent="+molID])
                    output_file.write(pos_line+"\n")

print("analysis complete")
                
                
