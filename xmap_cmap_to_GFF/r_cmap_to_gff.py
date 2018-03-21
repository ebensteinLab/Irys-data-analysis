chrom_lengths = {'chr1':249250621, 'chr2':243199373, 'chr3':198022430, 'chr4':191154276,
                 'chr5':180915260, 'chr6':171115067, 'chr7':159138663, 'chr8':146364022,
                 'chr9':141213431, 'chr10':135534747, 'chr11':135006516, 'chr12':133851895,
                 'chr13':115169878, 'chr14':107349540, 'chr15':102531392, 'chr16':90354753,
                 'chr17':81195210, 'chr18':78077248, 'chr19':59128983, 'chr20':63025520,
                 'chr21':48129895, 'chr22':51304566, 'chrX':155270560, 'chrY':59373566}

label_model = "CDS"
r_cmap_input_path = "MoleculeQualityReport_r.cmap"
gff_output_path = "bspq.reference.gff3"

prev_chrom = ""
with open(r_cmap_input_path) as input_file:
    with open(gff_output_path, 'w') as output_file:
        output_file.write("##gff-version 3\n")
        for line in input_file:
            if line[0] != "#":
                info = line.split()
                if info[0] == "23":
                    chrom = "chrX"
                elif info[0] == "24":
                    chrom = "chrY"
                else:
                    chrom = "chr"+info[0]
                pos = info[5]

                if chrom != prev_chrom:
                    chrom_len = chrom_lengths[chrom]+1
                    seq_region_line = "##sequence-region "+chrom+" 1 " + str(chrom_len)+"\n"
                    output_file.write(seq_region_line)
                    chrom_ID = "chromosome"+info[0]
                    chrom_line = "\t".join([chrom, ".", "transcript", "1", str(chrom_len), ".", ".", ".", "ID="+chrom_ID])
                    output_file.write(chrom_line+"\n")
                    prev_chrom = chrom


                posStart = float(pos) - 10
                posEnd = float(pos) + 10
                pos_line = "\t".join([chrom, ".", label_model, str(int(posStart)), str(int(posEnd)), ".", ".", ".", "Parent="+chrom_ID])
                output_file.write(pos_line+"\n")

print("analysis complete")
                
