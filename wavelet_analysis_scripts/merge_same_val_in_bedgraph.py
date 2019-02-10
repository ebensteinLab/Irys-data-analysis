bedgraph_input_file = "corr_vals_all_round.chrY.sorted.txt"
bedgraph_output_file = "corr_vals_all_round.chrY.sorted.merged.txt"

with open(bedgraph_input_file) as input_file:
    with open(bedgraph_output_file,'w') as output_file:
        #output_file.write("track type=bedGraph\n")
        line = input_file.readline()
        info = line.strip().split()
        prev_chrom = info[0]
        prev_start = info[1]
        prev_end = info[2]
        prev_val = info[3]

        line = input_file.readline()
        while line:
            info = line.strip().split()
            if info[0] != prev_chrom:
                #write to file and update
                try:
                    lineToPrint = prev_chrom+"\t"+prev_start+"\t"+prev_end+"\t"+prev_val+"\n"
                    output_file.write(lineToPrint)
                    prev_chrom = info[0]
                    prev_start = info[1]
                    prev_end = info[2]
                    prev_val = info[3]
                except:
                    print(type(prev_chrom), type(prev_start), type(prev_val))
            
            else:
                if info[3] != prev_val or info[1] != prev_end:
                    try:
                        #write to file and update
                        lineToPrint = prev_chrom+"\t"+prev_start+"\t"+prev_end+"\t"+prev_val+"\n"
                        output_file.write(lineToPrint)
                        prev_chrom = info[0]
                        prev_start = info[1]
                        prev_end = info[2]
                        prev_val = info[3]
                #else do nothing, just read the next line
                    except:
                        print(line)
                else:
                    prev_end = info[2]

            line = input_file.readline()
        lineToPrint = prev_chrom+"\t"+prev_start+"\t"+info[1]+"\t"+prev_val+"\n"
        output_file.write(lineToPrint)
                    
            

print("analysis complete")
