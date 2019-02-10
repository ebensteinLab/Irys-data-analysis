import os
import re

min_corr_lim = -2
root_dir = r'D:\misc\r2om2_from_scratch\wavelets\smooth_wgbs\automated_analysis'
subdir_name_regex = re.compile(r'(chr.*)_(\d+)_(\d+)')
output_bed_file_path = r'D:\misc\r2om2_from_scratch\wavelets\smooth_wgbs\automated_analysis\corr_vals_all_round.txt'
corr_file_name = "win.size.10000.corr.vals.smooth.txt"

with open(output_bed_file_path, 'w') as output_file:
    sub_dirs = [name for name in os.listdir(root_dir) if os.path.isdir(name)]
    dir_count = 0
    for dir_name in sub_dirs:
        print(dir_name)
        dir_count += 1
        if dir_count % 100 == 0:
            print(dir_count)
        name_regex_res = re.split(subdir_name_regex, dir_name)
        chrom = name_regex_res[1]
        start = int(name_regex_res[2])
        input_file_path = os.path.join(root_dir,dir_name,corr_file_name)
        with open(input_file_path) as input_file:
            for line in input_file:
                corr_val = round(float(line.strip()),1)
                if corr_val >= min_corr_lim and corr_val != 0.0:
                    output_file.write("\t".join([chrom,str(start),str(start+1),str(corr_val)])+"\n")
                start += 1
    
