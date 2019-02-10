import os
import subprocess
import numpy as np
from scipy.stats.stats import pearsonr


def call_bedtools_intersect(file_a_path, file_b_path, output_path):
    os.chdir(r'C:\cygwin64\bin')
    cmd = ["bash", "-c", '/home/Hila/bedtools.exe intersect -a ' +file_a_path+ ' -b ' +file_b_path]
    with open(output_path, 'w') as output_file:
        subprocess.call(cmd, stdout=output_file)

def vector_from_bg(bg_input_path, start, end, output_file_path):
    length = end - start
    #print(start, end, length)
    output = np.zeros(length)
    with open(bg_input_path) as bg_file:
        for line in bg_file:
            info = line.split()
            start_pos = int(info[1]) - start
            end_pos = int(info[2]) - start
            val = float(info[3])
            for pos in range(start_pos, end_pos):
                #print(pos)
                output[pos] = val
    np.savetxt(output_file_path, output, fmt="%1.5f", delimiter="\n")


def smooth_to_smooth_vect(smooth_file_input_path, output_file_path):
    coefficients = []
    with open(smooth_file_input_path) as input_file:
        for line in input_file:
            val = float(line.strip())
            coefficients.append(val)
    rev_coefficients = coefficients[::-1]
    with open(output_file_path, 'w') as output_file:
        for pos in range(len(coefficients)//2):
            av_val = (coefficients[pos] + rev_coefficients[pos]) / 2
            output_file.write(str(av_val)+"\n")

def calc_correlation_values(wgbs_cwt_file, irys_cwt_file, window_size, output_file_path):
    irys_cwt_matrix = np.loadtxt(irys_cwt_file, delimiter="\t")
    wgbs_cwt_matrix = np.loadtxt(wgbs_cwt_file, delimiter="\t")
    print("done reading files")
    coeff_values = np.zeros(irys_cwt_matrix.shape[0]-2*window_size+1)
    #coeff_matrix = np.zeros((irys_cwt_matrix.shape[1], irys_cwt_matrix.shape[0]-2*window_size+1))
    coeff_array_curr_pos = 0
    curr_pos = window_size
    while curr_pos <= irys_cwt_matrix.shape[0] - window_size:
        if curr_pos % 10000 == 0:
            print(curr_pos)
        irys_curr_values = irys_cwt_matrix[curr_pos-window_size:curr_pos+window_size]
        wgbs_curr_values = wgbs_cwt_matrix[curr_pos-window_size:curr_pos+window_size]
        corr_coeff = pearsonr(irys_curr_values, wgbs_curr_values)[0]
        coeff_values[coeff_array_curr_pos] = corr_coeff
        curr_pos += 1
        coeff_array_curr_pos += 1
    np.savetxt(output_file_path, coeff_values, fmt="%0.10f", delimiter="\n") 
    

regions_input_file = r'D:\misc\r2om2_from_scratch\wavelets\smooth_wgbs\automated_analysis\regions.txt'
input_output_dir = r'D:\misc\r2om2_from_scratch\wavelets\smooth_wgbs\automated_analysis'
irys_full_data_file_path = '/cygdrive/d/misc/r2om2_from_scratch/wavelets/smooth_wgbs/automated_analysis/min.conf.11.min.alignment.60.labels.channel.1.merged.plus.minus.500bp.min.5.mols.labels.to.mols.ratio.bedgraph'
wgbs_full_data_file_path = '/cygdrive/d/misc/r2om2_from_scratch/wavelets/smooth_wgbs/automated_analysis/wgbs_pbmc_invert_sorted.filtered.blacklist.bedgraph'

total_regions = 0
with open(regions_input_file) as regions_file:
    for line in regions_file:
        total_regions += 1
print("total regions: ", total_regions)

curr_region = 0
with open(regions_input_file) as regions_file:
    for line in regions_file:
        curr_region += 1
        print("region: "+ str(curr_region) +" out of: " + str(total_regions))
        info = line.strip().split()
        region_name = "_".join(info)
        print(region_name)
        new_dir_path = os.path.join(input_output_dir, region_name)
        if not os.path.exists(new_dir_path):
            os.makedirs(new_dir_path)
        ROI_file_path = os.path.join(new_dir_path, 'ROI.txt')
        with open(ROI_file_path, 'w') as ROI_file:
            ROI_file.write("\t".join(info)+"\n")
        irys_bg_intersection_output_path = os.path.join(new_dir_path, 'irys.mol.ratio.'+".".join(info)+".bedg")
        wgbs_bg_intersection_output_path = os.path.join(new_dir_path, 'wgbs.invert.'+".".join(info)+".bedg")
        ROI_intersection_path = os.path.join('/cygdrive/d/misc/r2om2_from_scratch/wavelets/smooth_wgbs/automated_analysis/', region_name+"/", 'ROI.txt')
        print("running BEDTools intersect on Irys data")
        call_bedtools_intersect(irys_full_data_file_path, ROI_intersection_path, irys_bg_intersection_output_path)
        print("running BEDTools intersect on raw wgbs data")
        call_bedtools_intersect(wgbs_full_data_file_path, ROI_intersection_path, wgbs_bg_intersection_output_path)
        os.chdir(input_output_dir)
        wgbs_raw_bg_to_vector_output = os.path.join(new_dir_path, 'wgbs.invert.'+".".join(info)+".vector.txt")
        print("converting wgbs raw bedgraph to vector")
        vector_from_bg(wgbs_bg_intersection_output_path, int(info[1]), int(info[2]), wgbs_raw_bg_to_vector_output)
        wavelet_smooth_output = os.path.join(new_dir_path, 'wgbs.invert.'+".".join(info)+".smooth.n4.txt")
        print("calculating MODWT smooth")
        subprocess.check_call(['Rscript', 'MODWT_calc.R', wgbs_raw_bg_to_vector_output,wavelet_smooth_output],shell=False)
        smooth_vect_path = os.path.join(new_dir_path, 'wgbs.invert.'+".".join(info)+".smooth.n4.vect.txt")
        print("calculating smooth vector from smooth coefficients")
        smooth_to_smooth_vect(wavelet_smooth_output, smooth_vect_path)
        irys_bg_to_vector_output = os.path.join(new_dir_path, 'irys.mol.ratio.'+".".join(info)+".vector.txt")
        print("converting irys bedgraph to vector")
        vector_from_bg(irys_bg_intersection_output_path, int(info[1]), int(info[2]), irys_bg_to_vector_output)
        wgbs_smooth_cwt_abs_path = os.path.join(new_dir_path, 'wgbs.invert.'+".".join(info)+".DOG.abs.txt")
        irys_smooth_cwt_abs_path = os.path.join(new_dir_path, 'irys.mol.ratio.'+".".join(info)+".DOG.abs.txt")
        print("calculating CWT DOG wavelet decomposition for both data sets")
        subprocess.check_call(['Rscript', 'CWT_DOG_calc.R', smooth_vect_path, irys_bg_to_vector_output,
                               wgbs_smooth_cwt_abs_path, irys_smooth_cwt_abs_path],shell=False)
        window_size = 10000
        correlations_output_path = os.path.join(new_dir_path, 'win.size.'+str(window_size)+".corr.vals.smooth.txt")
        print("calculating correlation values")
        calc_correlation_values(wgbs_smooth_cwt_abs_path, irys_smooth_cwt_abs_path, window_size,
                                correlations_output_path)
        
        
            
        

