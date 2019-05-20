import sys
import argparse
import os
import numpy as np
from scipy import interpolate

parser = argparse.ArgumentParser(description='Convert BioNano xmap and cmap to labels and molecules BED files')
parser.add_argument('xmap', help='input xmap file')
parser.add_argument('q_cmap', help='input q_cmap file')
parser.add_argument('r_cmap', help='input r_cmap file')
parser.add_argument('key', help='ref cmap key file')
parser.add_argument('-o', '--output', dest='output_dir', metavar='OUTPUT DIR', help='output files directory. default current working directory')
parser.add_argument('-c', '--conf', dest='min_confidence', type=int, default=12, help='min alignment confidence. default 12')
parser.add_argument('-a', '--alignment', dest='alignment_label_channel', metavar='ALIGNMENT LABEL CHANNEL', type=int, default=1, help='alignment label channel. default 1')
parser.add_argument('-q', '--query', dest='query_label_channel', metavar='QUERY LABEL CHANNEL', type=int, default=2, help='query label channel. default 2')
parser.add_argument('-lb', '--label_bed', default=False, action='store_true', help='save query labels BED file. default False')
parser.add_argument('-mb', '--molecule_bed', default=False, action='store_true', help='save molecules containing query labels. default False')
parser.add_argument('-p', '--prefix', default="", help='output files prefix')

args = parser.parse_args()

xmap_infos = {}
alignment_q_cmap_positions = {} #dict of positions where key=ID and val=list of tuples (posID,pos) for positions corresponding to alignment label channel
query_q_cmap_positions = {} #dict of positions where key=ID and val=list of tuples (posID,pos) for positions corresponding to query label channel
r_cmap_positions = {} #dict of positions where key=chrom and val=list of tuples (posID,pos)
correct_positions = {} #dict of correct genomic positions where key=chrom and val = positions list
molRegions = {} #dict of stretched mol regions. key=ID, val=tup(start,end)

def parseKeyFile(key_file_path):
    '''parse the chromosome key file created when digesting a fasta reference file into a cmap''' 
    chrom_lengths = {}
    chrom_names = {}
    with open(key_file_path) as key_file:
        for line in key_file:
            if not line.startswith("#") and not line.startswith("C"):
                chr_id, chr_name, chr_len=line.strip().split()
                chrom_names[chr_id] = chr_name
                chrom_lengths[chr_name] = int(chr_len)
    return (chrom_lengths, chrom_names)

def parseXmap (xmap_input_path, min_confidence):
    '''read alignment information for all molecules that aligned to the reference'''
    with open(xmap_input_path) as input_file:
        for line in input_file:
            if line[0] != "#":
                line = line.strip()
                info = line.split()
                if float(info[8]) >= min_confidence:
                    ID = info[1]
                    xmap_infos[ID] = info

def findRCmapPositions (cmap_input_path, label_channel, posDict):
    '''read reference label positions from r_cmap'''
    with open(cmap_input_path) as input_file:
        for line in input_file:
            if line[0] != "#":
                line = line.strip()
                info = line.split()
                if int(info[4]) == label_channel:
                    ID = info[0]
                    pos = info[5]
                    if ID in posDict:
                        posDict[ID].append(float(pos))
                    else:
                        posDict[ID] = [float(pos)]


def skipHeader(file_handle):
    line = file_handle.readline()
    while line.startswith("#"):
        line = file_handle.readline()
    return line

def addLabel(labels_dict, label_info):
    label_pos = float(label_info[5])
    label_channel = int(label_info[4])
    if label_channel != 0:
        labels_dict[label_channel].append(label_pos)
    

def getNextMoleculeLabels(file, line):
    labels_dict = {1:[], 2:[]}
    mol_id = line.split()[0]
    while line and line.split()[0] == mol_id:
        label_info = line.strip().split()
        addLabel(labels_dict, label_info)
        line = file.readline()
    return (labels_dict,line)

def parseAlignmentString(alignment_string):
    matches = alignment_string[1:len(alignment_string)-1].split(")(")
    tup_matches = []
    for match in matches:
        splitPos = match.find(",")
        tup = (match[:splitPos],match[splitPos+1:])
        tup_matches.append(tup)
    return tup_matches

def getRefAndQueryPositions(labels_dict, alignment_label_channel, alignment_string_matches, mol_orientation, ref_id):
    ref_positions = []
    mol_positions = []
    for tup in alignment_string_matches:
        ref_pos = r_cmap_positions[ref_id][int(tup[0])-1]
        mol_pos = labels_dict[alignment_label_channel][int(tup[1])-1]
        if mol_pos not in mol_positions:
            ref_positions.append(ref_pos)
            mol_positions.append(mol_pos)

    if mol_orientation == "-":
        ref_positions = ref_positions[::-1]
        mol_positions = mol_positions[::-1]

    return (mol_positions, ref_positions)

def getStretchedLabelPositions(ius, labels_dict, query_label_channel):
    corrected_positions = []
    for label in labels_dict[query_label_channel]:
        corrected_positions.append(ius(label))
    return corrected_positions

def getMolRegion(corrected_query_labels, corrected_alignment_labels):
    if len(corrected_query_labels) == 0:
        mol_start_pos = min(corrected_alignment_labels)
        mol_end_pos = max(corrected_alignment_labels)
    else:
        mol_start_pos = min(min(corrected_query_labels), min(corrected_alignment_labels))
        mol_end_pos = max(max(corrected_query_labels), max(corrected_alignment_labels))
    return (mol_start_pos, mol_end_pos)

def writeMolLabelsToFile(mol_label_list, ref_id, padding, output_file_handle):
    chrom_name = chrom_names[ref_id]
    for label in mol_label_list:
        start_pos = int(label) - padding
        end_pos = int(label) + padding
        if start_pos < 0:
            start_pos = 0
        if end_pos > chrom_lengths[chrom_name]:
            end_pos = chrom_lengths[chrom_name]
        line = "\t".join([chrom_name, str(start_pos), str(end_pos)]) + "\n"
        output_file_handle.write(line)

def writeMolRegionToFile(corrected_query_labels, corrected_alignment_labels, mol_id, ref_id, padding, output_file_handle):
    mol_start_pos, mol_end_pos = getMolRegion(corrected_query_labels, corrected_alignment_labels)
    chrom_name = chrom_names[ref_id]
    mol_start_pos = int(mol_start_pos) - padding
    mol_end_pos = int(mol_end_pos) + padding
    if mol_start_pos < 0:
        mol_start_pos = 0
    if mol_end_pos > chrom_lengths[chrom_name]:
        mol_end_pos = chrom_lengths[chrom_name]
    line = "\t".join([chrom_name, str(mol_start_pos), str(mol_end_pos), mol_id, "0", "+"]) + "\n"
    output_file_handle.write(line)
    


def getMoleculeStretchedAlignments(q_cmap_input_path, alignment_label_channel, query_label_channel, labels_output_path, molecules_output_path):
    counter = 0
    with open(q_cmap_input_path) as q_cmap_file:
        line = skipHeader(q_cmap_file)
        while line:
            mol_id = line.split()[0]
            labels_dict, line = getNextMoleculeLabels(q_cmap_file,line)
            if mol_id not in xmap_infos:
                continue
            counter += 1
            if counter % 1000 == 0:
                print(counter)
            alignment_string_matches = parseAlignmentString(xmap_infos[mol_id][13])
            mol_positions, ref_positions = getRefAndQueryPositions(labels_dict, alignment_label_channel, alignment_string_matches, xmap_infos[mol_id][7], xmap_infos[mol_id][2])

            ius = interpolate.InterpolatedUnivariateSpline(mol_positions,ref_positions, k=1)
            corrected_query_labels = getStretchedLabelPositions(ius, labels_dict, query_label_channel)
            corrected_alignment_labels = getStretchedLabelPositions(ius, labels_dict, alignment_label_channel)

            if labels_output_path != "":
                with open(labels_output_path, 'a') as labels_output_file:
                    writeMolLabelsToFile(corrected_query_labels, xmap_infos[mol_id][2], 1, labels_output_file)

            if molecules_output_path != "":
                with open(molecules_output_path, 'a') as molecules_output_file:
                    writeMolRegionToFile(corrected_query_labels, corrected_alignment_labels, mol_id, xmap_infos[mol_id][2], 10, molecules_output_file)
                    



if not args.label_bed and not args.molecule_bed:
    print("Missing Argument. At least one of -lb and -mb should be selected.")
    sys.exit(1)
    

chrom_lengths, chrom_names = parseKeyFile(args.key)
parseXmap(args.xmap, args.min_confidence)
print("done parsing xmap")

findRCmapPositions(args.r_cmap, 1, r_cmap_positions)
print("done finding r_cmap_positions")

if args.output_dir == None:
    output_dir = os.getcwd()
else:
    output_dir = args.output_dir

if args.label_bed:
    labels_file_name = args.prefix+"conf."+str(args.min_confidence)+".labels_channel."+str(args.query_label_channel)+".bed"
    labels_output_path = os.path.join(output_dir,labels_file_name)
else:
    labels_output_path = ""

if args.molecule_bed:
    molecules_file_name = args.prefix+"conf."+str(args.min_confidence)+".molecules.bed"
    molecule_output_path = os.path.join(output_dir,molecules_file_name)
else:
    molecule_output_path = ""

getMoleculeStretchedAlignments(args.q_cmap, args.alignment_label_channel, args.query_label_channel, labels_output_path, molecule_output_path)

print("analysis complete")
