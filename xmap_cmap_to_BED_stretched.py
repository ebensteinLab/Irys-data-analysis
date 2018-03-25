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

def parse_key_file(key_file_path):
    '''parse the chromosome key file created when digesting a fasta reference file into a cmap''' 
    chrom_lengths = {}
    chrom_names = {}
    with open(args.key) as key_file:
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
                    
def findQCmapPositions (cmap_input_path, label_channel, posDict):
    '''read label position information for all molecules that aligned to the reference'''
    molecule_counter = 0
    with open(cmap_input_path) as input_file:
        prevID = ""
        for line in input_file:
            if line[0] != "#":
                line = line.strip()
                info = line.split()
                ID = info[0]
                if ID in xmap_infos: #if molecule aligned
                    if ID != prevID: #check if new molecule is starting
                        counter = 0
                        molecule_counter += 1
                        if molecule_counter % 1000 == 0:
                            print(molecule_counter)
                    if int(info[4]) == label_channel: #if label channel is the one we are querying
                        pos = float(info[5])
                        counter+=1
                        siteID = str(counter)
                        tup = (pos,siteID)
                        if ID in posDict:
                            posDict[ID].append(tup)
                        else:
                            posDict[ID] = [tup]
                prevID = ID
                    

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


def correctMoleculeLabelPositions (molID, refID, chrom_lengths, chrom_names):
    #make the connection between the molecule label to the reference label
    alignment = xmap_infos[molID][13]
    matches = alignment[1:len(alignment)-1].split(")(")
    tup_matches = []
    for match in matches:
        splitPos = match.find(",")
        tup = (match[:splitPos],match[splitPos+1:])
        tup_matches.append(tup)

    refPositions = []
    queryPositions = []
    
    for tup in tup_matches:
        refPos = r_cmap_positions[refID][int(tup[0])-1]
        found = False
        i = 0
        while found == False and i<len(alignment_q_cmap_positions[molID]):
            if alignment_q_cmap_positions[molID][i][1] == tup[1]:
                queryPos = alignment_q_cmap_positions[molID][i][0]
                found = True
            else:
                i += 1

        if found == True:
            if xmap_infos[molID][7] == "+":
                if queryPos not in queryPositions and refPos not in refPositions:
                    refPositions.append(refPos)
                    queryPositions.append(queryPos)
            else:
                if queryPos not in queryPositions and refPos not in refPositions:
                    refPositions.insert(0,refPos)
                    queryPositions.insert(0,queryPos)

    ius = interpolate.InterpolatedUnivariateSpline(queryPositions,refPositions, k=1)

    min_pos = min(refPositions)
    max_pos = max(refPositions)
    
    for pos in query_q_cmap_positions[molID]:
        correctPos = ius(pos[0])
        chrom_name = chrom_names[refID]
        if correctPos >= 0 and correctPos <= chrom_lengths[chrom_name]:
            if refID in correct_positions:
                correct_positions[refID].append(correctPos)
            else:
                correct_positions[refID] = [correctPos]

            if correctPos > max_pos:
                max_pos = correctPos
            elif correctPos < min_pos:
                min_pos = correctPos

    tup = (chrom_name,min_pos,max_pos)
    molRegions[molID] = tup


def correctAllMoleculePositions (chrom_lengths, chrom_names):
    mol_counter = 0
    for ID in xmap_infos:
        if ID in query_q_cmap_positions:
            refID = xmap_infos[ID][2]
            correctMoleculeLabelPositions(ID, refID, chrom_lengths, chrom_names)
            mol_counter += 1
            if mol_counter % 1000 == 0:
                print(mol_counter)


def createBED (output_file_path, chrom_names):
    with open(output_file_path, 'w') as output_file:
        for ref in correct_positions:
            chrom_name = chrom_names[ref]
            for pos in correct_positions[ref]:
                start = int(pos)-1
                end=int(pos)+2
                lineToPrint = chrom_name+"\t"+str(start)+"\t"+str(end)+"\n"
                output_file.write(lineToPrint)


def printMolRegionsBed(output_file_path):
    with open(output_file_path, 'w') as output_file:
        for ID,tup in molRegions.items():
            lineToPrint = tup[0]+"\t"+str(int(tup[1])-10)+"\t"+str(int(tup[2])+10)+"\t"+ID+"\t"+str(0)+"\t+\n"
            output_file.write(lineToPrint)


if not args.label_bed and not args.molecule_bed:
    print("Missing Argument. At least one of -lb and -mb must be selected.")
    sys.exit(1)            
            
chrom_lengths, chrom_names = parse_key_file(args.key)
parseXmap(args.xmap, args.min_confidence)
print("done parsing xmap")

if args.alignment_label_channel != args.query_label_channel:
    findQCmapPositions(args.q_cmap, args.alignment_label_channel, alignment_q_cmap_positions)
    findQCmapPositions(args.q_cmap, args.query_label_channel, query_q_cmap_positions)
else:
    findQCmapPositions(args.q_cmap, args.alignment_label_channel, alignment_q_cmap_positions)
    query_q_cmap_positions = alignment_q_cmap_positions.copy()
print("done finding q_cmap positions")

findRCmapPositions(args.r_cmap, 1, r_cmap_positions)
print("done finding r_cmap_positions")

correctAllMoleculePositions(chrom_lengths, chrom_names)
print("done calculating genomic positions")

if args.output_dir == None:
    output_dir = os.getcwd()
else:
    output_dir = args.output_dir

if args.label_bed:
    labels_file_name = args.prefix+"conf."+str(args.min_confidence)+".labels_channel."+str(args.query_label_channel)+".bed"
    labels_output_path = os.path.join(output_dir,labels_file_name)
    createBED(labels_output_path, chrom_names)
    print("done writing BED file")

if args.molecule_bed:
    molecules_file_name = args.prefix+"conf."+str(args.min_confidence)+".molecules_channel."+str(args.query_label_channel)+".bed"
    molecule_output_path = os.path.join(output_dir,molecules_file_name)
    printMolRegionsBed(molecule_output_path)
    print("done writing mol regions BED")

print("analysis complete")
