import numpy as np
from scipy import interpolate

def parseXmap (xmap_input_path, min_confidence):
    #read alignment information for all molecules that aligned to the reference
    with open(xmap_input_path) as input_file:
        for line in input_file:
            if line[0] != "#":
                line = line.strip()
                info = line.split()
                if float(info[8]) >= min_confidence:
                    ID = info[1]
                    xmap_infos[ID] = info


def findQCmapPositions (cmap_input_path, label_channel, posDict):
    #read label position information for all molecules that aligned to the reference
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


def correctMoleculeLabelPositions (molID, refID, alignment_dict, query_dict):
    #return the correct positions, start and end
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
        while found == False and i<len(alignment_dict[molID]):
            if alignment_dict[molID][i][1] == tup[1]:
                queryPos = alignment_dict[molID][i][0]
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

    correctPositions = []

    for pos in query_dict[molID]:
        correctPos = ius(pos[0])
        if refID == "23":
            chrom = "chrX"
        elif refID == "24":
            chrom = "chrY"
        else:
            chrom = "chr"+refID
        if correctPos >= 0 and correctPos <= chrom_lengths[chrom]:
            correctPositions.append(int(correctPos))

            if correctPos > max_pos:
                max_pos = correctPos
            elif correctPos < min_pos:
                min_pos = correctPos

    return(correctPositions, min_pos, max_pos)


def processAllMolecules (output_file_path):
    with open(output_file_path, 'w') as output_file:
        mol_counter = 0
        for ID in xmap_infos:
            refID = xmap_infos[ID][2]
            #find the correct positions of the alignment label channel
            if ID in alignment_q_cmap_positions:
                alignment_correct_positions, alignment_min, alignment_max = correctMoleculeLabelPositions(ID, refID, alignment_q_cmap_positions, alignment_q_cmap_positions)
            #find the correct positions of the query label channel
            if ID in query_q_cmap_positions:
                query_correct_positions, query_min, query_max = correctMoleculeLabelPositions(ID, refID, alignment_q_cmap_positions, query_q_cmap_positions)
            else:
                query_min = alignment_min
                query_max = alignment_max
                query_correct_positions = []

            #find the molecule start and end positions - start = min of both, end = max of both
            if alignment_min <= query_min:
                min_pos = int(alignment_min) - 10
            else:
                min_pos = int(query_min) - 10

            if alignment_max >= query_max:
                max_pos = int(alignment_max) + 10
            else:
                max_pos = int(query_max) + 10

            #write to file
            alignment_positions = ",".join([str(pos) for pos in alignment_correct_positions])
            query_positions = ",".join([str(pos) for pos in query_correct_positions])
            line = "\t".join([ID, refID, str(min_pos), str(max_pos), alignment_positions, query_positions])+"\n"
            output_file.write(line)

            #update counter     
            mol_counter += 1
            if mol_counter % 1000 == 0:
                print(mol_counter)
                

                
xmap_input_path = "..//..//..//data//xmap_cmap//new_versions//april_2016//MoleculeQualityReport.confidence.12.min.alignment.60.no.duplicates.xmap"
q_cmap_input_path = "..//..//..//data//xmap_cmap//new_versions//april_2016//MoleculeQualityReport_q_cmap.confidence.12.min.alignment.60.no.duplicates.cmap"
r_cmap_input_path = "..//..//..//data//xmap_cmap//new_versions//april_2016//MoleculeQualityReport_r.cmap"
output_path = "april.conf.12.min.alignment.60.no.duplicates.positions.list.txt"

min_confidence = 12
alignment_label_channel = 1
query_label_channel = 2

chrom_lengths = {'chr1':249250621, 'chr2':243199373, 'chr3':198022430, 'chr4':191154276,
                 'chr5':180915260, 'chr6':171115067, 'chr7':159138663, 'chr8':146364022,
                 'chr9':141213431, 'chr10':135534747, 'chr11':135006516, 'chr12':133851895,
                 'chr13':115169878, 'chr14':107349540, 'chr15':102531392, 'chr16':90354753,
                 'chr17':81195210, 'chr18':78077248, 'chr19':59128983, 'chr20':63025520,
                 'chr21':48129895, 'chr22':51304566, 'chrX':155270560, 'chrY':59373566}

xmap_infos = {}
alignment_q_cmap_positions = {} #dict of positions where key=ID and val=list of tuples (posID,pos) for positions corresponding to alignment label channel
query_q_cmap_positions = {} #dict of positions where key=ID and val=list of tuples (posID,pos) for positions corresponding to query label channel
r_cmap_positions = {} #dict of positions where key=chrom and val=list of tuples (posID,pos)

#main
parseXmap(xmap_input_path, min_confidence)
print("done parsing xmap")

findQCmapPositions(q_cmap_input_path, alignment_label_channel, alignment_q_cmap_positions)
print("done parsing q_cmap for alignment channel")
findQCmapPositions(q_cmap_input_path, query_label_channel, query_q_cmap_positions)
print("done parsing q_cmap for query channel")
findRCmapPositions (r_cmap_input_path, alignment_label_channel, r_cmap_positions)
print("done finding r_cmap_positions")
processAllMolecules (output_path)
print("analysis complete")
