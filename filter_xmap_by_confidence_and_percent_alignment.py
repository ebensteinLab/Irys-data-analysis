import argparse
import os

parser = argparse.ArgumentParser(description='Filter BioNano xmap and q_cmap by alignment confidence and percent alignment out of molecule length')
parser.add_argument('xmap', help='input xmap file')
parser.add_argument('q_cmap', help='input q_cmap file')
parser.add_argument('-o', '--output', dest='output_dir', metavar='OUTPUT DIR', help='output files directory. default current working directory')
parser.add_argument('-c', '--conf', dest='min_confidence', metavar='MIN CONFIDENCE', type=int, default=12, help='min alignment confidence. default 12')
parser.add_argument('-a', '--alignment', dest='min_percent_alignment', metavar='MIN PERCENT ALIGNMENT', type=float, default=60, help='min percent of molecule length aligned. default 60')
parser.add_argument('-p', '--prefix', default="", help='output files prefix')

args = parser.parse_args()

mol_IDs = set()

def filter_xmap(xmap_input_path, xmap_output_path):
    total_molecules = 0
    molecules_kept = 0
    with open(xmap_input_path) as input_file:
        with open(xmap_output_path, 'w') as output_file:
            for line in input_file:
                if line.startswith('#'):
                    output_file.write(line)
                else:
                    total_molecules += 1
                    info = line.split()
                    qry_start = float(info[3])
                    qry_end = float(info[4])
                    qry_length = float(info[10])
                    percent_alignment = (abs(qry_end - qry_start) / qry_length) * 100
                    alignment_confidence = float(info[8])
                    if alignment_confidence >= args.min_confidence and percent_alignment >= args.min_percent_alignment:
                        output_file.write(line)
                        mol_IDs.add(info[1])
                        molecules_kept += 1

    print("finished parsing xmap, kept " + str(molecules_kept) + " out of " + str(total_molecules))


def filter_q_cmap(qcmap_input_path, qcmap_output_path):
    last_ID = ""
    counter = 0
    with open(qcmap_input_path) as input_file:
        with open(qcmap_output_path, 'w') as output_file:
            for line in input_file:
                if line.startswith('#'):
                    output_file.write(line)
                else:
                    info = line.split()
                    if info[0] in mol_IDs:
                        output_file.write(line)
                        if info[0] != last_ID:
                            counter += 1
                            if counter % 1000 == 0:
                                print("mol "+str(counter)+" out of " + str(len(mol_IDs)))
                            last_ID = info[0]


if args.output_dir == None:
    output_dir = os.getcwd()
else:
    output_dir = args.output_dir

xmap_output_file_name = args.prefix+"min.conf."+str(args.min_confidence)+".min.alignment."+str(args.min_percent_alignment)+".xmap"
q_cmap_output_file_name = args.prefix+"min.conf."+str(args.min_confidence)+".min.alignment."+str(args.min_percent_alignment)+"_q.cmap"
xmap_output_path = os.path.join(output_dir, xmap_output_file_name)
q_cmap_output_path = os.path.join(output_dir, q_cmap_output_file_name)


filter_xmap(args.xmap, xmap_output_path)
filter_q_cmap(args.q_cmap, q_cmap_output_path)
print("analysis complete")
                        
                    
                    
