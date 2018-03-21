import argparse
import os

parser = argparse.ArgumentParser(description='Merge channel 1 and channel 2 molecule region BED files')
parser.add_argument('regions_ch1', help='molecule regions channel 1 input path')
parser.add_argument('regions_ch2', help='molecule regions channel 2 input path')
parser.add_argument('-o', '--output', dest='output_dir', metavar='OUTPUT DIR', help='output file directory. default current working directory')
parser.add_argument('-p', '--prefix', default="", help='output file prefix')

args = parser.parse_args()


def parse_mol_regions_bed(mol_regions_path):
    mol_regions_dict = {}
    with open(mol_regions_path) as mol_file:
        for line in mol_file:
            info = line.split()
            mol_ID = info[3]
            tup = (info[0],int(info[1]),int(info[2])) #chrom,start,end
            mol_regions_dict[mol_ID] = tup

    return mol_regions_dict


def combine_mol_regions(mol_regions_ch1_dict, mol_regions_ch2_dict, output_path):
    with open(output_path, 'w') as output_file:
        for mol_ID in mol_regions_ch1_dict:
            if mol_ID in mol_regions_ch2_dict:
                start = min(mol_regions_ch1_dict[mol_ID][1],
                            mol_regions_ch2_dict[mol_ID][1])
                end = max(mol_regions_ch1_dict[mol_ID][2],
                          mol_regions_ch2_dict[mol_ID][2])
                info = [mol_regions_ch2_dict[mol_ID][0],
                        str(start),str(end),mol_ID,"0","+"]
            else:
                info = [mol_regions_ch1_dict[mol_ID][0],
                        str(mol_regions_ch1_dict[mol_ID][1]),
                        str(mol_regions_ch1_dict[mol_ID][2]),
                        mol_ID, "0", "+"]
            line_to_print = "\t".join(info)+"\n"
            output_file.write(line_to_print)


if args.output_dir == None:
    output_dir = os.getcwd()
else:
    output_dir = args.output_dir

output_file_name = args.prefix+".molecule.regions.ch1.ch2.merged.bed"
output_file_path = os.path.join(output_dir, output_file_name)

molecules_dict_ch1 = parse_mol_regions_bed(args.regions_ch1)
molecules_dict_ch2 = parse_mol_regions_bed(args.regions_ch2)
combine_mol_regions(molecules_dict_ch1, molecules_dict_ch2, output_file_path)
