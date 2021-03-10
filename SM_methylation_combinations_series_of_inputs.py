# this script gets as input: 1) directory with several molecules labels list files
# 2) promoters list. 3) list of enhancers with their location and the gene they are linked to. 4) output directory.
# additional parameters: number of chromosomes to process, minimal number of molecules per pair, minimal E-P distance.
# The script looks for molecules that contain an enhancer-promoter pair, and records the
# methylation states of the enhancer and promoter on the molecule. then the methylation states of that pair from all
# the molecules containing it are summed. the final output is a tab separated txt file with the fields: Gene_sym, enh_ID,
# number of molecules with: enh_m_prom_m, enh_um_prom_um, enh_um_prom_m, ehn_m_prom_um, chromosome, total number of
# molecules, and E-P distance. this version is for cmd.

import os
import pandas as pd
from time import perf_counter
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='get distributions of methylation combinations in enhancer-promoters pairs')
parser.add_argument('-pr', '--prom', help='promoters file')
parser.add_argument('-en', '--enh', help='enhancers file')
parser.add_argument('-m', '--molcs', dest='input_dir', metavar='INPUT DIR', help='input files directory')
parser.add_argument('-o', '--output', dest='output_dir', metavar='OUTPUT DIR', help='output files directory')
parser.add_argument('-n', '--chromnum', type=int, default=24, help='number of chromosomes to process')
parser.add_argument('-mc', '--mincov', type=int, default=10, help='minimal molecules number per pair')
parser.add_argument('-d', '--distance', type=int, default=5000, help='minimal distance between enh and prom [bp]')
args = parser.parse_args()

class Promoters:
    def __init__(self, chrom, start, end, gene_sym, TCGA_count):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.gene_sym = gene_sym
        self.TCGA_count = TCGA_count


class Enhancers:
    def __init__(self, chrom, start, end, gene_sym, enh_ID, TCGA_count):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.gene_sym = gene_sym
        self.enh_ID = enh_ID
        self.TCGA_count = TCGA_count


class Molecules:
    def __init__(self, mol_ID, chrom, start, end, labels_list):
        self.mol_ID = mol_ID
        self.chrom = chrom
        self.start = start
        self.end = end
        self.labels_list = labels_list


class Labels:
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end


def promoters_file_to_dict(prom_file_path, num_of_desired_chroms):
    print("promoters_file_to_dict", flush=True)
    chromosomes = canonical_chromosomes_list(num_of_desired_chroms)
    prom_dict = {chromosome: {} for chromosome in chromosomes}
    with open(prom_file_path, 'r') as prom_file:
        for prom_line in prom_file:
            chrom, start, end, gene_sym, TCGA_count = prom_line.split("\t")
            start = int(start)
            end = int(end)
            prom = Promoters(chrom, start, end, gene_sym, TCGA_count)
            prom_dict[chrom].update({gene_sym: prom})
        return prom_dict


def enhancers_file_to_dict(enhancers_file_path, num_of_desired_chroms):
    print("enhancers_file_to_dict", flush=True)
    chromosomes = canonical_chromosomes_list(num_of_desired_chroms)
    enh_dict = {chromosome: {} for chromosome in chromosomes}
    # pairs_dict = {chromosome: [] for chromosome in chromosomes}
    pairs_dict = defaultdict(set)
    with open(enhancers_file_path, 'r') as enh_file:
        for enh_line in enh_file:
            chrom, start, end, gene_sym, enh_ID, TCGA_count = enh_line.split("\t")
            start = int(start)
            end = int(end)
            enh = Enhancers(chrom, start, end, gene_sym, enh_ID, TCGA_count)
            enh_dict[chrom].update({enh_ID: enh})
            pair_name = enh_ID + "_" + gene_sym
            pairs_dict[chrom].add(pair_name)
    return enh_dict, pairs_dict


def canonical_chromosomes_list(actual_num_chroms):
    chroms = 24 * [0]
    for i in range(1, 23):
        char = "chr" + str(i)
        chroms[i - 1] = char
    chroms[22] = "chrX"
    chroms[23] = "chrY"
    return chroms[0:actual_num_chroms]


def molecules_file_to_dict(molcs_file_path, num_of_desired_chroms):
    print("molecules_file_to_dict", flush=True)
    chromosomes = canonical_chromosomes_list(num_of_desired_chroms)
    molcs_dict = {chromosome: {} for chromosome in chromosomes}
    with open(molcs_file_path) as molcs_file:
        for molc_line in molcs_file:
            mol_ID, chrom, start, end, labels_string = molc_line.split("\t")
            start = int(start) - 500
            end = int(end) + 500
            labels_list = labels_string_to_labels_class(labels_string, chrom)
            molc = Molecules(mol_ID, chrom, start, end, labels_list)
            molcs_dict[chrom].update({molc.mol_ID: molc})
    return molcs_dict


def labels_string_to_labels_class(labels_string, chrom):
    labels_list = []
    if labels_string != '\n':
        for i in labels_string.split(","):
            start = int(i) - 500
            end = int(i) + 500
            label = Labels(chrom, start, end)
            labels_list.append(label)
    return labels_list


def find_overlaps(query_molcs_dict_chrom, subject_dict_chrom, molcs_elements_overlaps_dict, chrom):
    print("find_overlaps", flush=True)
    start_time = perf_counter()
    overlaps_list_of_tuples = []
    for q in query_molcs_dict_chrom:
        for s in subject_dict_chrom:
            if query_molcs_dict_chrom[q].start <= subject_dict_chrom[s].start and query_molcs_dict_chrom[q].end >= \
                    subject_dict_chrom[s].end:
                overlaps_list_of_tuples.append((q, s))
    overlaps_pd_df = pd.DataFrame(overlaps_list_of_tuples, columns=['mol', 'subject'])
    molcs_elements_overlaps_dict[chrom] = overlaps_pd_df
    print(perf_counter() - start_time, flush=True)
    return molcs_elements_overlaps_dict


def match_EP_on_the_same_mol(molcs_promoters_overlaps_chr, molcs_enhancers_ovelaps_chr, EP_on_same_mol_dict, chrom):
    print("match_EP_on_the_same_mol", flush=True)
    start_time = perf_counter()
    EP_on_same_mol_chr = molcs_promoters_overlaps_chr.merge(molcs_enhancers_ovelaps_chr, on='mol')
    EP_on_same_mol_dict[chrom] = EP_on_same_mol_chr
    print(perf_counter() - start_time, flush=True)
    return EP_on_same_mol_dict


def filter_pairs_on_molcs(same_molecule_EP, pairs_list, filtered_EP_on_same_mol_dict, chrom):
    print("filter_pairs_on_molcs", flush=True)
    start_time = perf_counter()
    same_molecule_EP['pair_name'] = None
    same_molecule_EP['pair_name'] = (same_molecule_EP['subject_y'] + "_" + same_molecule_EP['subject_x'])
    same_molecule_EP.query('pair_name in @pairs_list', inplace=True)
    filtered_EP_on_same_mol_dict[chrom] = same_molecule_EP
    # print(same_molecule_EP.head(3))
    print(perf_counter() - start_time, flush=True)
    return filtered_EP_on_same_mol_dict


def overlaps_labels_elements(mol_ID, element_ID, molcs_dict, elements_dict):
    labels_count = 0
    molecule = molcs_dict[mol_ID]
    element = elements_dict[element_ID]
    for lab in molecule.labels_list:
        if ((element.start <= lab.start <= element.end and element.start <= lab.end <= element.end) or
                (lab.start < element.start and (element.start + 500) <= lab.end < element.end) or
                (lab.end > element.end and element.start < lab.start <= element.end - 500)):
            labels_count += 1
    return labels_count


def get_distance(gene_sym, enh_ID, prom_dict, enh_dict):
    prom_start = prom_dict[gene_sym].start
    prom_end = prom_dict[gene_sym].end
    enh_start = enh_dict[enh_ID].start
    enh_end = enh_dict[enh_ID].end
    distance = 0
    if prom_start > enh_end:
        distance = prom_start - enh_end
    elif enh_start > prom_end:
        distance = enh_start - prom_end
    return distance


def get_methylation(filtered_EP_on_same_mol, prom_dict, enh_dict, molcs_dict, chrom):
    print("get_methylation", flush=True)
    start_time = perf_counter()
    filtered_EP_on_same_mol['enh_m_prom_m'] = 0
    filtered_EP_on_same_mol['enh_um_prom_um'] = 0
    filtered_EP_on_same_mol['enh_um_prom_m'] = 0
    filtered_EP_on_same_mol['enh_m_prom_um'] = 0
    filtered_EP_on_same_mol['prom_labels'] = [overlaps_labels_elements(mol, gene_sym, molcs_dict, prom_dict)
                                              for (mol, gene_sym) in zip(filtered_EP_on_same_mol['mol'],
                                                                         filtered_EP_on_same_mol['subject_x'])]
    filtered_EP_on_same_mol['enh_labels'] = [overlaps_labels_elements(mol, enh_ID, molcs_dict, enh_dict) for
                                             (mol, enh_ID) in zip(filtered_EP_on_same_mol['mol'],
                                                                  filtered_EP_on_same_mol['subject_y'])]
    filtered_EP_on_same_mol.loc[(filtered_EP_on_same_mol['prom_labels'] == 0) &
                                (filtered_EP_on_same_mol['enh_labels'] == 0), 'enh_m_prom_m'] = 1
    filtered_EP_on_same_mol.loc[(filtered_EP_on_same_mol['prom_labels'] > 0) &
                                (filtered_EP_on_same_mol['enh_labels'] > 0), 'enh_um_prom_um'] = 1
    filtered_EP_on_same_mol.loc[(filtered_EP_on_same_mol['prom_labels'] == 0) &
                                (filtered_EP_on_same_mol['enh_labels'] > 0), 'enh_um_prom_m'] = 1
    filtered_EP_on_same_mol.loc[(filtered_EP_on_same_mol['prom_labels'] > 0) &
                                (filtered_EP_on_same_mol['enh_labels'] == 0), 'enh_m_prom_um'] = 1
    filtered_EP_on_same_mol = filtered_EP_on_same_mol.rename(columns={"subject_x": "gene_sym", "subject_y": "enh_ID"})
    filtered_EP_on_same_mol = filtered_EP_on_same_mol.drop(columns=['mol', 'pair_name', 'prom_labels', 'enh_labels'])
    # print(filtered_EP_on_same_mol.head)
    results_by_pair = filtered_EP_on_same_mol.groupby(by=['gene_sym', 'enh_ID'], as_index=False).sum()

    results_by_pair['total_molcs'] = None
    results_by_pair['total_molcs'] = (results_by_pair['enh_m_prom_m'] + results_by_pair['enh_um_prom_um'] +
                                      results_by_pair['enh_um_prom_m'] + results_by_pair['enh_m_prom_um'])
    results_by_pair['chrom'] = None
    results_by_pair['chrom'] = chrom
    #print(results_by_pair.head)
    results_by_pair['EP_distance'] = None
    results_by_pair['EP_distance'] = [get_distance(gene_sym, enh_ID, prom_dict, enh_dict) for
                                      (gene_sym, enh_ID) in zip(results_by_pair['gene_sym'], results_by_pair['enh_ID'])]
    print(results_by_pair.head, flush=True)
    print(perf_counter() - start_time, flush=True)
    return results_by_pair


def main():
    start_time1 = perf_counter()

    prom_file_path = args.prom
    enhancers_file_path = args.enh
    inputs_dir = args.input_dir
    out_dir = args.output_dir
    min_EP_distance = args.distance
    min_pair_coverage = args.mincov
    num_of_desired_chroms = args.chromnum

    chromosomes = canonical_chromosomes_list(num_of_desired_chroms)

    prom_dict = promoters_file_to_dict(prom_file_path, num_of_desired_chroms)
    enh_dict, enh_dict_pairs = enhancers_file_to_dict(enhancers_file_path, num_of_desired_chroms)

    input_dir_list = os.listdir(inputs_dir)
    for input_file in input_dir_list:
        input_file_path = os.path.join(inputs_dir, input_file)
        molcs_dict = molecules_file_to_dict(input_file_path, num_of_desired_chroms)
        molcs_promoters_overlaps_dict = {chromosome: {} for chromosome in chromosomes}
        molcs_enhancers_ovelaps_dict = {chromosome: {} for chromosome in chromosomes}
        EP_on_same_mol_dict = {chromosome: {} for chromosome in chromosomes}
        filtered_EP_on_same_mol_dict = {chromosome: {} for chromosome in chromosomes}

        results_all_chroms = pd.DataFrame([], columns=['gene_sym', 'enh_ID', 'enh_m_prom_m', 'enh_um_prom_um',
                                                       'enh_um_prom_m', 'enh_m_prom_um', 'chrom', 'total_molcs',
                                                       'EP_distance'])

        for chrom in chromosomes:
            print(chrom, flush=True)
            if molcs_dict[chrom] and prom_dict[chrom] and enh_dict[chrom]:
                molcs_promoters_overlaps_dict = find_overlaps(molcs_dict[chrom], prom_dict[chrom],
                                                              molcs_promoters_overlaps_dict, chrom)
                molcs_enhancers_ovelaps_dict = find_overlaps(molcs_dict[chrom], enh_dict[chrom],
                                                             molcs_enhancers_ovelaps_dict, chrom)
                EP_on_same_mol_dict = match_EP_on_the_same_mol(molcs_promoters_overlaps_dict[chrom],
                                                               molcs_enhancers_ovelaps_dict[chrom], EP_on_same_mol_dict, chrom)
                filtered_EP_on_same_mol_dict = filter_pairs_on_molcs(EP_on_same_mol_dict[chrom], enh_dict_pairs[chrom],
                                                                     filtered_EP_on_same_mol_dict, chrom)
                results_all_chroms = pd.concat([results_all_chroms, get_methylation(filtered_EP_on_same_mol_dict[chrom],
                                                                                    prom_dict[chrom], enh_dict[chrom],
                                                                                    molcs_dict[chrom], chrom)])
        results_all_chroms.query('EP_distance > @min_EP_distance and total_molcs > @min_pair_coverage', inplace=True)

        in_name = os.path.splitext(input_file)[0]
        out_name = "SM_results_" + in_name + ".txt"
        results_file_path = os.path.join(out_dir, out_name)

        results_all_chroms.to_csv(results_file_path, sep='\t', index=False)
        print("total time:", flush=True)
        print(perf_counter() - start_time1, flush=True)


if __name__ == '__main__':
    main()
