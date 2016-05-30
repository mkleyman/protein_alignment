import pandas
import csv

def read_in_expression(filename):
    '''

    :param filename:
    :return: table of expression with time values sorted
    '''
    table = pandas.read_csv(filename, header= None, index_col=0)
    sort_by = table.index[0]
    time_values = list(table.iloc[0])
    time_values.sort()
    sorted_table = table.T.sort_values(by=sort_by).T
    return sorted_table.drop(table.index[0]),  time_values

def read_in_homologs(filename, homolog_col, organism_col, protein_col, reference,
                     comparison):
    '''

    :param filename:
    :param homolog_col: index of column with homolog id
    :param organism_col: index of column with organism name
    :param protein_col: index of column with protein id
    :param reference: name of reference organism
    :param comparison: name of comparison organism
    :return: ref_dict: dictionary mapping reference protein to homolog id,
    comparison_dict: dictionary mapping homolog ids to comparison protein
    '''
    #maps refrence protiens to homologs
    ref_dict = {}
    #maps the homologs to the comparions protiens
    comparison_dict = {}
    with open(filename, 'rb') as homolog_file:
        homolog_reader = csv.reader(homolog_file, delimiter = "\t")
        homolog_reader.next()
        for row in homolog_reader:
            if row[organism_col] == reference:
                ref_dict[row[protein_col]] = row[homolog_col]
            elif row[organism_col] == comparison:
                comparison_dict[row[homolog_col]] = row[protein_col]
    return ref_dict, comparison_dict


def reference_to_comparison_map(ref_dict, comparison_dict, ref_protien_list,
                                comparison_protien_set):
    '''

    :param ref_dict: dictionary mapping reference protein to homolog id
    :param comparison_dict: dictionary mapping homolog ids to comparison protein
    :param ref_protien_list: list of relevant proteins from the reference
    :param comparison_protien_set: set of relevant proteins from the comparison
    :return: homolog dict: dictionary mapping reference protien to single comaprison
    homolog, homologs, set of comparison homologs
    '''
    homologs = set()
    homolog_dict = {}
    for ref_protein in ref_protien_list:
        if ref_protein in ref_dict:
            homolog = ref_dict[ref_protein]
            if homolog in comparison_dict:
                comp_protein = comparison_dict[homolog]
                if comp_protein in comparison_protien_set:
                    homolog_dict[ref_protein] = comp_protein
                    homologs.add(comp_protein)
    return homolog_dict, homologs

def reference_to_comparison_listmap(ref_dict, comparison_dict, ref_protien_list,
                                    comparison_protien_set):
    '''

    :param ref_dict: dictionary mapping reference protein to homolog id
    :param comparison_dict: dictionary mapping homolog ids to comparison protein
    :param ref_protien_list: list of relevant proteins from the reference
    :param comparison_protien_set: set of relevant proteins from the comparison
    :return: homolog dict: dictionary mapping reference protien to list of comaprison
    homolog, homologs, set of comparison homologs
    '''
    homolog_dict = {}
    homologs = set()
    for ref_protein in ref_protien_list:
        if ref_protein in ref_dict:
            homolog = ref_dict[ref_protein]
            if homolog in comparison_dict:
                comp_protein = comparison_dict[homolog]
                if comp_protein in comparison_protien_set:
                    homologs.add(comp_protein)
                    if ref_protein in homolog_dict:
                        homolog_dict[ref_protein].append(comp_protein)
                    else: homolog_dict[ref_protein] = [comp_protein]
    return homolog_dict, homologs

def read_in_GO(filename, proteins):
    '''

    :param filename:
    :param proteins: collection of relevant protiens
    :return: dictionary mapping protein id to lsit of corresponding GO ids
    '''
    protien_go_dict = {}
    go_protein_dict = {}
    with open(filename, 'rb') as homolog_file:
        go_reader = csv.reader(homolog_file, delimiter="\t")
        for row in go_reader:
            if row[0] in proteins:
                go_list = row[6].split(";")
                protien_go_dict[row[0]] = go_list
                for go in go_list:
                    if go not in go_protein_dict:
                        go_protein_dict[go] = set()
                    go_protein_dict[go].add(row[0])


    return protien_go_dict, go_protein_dict

def merge_GO(reference_protien_go_dict, reference_go_protein, homolog_dict, comparison_protein_go_dict):
    go_protein_dict = {}
    for protein in reference_protien_go_dict:
        homolog  = homolog_dict[protein]
        reference_protien_go_dict[protein] = reference_protien_go_dict[protein].union(comparison_protein_go_dict[homolog])
        for go in reference_protien_go_dict[protein]:
            if go not in go_protein_dict:
                go_protein_dict[go] = set()
            go_protein_dict[go].add(protein)
    return reference_protien_go_dict, go_protein_dict








