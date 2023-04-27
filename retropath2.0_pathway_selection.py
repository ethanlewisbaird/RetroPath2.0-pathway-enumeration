import argparse as ap
import pandas as pd
import csv
import time

start_time = time.time()

# function that processes data into useful dictionaries and lists
def process_data(data):
    rxns = {}
    products_sets_list = []
    substrates_list = []
    for i, row in data.iterrows():
        print('creating a dictionary of reaction ID and reaction SMILES, a list of product sets and a list of substrates')
        if row['Transformation ID'] not in rxns:
            rxns[row['Transformation ID']] = []
        rxns[row['Transformation ID']].append(row['Reaction SMILES'])
        products_sets_list.append(row['Reaction SMILES'].split('>>')[1])
        substrates_list.append(row['Reaction SMILES'].split('>>')[0])
    products_sets_list = list(set(products_sets_list))
    substrates_list = list(set(substrates_list))
    return rxns, products_sets_list, substrates_list

# function that changes the value of specific datframe location in the
# hypergraph matrix to represent the existence of a reaction between the
# substrates (matrix row labels) and the product (matrix column label)
def update_matrix(rxns, matrix):
    for i, (trans_id, rxn_smiles) in enumerate(rxns.items()):
        for rxn in rxn_smiles:
            print('updating matrix with binary information displaying what product sets are required for each substrate')
            substrate, products = rxn.split('>>')
            products = '.'.join(sorted(products.split('.')))
            matrix.at[products, substrate] = 1

# depth-first search function that starts with the target metabolite  and
# searches within the column in matrix that corresponds to the current metabolite
# for a 1. Then recursively calls itself (this function) to the substrates of the
# reaction that produces the current metabolite until all paths to the target metabolome
# have been explored
def dfs_paths(max_length, matrix, metabolite, target_metabolome, current_path=None):
    if current_path is None:
        current_path = []
    current_path.append(metabolite)
    if metabolite in target_metabolome:
        print('Path found')
        with open('paths.csv', mode='a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(current_path)
    if len(current_path) > max_length:
        print('Pathway skipped as it is greater than the max length')
    else:
        for index, row in matrix.iterrows():
            if row[metabolite] == 1:
                print('Exploring next node')
                for product in row.name.split('.'):
                    if 'C' in product or 'c' in product and product is not 'O=C=O':
                        if product in matrix.columns and product not in current_path:
                            dfs_paths(max_length, matrix, product, target_metabolome, current_path.copy())
                        else:
                            print(f'''Skipping node because it is not a precursor to any 
other nodes or has already been visited in the current path''')
                    else:
                        print(f'''Skipping node {product} because it is a currency metabolite''')
    return current_path


def main():
#   command line parser    
    parser = ap.ArgumentParser(description='Pathway selection for RetroPath2.0')
    parser.add_argument("fileName", help="name of file containing data")
    parser.add_argument("max_length", help="maximum number of reactions in a pathway", default=10)
    args = parser.parse_args()

#   create a dataframe for the retropath2.0 output data
    print('creating dataframe for retropath2.0 output data')
    data = pd.read_csv(args.fileName)

#   generating list of nodes
    print('generating list of node and hyperedge names')
    substrates = data['Substrate SMILES'].unique().tolist()
    products = data['Product SMILES'].unique().tolist()
    nodes = list(set(substrates + products))

#   process the data to generate useful dictionaries and lists
    rxns, products_sets_list, substrates_list = process_data(data)

#   create a matrix to represent the hypergraph and initialise all values as 0
    print('creating a matrix to represent the metabolic network as a hypergraph')
    matrix = pd.DataFrame(0, index=products_sets_list, columns=nodes)

#   update the matrix with the binary information giving the hypergraph directionality
    update_matrix(rxns, matrix)
    print('matrix generated')
    matrix.to_csv('matrix.csv')

#   call the depth-first search function to enumarate all possible pathways
    starting_metabolite = data.loc[0, 'Substrate SMILES']
    target_metabolome = set(data.loc[data['In Sink'] == 1, 'Product SMILES'])
    max_length = int(args.max_length)
    dfs_paths(max_length, matrix, starting_metabolite, target_metabolome)

    print('current method complete. results available in paths.csv in the current directory')


    end_time = time.time()
    total_time = end_time - start_time

    print('Total execution time:', total_time, 'seconds')

if __name__ == '__main__':
    main()

