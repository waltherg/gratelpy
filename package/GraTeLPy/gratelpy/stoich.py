import numpy as np
import networkx as nx

def get_substance_adjacency(alpha, beta):
    # alpa:
    # beta:
    
    # number of reactions = number of columns of alpha
    no_rxn = alpha.shape[1]
    # number of substance = number of rows of alpha
    no_sub = alpha.shape[0]

    # check
    if no_rxn != beta.shape[1] or no_sub != beta.shape[0]:
        raise

    # substance adjacency matrix
    subs_adj = []

    # build subs_adj matrix row-by-row
    for sub_i in range(no_sub):
        subs_adj_row = np.zeros((no_sub,), dtype=np.int)
        for rxn_i in range(no_rxn):
            if alpha[sub_i][rxn_i] > 0:
                for adj_sub_i in range(no_sub):
                    if beta[adj_sub_i][rxn_i] > 0:
                        subs_adj_row[adj_sub_i] = 1
        subs_adj.append(subs_adj_row)

    # return numpy array of adjacency matrix
    return np.array(subs_adj)
        
def get_graph_stoich(alpha, beta, complex_names = None, constant_names = None):

    # stoichiometry matrix
    stoich = beta - alpha
    # rank of stoichiometry matrix
    stoich_rank = np.linalg.matrix_rank(stoich)
    # number of reactions = number of columns of stoich
    no_rxn = stoich.shape[1]
    # number of substance = number of rows of stoich
    no_sub = stoich.shape[0]

    # make directed graph
    G = nx.DiGraph()

    # add nodes to graph: add all reactions 'w' and substances 's'
    if constant_names == None:
        G.add_nodes_from(['w'+str(n) for n in range(1,no_rxn+1)], bipartite=1)
    else:
        G.add_nodes_from([constant_names[n-1] for n in range(1,no_rxn+1)], bipartite=1)

    if complex_names == None:
        G.add_nodes_from(['s'+str(n) for n in range(1,no_sub+1)], bipartite=0)
    else:
        G.add_nodes_from([complex_names[n-1] for n in range(1,no_sub+1)], bipartite=0)

    # add ('s','w') edges, i.e. positive element in alpha
    for row_i in range(no_sub):
        curr_row = alpha[row_i]
        for col_i in range(no_rxn):
            curr_coeff = curr_row[col_i]
            if curr_coeff > 0:
                if constant_names == None and complex_names == None:
                    G.add_edge('s'+str(row_i+1), 'w'+str(col_i+1), coeff=curr_coeff)
                elif constant_names == None and complex_names != None:
                    G.add_edge(complex_names[row_i], 'w'+str(col_i+1), coeff=curr_coeff)
                elif constant_names != None and complex_names == None:
                    G.add_edge('s'+str(row_i+1), constant_names[col_i], coeff=curr_coeff)
                elif constant_names != None and complex_names != None:
                    G.add_edge(complex_names[row_i], constant_names[col_i], coeff=curr_coeff)
                else:
                    raise

    # add ('w','s') edges, i.e. positive element in beta
    for row_i in range(no_sub):
        curr_row = beta[row_i]
        for col_i in range(no_rxn):
            curr_coeff = curr_row[col_i]
            if curr_coeff > 0:
                if constant_names == None and complex_names == None:
                    G.add_edge('w'+str(col_i+1), 's'+str(row_i+1), coeff=curr_coeff)
                elif constant_names == None and complex_names != None:
                    G.add_edge('w'+str(col_i+1), complex_names[row_i], coeff=curr_coeff)
                elif constant_names != None and complex_names == None:
                    G.add_edge(constant_names[col_i], 's'+str(row_i+1), coeff=curr_coeff)
                elif constant_names != None and complex_names != None:
                    G.add_edge(constant_names[col_i], complex_names[row_i], coeff=curr_coeff)
                else:
                    raise

    return G, stoich, stoich_rank