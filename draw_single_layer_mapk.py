from stoich import get_graph_stoich
from parse_mechanism import get_network_from_mechanism
from drawing import gratelpy_draw, gratelpy_dot
from matplotlib import pyplot as plt
from fragments import get_sensible_fragments
from subgraphs import get_all_valid_subgraphs

#mechanism_file = '/usr/users/cbu/waltherg/JIC/Dropbox/GratelPy/mechanisms/reversible_substrate_inhibition.txt'
mechanism_file = '../mechanisms/single_layer_mapk_mechanism.txt'

base_name = 'single_layer_mapk'

alpha, beta, dict_complexes, dict_constants, dict_complexes_reverse, dict_constants_reverse = get_network_from_mechanism(mechanism_file, 9)
print dict_complexes
print dict_constants
G, stoich, stoich_rank = get_graph_stoich(alpha, beta)

dx = 1
pos_mech = {
    '[E1]': [3,4],
    'k1': [0.5,3.5], 'k4': [4.5, 3.5],
    '[AE1]': [1,3], 'k3': [2,3], '[Ap]': [3,3], '[ApE1]': [5,3], 'k6': [6,3],
    'k2': [0,3],'k5': [3.5, 2.5],
    '[A]': [0,2], '[App]': [6,2],
    'k10': [1.5, 1.5], 'k7': [4.5, 1.5],
    'k12': [0,1], '[ApE2]': [1,1], 'k9': [3,1], '[AppE2]': [4,1],
    'k11': [1.5, 0.5], 'k8': [5.5, 0.5],
    '[E2]': [3,0]
    }

pos = {}
for s in dict_complexes:
    pos['s'+str(dict_complexes[s]+1)] = pos_mech[s]
for w in dict_constants:
    pos['w'+str(dict_constants[w]+1)] = pos_mech[w]

dot_graph = gratelpy_dot(G, positions=pos_mech, dictionary_complexes=dict_complexes_reverse, dictionary_reactions=dict_constants_reverse,filename='single_layer_mapk.dot')
