import networkx as nx
from networkx.algorithms import bipartite
from networkx.drawing import nx_pydot
import pydot
import os

from matplotlib import pyplot as plt
import matplotlib

def gratelpy_dot(G, positions=None, dictionary_complexes=None, dictionary_reactions=None, subgraph_components=None, filename=None, filename_print=None):
    
    if dictionary_complexes is not None:
        if all(type(k) is not type(int()) for k in dictionary_complexes.keys()):
            new_dict = {v:k for k,v in dictionary_complexes.iteritems()}
            dictionary_complexes = new_dict

    if dictionary_reactions is not None:
        if all(type(k) is not type(int()) for k in dictionary_reactions.keys()):
            new_dict = {v:k for k,v in dictionary_reactions.iteritems()}
            dictionary_reactions = new_dict

    pydot_graph = nx_pydot.to_pydot(G)

    if positions is not None:
        if not any('[' in k for k in positions.keys()):
            for node in pydot_graph.get_node_list():
                pos_list = positions[node.get_name()]
            
                node.set_pos('\"'+str(pos_list[0])+','+str(pos_list[1])+'!\"')

    for node in pydot_graph.get_node_list():
        if 'w' == node.get_name()[0]:
            node.set_shape('\"box\"')
        elif 's' == node.get_name()[0]:
            node.set_shape('\"circle\"')

    if dictionary_reactions is not None:
        for node in pydot_graph.get_node_list():
            node_name = node.get_name()
            node_index = int(node_name[1:])
            if 'w' == node_name[0]:
                node.set_name(dictionary_reactions[node_index-1].translate(None, '[]'))

    if dictionary_complexes is not None:
        for node in pydot_graph.get_node_list():
            node_name = node.get_name()
            node_index = int(node_name[1:])
            if 's' in node_name:
                node.set_name(dictionary_complexes[node_index-1].translate(None, '[]'))

    if positions is not None:
        if any('[' in k for k in positions.keys()):
            for node in pydot_graph.get_node_list():
                if node.get_shape() == '\"circle\"':
                    pos_list = positions['['+node.get_name()+']']
                elif node.get_shape() == '\"box\"':
                    pos_list = positions[node.get_name()]
                else:
                    print node.get_shape()
                    raise
            
                node.set_pos('\"'+str(pos_list[0])+','+str(pos_list[1])+'!\"')

    if dictionary_reactions is not None and dictionary_complexes is not None:
        for edge in pydot_graph.get_edge_list():
            edge_dest = edge.get_destination()
            edge_dest_index = int(edge_dest[1:])

            edge_orig = edge.get_source()
            edge_orig_index = int(edge_orig[1:])

            if edge_dest[0] == 'w':
                edge_dest_new = dictionary_reactions[edge_dest_index-1].translate(None, '[]')
            elif edge_dest[0] == 's':
                edge_dest_new = dictionary_complexes[edge_dest_index-1].translate(None, '[]')
            else:
                raise

            if edge_orig[0] == 'w':
                edge_orig_new = dictionary_reactions[edge_orig_index-1].translate(None, '[]')
            elif edge_orig[0] == 's':
                edge_orig_new = dictionary_complexes[edge_orig_index-1].translate(None, '[]')
            else:
                raise

            edge_new = pydot.Edge(src=edge_orig_new, dst=edge_dest_new)
            pydot_graph.del_edge(edge_orig, dst=edge_dest)
            pydot_graph.add_edge(edge_new)

    # write dot file
    if filename is not None:
        try:
            with open(filename) as df:
                print ''
                print 'WARNING. gratelpy_dot: dot file',df.name,'exists already. Will not overwrite hence doing nothing.'
                print ''
        except IOError as e:
            with open(filename, 'w') as df:
                df.write(pydot_graph.to_string())
                if filename_print is not None:# and filename_print.split('.')[-1].lower() == 'pdf': 
                    os.system('dot -Kfdp -n -Tpdf -o '+filename_print+' '+df.name)

    return pydot_graph

def gratelpy_draw(G, positions=None, dictionary_complexes=None, dictionary_reactions=None, filename=None, subgraph=None, rnsize = 1600, cnsize = 1400, dotfile=None):
    # draws entire graph or subgraph (subgraph being a fragment)
    # squares for reaction nodes, circles for complex nodes
    # inscribes s, w (complex and reaction nodes) labels into nodes

    # positions dictionary expected of the form {'si': [x_pos, y_pos], 'sj': ...}
    # supplied dictionary_complexes expected of the form {i: 'descriptive name of s(i+1)'}
    # supplied dictionary_reactions expected of the form {i: 'descriptive name of w(i+1)'}
    # filename expected as string including prefix

    font = {'family' : 'sans-serif','sans-serif':['Helvetica'],
            'weight' : 'normal',
            'size'   : 26}

    matplotlib.rc('font', **font)

    # generate positions of nodes if not supplied
    if positions is None:
        positions = nx.spring_layout(G)

    # generate and modify figure
    fig = plt.figure(num=None, figsize=(20,10), dpi=80, facecolor='w', edgecolor='k')
    fig_axis = fig.add_subplot(111)
    fig_axis.axis('off')

    # generate separate graphs for both complexes and reactions so we can draw them differently
    substance_nodes = []
    reaction_nodes = []
    for n in G.nodes():
        if G.node[n]['bipartite']==0 and n not in substance_nodes:
            substance_nodes.append(n)
        if G.node[n]['bipartite']==1 and n not in reaction_nodes:
            reaction_nodes.append(n)
    substance_graph = nx.DiGraph()
    substance_graph.add_nodes_from(substance_nodes)
    reaction_graph = nx.DiGraph()
    reaction_graph.add_nodes_from(reaction_nodes)

    if dotfile is not None:
        dot_graph = nx.DiGraph()
        dot_graph.add_nodes_from(substance_nodes+reaction_nodes)
    else:
        dot_graph = None

    # if drawing subgraph, then generate specifically edges that are to be displayed
    if subgraph is not None:
        edges_graph = nx.DiGraph()
        edges_graph.add_nodes_from(substance_nodes+reaction_nodes)
        for el in subgraph:
            if len(el) == 2:
                # edge
                edges_graph.add_edge(el[0], el[1])
            else:
                if el[-1] == 'p':
                    edges_graph.add_edge(el[0], el[1])
                    edges_graph.add_edge(el[1], el[2])
                elif el[-1] == 'n':
                    edges_graph.add_edge(el[0], el[1])
                    edges_graph.add_edge(el[2], el[1])
                else:
                    raise

        if dotfile is not None:
            dot_graph.add_edges_from(edges_graph.edges())
    else:
        edges_graph = None
        if dotfile is not None:
            dot_graph.add_edges_from(edge for edge in G.edges() if all(node in substance_nodes+reaction_nodes for node in edge))

    # generate complex labels
    if dictionary_complexes is None:
        complex_labels = {}
        for n in substance_graph.nodes():
            complex_labels[n] = str(n)
    else:
        complex_labels = {}
        for n in substance_graph.nodes():
            complex_labels[n] = dictionary_complexes[int(n[1:])-1].translate(None, '[]')

    # generate reaction labels
    if dictionary_reactions is None:
        reaction_labels = {}
        for n in reaction_graph.nodes():
            reaction_labels[n] = str(n)
    else:
        reaction_labels = {}
        for n in reaction_graph.nodes():
            reaction_labels[n] = dictionary_reactions[int(n[1:])-1].translate(None, '[]')
    
    # draw substance and reaction nodes
    nx.draw_networkx_nodes(substance_graph, positions, node_shape='o', node_size=cnsize, node_color='white')
    nx.draw_networkx_nodes(reaction_graph, positions, node_shape='s', node_size=rnsize, node_color='white')

    if subgraph is None:
        nx.draw_networkx_edges(G, positions, width=2)
    else:
        nx.draw_networkx_edges(edges_graph, positions, width=2)

    nx.draw_networkx_labels(substance_graph, positions, complex_labels, font_size=26)
    nx.draw_networkx_labels(reaction_graph, positions, reaction_labels, font_size=26)

    # write dot file
    if dotfile is not None:
        try:
            with open(dotfile) as df:
                print ''
                print 'WARNING. gratelpy_draw: dot file',df.name,'exists already. Will not overwrite hence doing nothing.'
                print ''
        except IOError as e:
            nx.write_dot(dot_graph, dotfile)

    return fig
