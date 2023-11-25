import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

tf_list = []
with open('./TFlist.txt', 'r') as f:
    for line in f:
        tf_list.append(line.strip())

highlight_list = []
with open('./Highlight_TF.txt', 'r') as f:
    for line in f:
        highlight_list.append(line.strip().split('\t'))


save_folder = './'


df = pd.read_csv('./replicates.txt', sep='\t', index_col=0)


G = nx.Graph()


for tf in tf_list:
    G.add_node(tf)


for i in range(len(df.index)):
    for j in range(len(df.columns)):
        if df.iloc[i, j] == 1:
            G.add_edge(df.index[i], df.columns[j])


angles = [2 * i * np.pi / len(tf_list) for i in range(len(tf_list))]
pos = {tf_list[i]: (np.cos(angles[i]), np.sin(angles[i])) for i in range(len(tf_list))}



plt.figure(figsize=(10, 10))
nx.draw(G, pos=pos, node_color='#cccccc', font_color='k', node_size=100, with_labels=False, font_size=1.6, edge_color='#cccccc', width=1, alpha=0.3)


degrees = dict(G.degree())
zero_degrees = [node for node in degrees.keys() if degrees[node] == 0]
nx.draw_networkx_nodes(G, pos=pos, nodelist=zero_degrees, node_color='white', node_size=100)


for i in range(len(highlight_list)):
    if i == 0:
        color = 'r'
    else:
        color = 'b'
    for tf in highlight_list[i]:
        if tf in G.nodes():
            nx.draw_networkx_nodes(G, pos=pos, nodelist=[tf], node_color=color, node_size=100)
            nx.draw_networkx_edges(G, pos=pos, edgelist=list(G.edges(tf)), edge_color=color, width=1, alpha=0.3)


plt.title('TF interaction network')
        
save_path = os.path.join(save_folder, 'network.png')
plt.savefig(save_path)
        
plt.close()
