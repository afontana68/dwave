import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

j1=1
j2=2
adjacency_matrix = np.array([[ 0,j1, 0, 0,j1,j2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                             [j1, 0,j1, 0,j2,j1,j2, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                             [ 0,j1, 0,j1, 0,j2,j1,j2, 0, 0, 0, 0, 0, 0, 0, 0], 
                             [ 0, 0,j1, 0, 0, 0,j2,j1, 0, 0, 0, 0, 0, 0, 0, 0], 
                             [j1,j2, 0, 0, 0,j1, 0, 0,j1,j2, 0, 0, 0, 0, 0, 0], 
                             [j2,j1,j2, 0,j1, 0,j1, 0,j2,j1,j2, 0, 0, 0, 0, 0], 
                             [ 0,j2,j1,j2, 0,j1, 0,j1, 0,j2,j1,j2, 0, 0, 0, 0], 
                             [ 0, 0,j2,j1, 0, 0,j1, 0,j2,j1,j2,j1, 0, 0, 0, 0], 
                             [ 0, 0, 0, 0,j1,j2, 0,j2, 0,j1, 0, 0,j1,j2, 0, 0], 
                             [ 0, 0, 0, 0,j2,j1,j2,j1,j1, 0,j1, 0,j2,j1,j2, 0], 
                             [ 0, 0, 0, 0, 0,j2,j1,j2, 0,j1, 0,j1, 0,j2,j1,j2], 
                             [ 0, 0, 0, 0, 0, 0,j2,j1, 0, 0,j1, 0, 0, 0,j2,j1], 
                             [ 0, 0, 0, 0, 0, 0, 0, 0,j1,j2, 0, 0, 0,j1, 0, 0], 
                             [ 0, 0, 0, 0, 0, 0, 0, 0,j2,j1,j2, 0,j1, 0,j1, 0], 
                             [ 0, 0, 0, 0, 0, 0, 0, 0, 0,j2,j1,j2, 0,j1, 0,j1],
                             [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,j2,j1, 0, 0,j1, 0]])
print(adjacency_matrix)

G = nx.Graph(adjacency_matrix, nodetype=int)
G.edges()

plt.axis('off')
#pos = nx.random_layout(G)
pos = np.array([[-1, 0.5], [-0.5, 0.5], [0, 0.5], [0.5,0.5], 
                [-1, 0], [-0.5, 0], [0, 0], [0.5,0],
                [-1, -0.5],[-0.5, -0.5],[0, -0.5], [0.5,-0.5],
                [-1, -1],[-0.5, -1],[0, -1], [0.5,-1]])
#print(pos)
nx.draw_networkx(G,pos)
plt.show()

from dwave.system import DWaveSampler, EmbeddingComposite, DWaveCliqueSampler
import dimod
import dwave.inspector

bqm = dimod.from_networkx_graph(G, vartype='BINARY', edge_attribute_name='weight', node_attribute_name='weight')

sampleset = EmbeddingComposite(DWaveSampler()).sample(bqm, num_reads=10000, label='cfrustrated')

sampleset.variables
sampleset.record

res = np.zeros(len(sampleset.record)) 
s = 0

for k, rec in enumerate(sampleset.record):
    states = rec[0]
    print(states)
    for i, sol in enumerate(states):
        p = int(sampleset.variables[i])

r=sampleset.record.num_occurrences
idx=np.argmax(r)
print(idx,r)
solution=sampleset.record
#print(solution)

plt.bar(range(len(solution)),r)

dwave.inspector.show(sampleset)
