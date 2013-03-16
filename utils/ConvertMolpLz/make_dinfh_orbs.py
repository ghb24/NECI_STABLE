import networkx as nx
import numpy


def match_nodes(node1, node2):
    return True
def match_edges(edge1, edge2):
    if (edge1['weight'] == edge2['weight']):
        return True;
    else:
        return False;
def match_edges2(edge1, edge2):
    if (abs(edge1['weight']) == abs(edge2['weight'])):
        return True;
    else:
        return False;



G = nx.Graph()


f = open("qcdmrg.int1","r")
lines = f.readlines()
norbs = int(lines[0].split()[0])

int1 = numpy.zeros(shape=(norbs, norbs));
G.add_nodes_from(range(norbs))

for line in lines[1:]:
    tokens = line.split()
    i = int(tokens[0])
    j = int(tokens[1])
    int1[int(tokens[0]),int(tokens[1])] = float(tokens[2])
    int1[int(tokens[1]),int(tokens[0])] = float(tokens[2])

    if (abs(int1[i,j]) > 1e-8):
        G.add_weighted_edges_from([(i,j, int(int1[i,j]*1e6))])

subgraph_indices= nx.connected_components(G)

disjointgraphs = []

for indices in subgraph_indices :
    disjointgraphs.append(G.subgraph(indices))

print "DISJOINTED GRAPHS ***********"
for graphs in disjointgraphs:
    print graphs.nodes(), "  ->  ", len(graphs.nodes())
print "**********\n"

index_maps=[]
usedgraphs=[]
GM=nx.algorithms.isomorphism.GraphMatcher

negatives=[]



screen_error = 1e-9
index_maps=[]
usedgraphs=[]
index = 0
for i in range(len(disjointgraphs)):
    if (i in usedgraphs):
        continue
    usedgraphs.append(i)
    g1 = disjointgraphs[i]
    for j in range(i,len(disjointgraphs)):
        if (j in usedgraphs):
            continue
        g2 = disjointgraphs[j]
        GM=nx.algorithms.isomorphism.GraphMatcher(g1, g2, match_nodes, match_edges)
        if (GM.is_isomorphic()):
            break
    if (GM.is_isomorphic()):
        usedgraphs.append(j)
        index_maps.append([[],[]])
        #print GM.mapping
        for imap in GM.mapping:
            index_maps[index][0].append(imap)
            index_maps[index][1].append(GM.mapping[imap])
    else:
        index_maps.append(g1.nodes())

    index+=1

for imaps in index_maps:
    print imaps



screen_error = 1e-9
index_maps=[]
usedgraphs=[]
index = 0
for i in range(len(disjointgraphs)):
    if (i in usedgraphs):
        continue
    usedgraphs.append(i)
    g1 = disjointgraphs[i]
    for j in range(i,len(disjointgraphs)):
        if (j in usedgraphs):
            continue
        g2 = disjointgraphs[j]
        GM=nx.algorithms.isomorphism.GraphMatcher(g1, g2, match_nodes, match_edges2)
        if (GM.is_isomorphic()):
            break

    if (GM.is_isomorphic()):
        usedgraphs.append(j)
        index_maps.append([[],[]])
        #print GM.mapping
        firsta = -1
        firstb = -1
        for imap in GM.mapping:
            index_maps[index][0].append(imap)
            index_maps[index][1].append(GM.mapping[imap])
            if (firsta == -1):
                firsta = imap
                firstb = GM.mapping[imap]
            elif ( abs(int1[firsta,imap]) > screen_error and abs(int1[firsta, imap]+int1[firstb,GM.mapping[imap]]) < screen_error*0.5) :
                negatives.append(GM.mapping[imap])
            #print int1[firsta, imap], firsta, imap, firstb, GM.mapping[imap], int1[firstb,GM.mapping[imap]]
    else:
        index_maps.append(g1.nodes())

    index+=1

#for imaps in index_maps:
#    print imaps

print negatives

print "Check if negatives work"
index = 0
for i in range(len(disjointgraphs)):
    for edges in disjointgraphs[i].edges(data=True):
        scale = 1.0
        if edges[0] in negatives:
            scale*= -1.0
        if edges[1] in negatives:
            scale*=-1.0
        edges[2]['weight'] =edges[2]['weight']*scale

for i in range(len(disjointgraphs)):
    if (i in usedgraphs):
        continue
    usedgraphs.append(i)
    g1 = disjointgraphs[i]
    for j in range(i,len(disjointgraphs)):
        if (j in usedgraphs):
            continue
        g2 = disjointgraphs[j]
        GM=nx.algorithms.isomorphism.GraphMatcher(g1, g2, match_nodes, match_edges)
        if (GM.is_isomorphic()):
            break

    if (GM.is_isomorphic()):
        usedgraphs.append(j)
        index_maps.append([[],[]])
        for imap in GM.mapping:
            index_maps[index][0].append(imap)
            index_maps[index][1].append(GM.mapping[imap])
    else:
        index_maps.append(g1.nodes())

    index+=1

for imaps in index_maps:
    print imaps,"  ->  ", len(imaps)


print "**************************"


