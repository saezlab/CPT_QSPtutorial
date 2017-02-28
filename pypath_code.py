# node1 and node2 are gene symbols (HUGO)  
# pa is the network loaded from databases with pa.init_network()
# get all one-way directed interactions between two nodes, defined by "=+=>", "=-=>", "<=+=" or "<=-="
# add also the number of references for all the interactions, and the number of sources for each interacion
# output: listinteractions=[(edge,(source,arrow,target,nbsources,nbrefs)),...] 
def list_interactions(pa,node1,node2):
    if pa.gs(node1)==None:
        print node1+" not found in Omnipath"
        return(None)
    if pa.gs(node2)==None:
        print node2+" not found in Omnipath"
        return(None)
    edge = pa.graph.es.select(_between = ([pa.gs(node1).index], [pa.gs(node2).index])) # select interaction if it exists
    listinteractions=list()
    if len(edge)>0:
        lines = str(edge[0]['dirs']).split("\n") # Omnipath returns a text for each interaction with additional information
	# select lines from the text giving details on signed interactions:
	signedinteractions = [l for l in lines if "=+=>" in l or "<=+=" in l or "=-=>" in l or "<=-=" in l]
        interactions = [line.split("::")[0] for line in signedinteractions] # signed interactions
        databases = [line.split("::")[1] for line in signedinteractions] # associated databases
	# in the texts for signed interactions replace gene names (Omnipath ID) with gene symbols (HUGO):
        newinteractions=[el.replace(pa.gs(node1)['name'],node1).replace(pa.gs(node2)['name'],node2).replace("\t","") for el in interactions]
        for i,e in enumerate(newinteractions): # reverse interactions where target is on the left
            splitint = [f for f in e.split(" ") if f!=""]
            if splitint[1]=="<=+=":
                splitint=[splitint[2],"=+=>",splitint[0]]
            elif splitint[1]=="<=-=":
                splitint=[splitint[2],"=-=>",splitint[0]]
            splitint.append(str(len(edge['references'][0])))
            splitint.append(str(databases[i]))
            listinteractions.append(tuple(splitint))

    listinteractions = ["\t".join(list(e)) for e in listinteractions]
    print("\n".join(listinteractions))
    return(edge,listinteractions)

    

# plot all paths between node and some nodes in the list with the shortest length
# listofnodes is a list of gene symbols (HUGO)
def plot_paths_to_listofnodes(pa,node,listofnodes,size=500):
    if pa.gs(node)==None:
        print(node+" was not found")
        return None
    nodesindices = [pa.gs(n).index for n in listofnodes if pa.gs(n) is not None]
    
    # find the length of the shortest path between node and the list
    dist = pa.shortest_path_dist(pa.graph,([pa.gs(node).index],nodesindices))
    mindist = min(dist)-1
    print("Distance: "+str(mindist))
    # find all shortest paths between node and the list
    paths = pa.find_all_paths(pa.gs(node).index,nodesindices, mode = 'ALL', maxlen = mindist)
    
    # select paths composed of signed interactions -> list signedpaths
    signedpaths = list()
    for p in paths:
        keeppath=True
        for i in range(len(p)-1):
            signedinter = list_interactions(pa,pa.graph.vs[p[i]]['label'],pa.graph.vs[p[i+1]]['label'])[1]
            if len(signedinter)==0: # unsigned interaction
                keeppath=False
                break
        if keeppath==True:
            signedpaths.append(p)
    setgenes = set([item for sublist in signedpaths for item in sublist])
    
    print('Number of nodes in subgraph: {}'.format(len(setgenes)))
    sgraph = pa.graph.induced_subgraph(setgenes)
    
    # turn undirected subgraph to directed subgraph	
    dgenes = set([])
    for i,gene in enumerate([v['label'] for v in sgraph.vs]):
        if pa.dgs(gene)!=None:
            dgenes.add(pa.dgs(gene))
        else:
            print(gene+" not found in dgs") # dgs=directed graph
    dgraph = pa.dgraph.induced_subgraph(dgenes)    
    
    # prepare color and size vectors for nodes
    vcolors=list()
    vsizes=list()
    for v in dgraph.vs:
        if v['label'] == node:
            vcolors = vcolors+['#A9D0F5'] # main node in blue
            vsizes =  vsizes+[40]
        elif v['label'] in nodesinmodel: 
            vcolors = vcolors+['#CEF6CE'] # nodes in the list in green
            vsizes =  vsizes+[30]
        else:
            vcolors = vcolors+['#F79F81'] # orange for additional nodes
            vsizes =  vsizes+[20]                 

    # prepare color vector for edges    
    ecolors=list()
    for edge in dgraph.es:
        neg = max([len(c) for c in edge['dirs'].negative_sources.values()])
        pos = max([len(c) for c in edge['dirs'].positive_sources.values()])
        if neg>0 and pos>0:
           ecolors.append("blue") # dual edge
        elif neg>0:
           ecolors.append("red") # inhibition edge: red
        elif pos>0:
           ecolors.append("green") # activation edge      
        else:
           ecolors.append("black") # directed but unsigned edge

    # prepare width vector for edges    
    coef=0.2 # width = nrefs*0.2     
    ewidth=list()
    for edge in dgraph.es:
        if isinstance(edge['dirs'],list) or edge['dirs']==None: # no reference
            ewidth.append(coef)
        else:
            nref=len(edge['references'])
            ewidth.append(min(nref*coef,10) )
    
    # plot directed graph
    layout= dgraph.layout_fruchterman_reingold(repulserad = dgraph.vcount() ** 2.8, maxiter = 1000, area = dgraph.vcount() ** 2.3)    
    igraph.plot(dgraph, node+'_to_list_'+str(mindist)+'_directed.png', 
               vertex_color = vcolors, vertex_size=vsizes, vertex_frame_width = 0, vertex_label_size = 15, 
               edge_color = ecolor, edge_width=ewidth, layout=layout, **visual_style)   

    # write list of formatted interactions in a file
    listinteractions = list()                 
    for pair in get_directed_edges_labels(pa.graph,sgraph.es):
        if len(pair)>0:
            listinteractions.extend(list_interactions(pa,pair[0],pair[1])[1])
    listinteractions = ["\t".join(list(e)) for e in listinteractions]
    fo=open(node+'_to_list_dist_'+str(mindist)+'_alldirectedinteractions.txt','w+')
    fo.write("\n".join(listinteractions))
    fo.close()
    
    # Print paths
    print('Number of paths: {}'.format(len(signedpaths)))
    gspaths = [[pa.graph.vs[i]['label'] for i in p] for p in signedpaths]
    print(gspaths)

    return(sgraph)      
    
    

