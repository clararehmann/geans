import pandas as pd,  numpy as np
import os, pickle
from goatools import obo_parser
from goatools.godag.go_tasks import get_go2parents, get_go2children

gopath='/home/crehmann/vectorcomp/data/go-basic.obo' # path to GO OBO file
statpath='/home/crehmann/vectorcomp/out/gamb_colu_arab_gene_stats_wide_GO_IDs.txt' # path to gene stats file
savepath='/home/crehmann/vectorcomp/out/gamb_colu_arab'


def updatedict(goID, goDAG, goDICT, agSLICE):
    transcript_ID = agSLICE['Transcript']
    gene_ID = agSLICE['Gene']
    rec = goDAG[goID]
    parents = rec.get_all_parents()
    parents = {goDAG[x].id: {'level': goDAG[x].level,
                              'depth': goDAG[x].depth} for x in parents}
    children = rec.get_all_children()
    children = {goDAG[x].id: {'level': goDAG[x].level,
                              'depth': goDAG[x].depth} for x in children}
    if goDICT.get(goID) is None: # if component ID not in dictionary
        goDICT[goID] = {'gene': {gene_ID: {transcript_ID: dict(agSLICE.filter(regex='^(gene_|aa_|cds_|ns_|ss_|div_)'))}}, 
                            'parents': parents, 
                            'children': children} # add to dictionary
    else:
        if goDICT[goID]['gene'].get(gene_ID) is None: # if gene not in component ID
            goDICT[goID]['gene'][gene_ID] = {transcript_ID: dict(agSLICE.filter(regex='^(gene_|aa_|cds_|ns_|ss_|div_)'))} # add gene and transcript
        else:
            goDICT[goID]['gene'][gene_ID][transcript_ID] = dict(agSLICE.filter(regex='^(gene_|aa_|cds_|ns_|ss_|div_)'))
    return goDICT

def build_go_dict(goDAG, ag):
    compdict = {}
    funcdict = {}
    procdict = {}

    ivals = range(len(ag))
    for i in ivals: 
        compIDs = np.nan if isinstance(ag.iloc[i]['Curated GO Component IDs'], float) else ag.iloc[i]['Curated GO Component IDs'].split(';')
        funcIDs = np.nan if isinstance(ag.iloc[i]['Curated GO Function IDs'], float) else ag.iloc[i]['Curated GO Function IDs'].split(';')
        procIDs = np.nan if isinstance(ag.iloc[i]['Curated GO Process IDs'], float) else ag.iloc[i]['Curated GO Process IDs'].split(';')

        if not isinstance(compIDs, float):
            jvals = range(len(compIDs))
            for j in jvals:
                try:
                    cid = compIDs[j] # get first component ID
                    #print(cid)
                    compdict = updatedict(cid, goDAG, compdict, ag.iloc[i])
                except KeyError:
                    continue
            
        if not isinstance(funcIDs, float):
            jvals = range(len(funcIDs))
            for j in jvals:
                try:
                    cid = funcIDs[j] # get first function ID
                    funcdict = updatedict(cid, goDAG, funcdict, ag.iloc[i])
                except KeyError:
                    continue
        if not isinstance(procIDs, float):
            jvals = range(len(procIDs))
            for j in jvals:
                try:
                    cid = procIDs[j] # get first process ID
                    procdict = updatedict(cid, goDAG, procdict, ag.iloc[i])
                except KeyError:
                    continue
    print(compdict)
    GOdict = {'Component': compdict, 'Function': funcdict, 'Process': procdict}
    return GOdict

def get_stats(goDICT, goTERM):
    """
    Get gene statistics for a specific GO term.
    """
    if goDICT.get(goTERM) is None:
        print(f'GO term {goTERM} not found in dictionary.')
        return None
    else:
        gene_stats = []
        for gene, transcripts in goDICT[goTERM]['gene'].items():
            for transcript, stats in transcripts.items():
                gene_stats.append(list(stats.values()))
        gene_stats_array = np.array(gene_stats)
        if gene_stats_array.size == 0:
            print(f'No gene statistics found for GO term {goTERM}.')
            return None
        else:
            #gene_stats_mean = np.nanmean(gene_stats_array, axis=0)
            return gene_stats_array

def getparents(gterm, go2parents_isa):
    return go2parents_isa.get(gterm, [])

def getchildren(gterm, go2children_isa):
    return go2children_isa.get(gterm, [])

def stackparents(gterm, goparents):
    parents = getparents(gterm, goparents)
    stacked = []
    stacked.extend(f"{p};{gterm}" for p in parents)
    while parents:
        stacked.extend(f"{p};{parent}" for parent in parents for p in getparents(parent, goparents))
        parents = [p for parent in parents for p in getparents(parent, goparents)]
    return stacked

def stackchildren(gterm, gochildren):
    children = getchildren(gterm, gochildren)
    stacked = []
    stacked.extend(f"{gterm};{c}" for c in children)
    while children:
        stacked.extend(f"{c};{child}" for c in children for child in getchildren(c, gochildren))
        children = [child for c in children for child in getchildren(c, gochildren)]
    return stacked

def go_edge_matrix(goDICT, goDAG, genedf, genedfcol, goparents, gochildren):
    """
    Create a GO edge matrix from the GO dictionary and gene DataFrame.
    Rows represent (parent, child) relationships present for a given gene,
    columns represent GO terms.

    Parameters:
    goDICT (dict): Dictionary containing GO terms and their relationships.
    goDAG (GODag): GO Directed Acyclic Graph.
    genedf (DataFrame): DataFrame containing gene information.
    genedfcol (str): Column name in genedf containing GO terms.
    """
    GENES = genedf['Gene'].unique()
    # for each gene, construct a list of (parent, child) relationships both up and down the GO hierarchy
    edgedict = {}
    for gene in GENES:
        try: gene_go_terms = genedf[genedf['Gene'] == gene][genedfcol].dropna().unique()[0].split(';')
        except IndexError: continue
        for go_term in gene_go_terms:
            if go_term in goDICT:
                relationships = stackparents(go_term, goparents) + stackchildren(go_term, gochildren)
                if edgedict.get(gene) is None:
                    edgedict[gene] = relationships
                else:
                    edgedict[gene] = edgedict.get(gene, []) + relationships
    # create a list of unique edges
    edges = np.unique(sum([i for i in edgedict.values()], []))
    # initialize an edge table with zeros
    edgetable = np.zeros((len(edges), len(edgedict.keys())), dtype=int)
    # fill the edge table with ones where edges are present
    for gene, relationships in edgedict.items():
        for edge in relationships:
            if edge in edges:
                edgetable[edges.tolist().index(edge), list(edgedict.keys()).index(gene)] = 1
    return edges, edgetable, edgedict


def main():
    # Load GO OBO file
    goDAG = obo_parser.GODag(gopath)
    # set up relationships for building GO DAG
    # optional_relationships = {'is_a', 'part_of', 'regulates', 'posit...
    optional_relationships = set()
    go2parents_isa = get_go2parents(goDAG, optional_relationships)
    go2children_isa = get_go2children(goDAG, optional_relationships)
    # Load gene statistics data
    ag = pd.read_csv(statpath, sep='\t')
    print(ag.head())
    goDICT = build_go_dict(goDAG, ag)
    #with open(f'{savepath}_GOterms_stats.pkl', 'wb') as f:
    #    pickle.dump(goDICT, f)
    #    print(f'Saved GO terms statistics to {savepath}_GOterms_stats.pkl')
    
    #print(goDICT['Component'])  # Print the keys of the GO dictionary to verify structure


    # create array of gene statistics
    for goTYPE in goDICT.keys():
        # Collect all data for this GO type
        """
        all_data = []
        
        # Extract values from the GO dictionary
        go_ids = list(goDICT.get(goTYPE, {}).keys())
        
        for go_id in go_ids:
            for gene_id, transcripts in goDICT[goTYPE][go_id]['gene'].items():
                for transcript_id, stats in transcripts.items():
                    # Create a row with GO_ID, Gene, Transcript, and all stats
                    row_data = {
                        'GO_ID': go_id,
                        'Gene': gene_id,
                        'Transcript': transcript_id,
                        **stats  # Unpack the statistics dictionary
                    }
                    all_data.append(row_data)

        # Convert to DataFrame
        print(f'Converting statistics to DataFrame for {goTYPE} terms...')
        df = pd.DataFrame(all_data)
        print(f'DataFrame shape: {df.shape}')
        
        # Add level and depth for each GO term
        df['Level'] = df['GO_ID'].apply(lambda x: goDAG[x].level if x in goDAG else None)
        df['Depth'] = df['GO_ID'].apply(lambda x: goDAG[x].depth if x in goDAG else None)
        
        # Save to file
        df.to_csv(f'{savepath}_{goTYPE}_terms_stats.txt', sep='\t', index=False)
        print(f'Saved {goTYPE} terms statistics to {savepath}_{goTYPE}_terms_stats.txt')
        """

        goterms = goDICT.get(goTYPE, {})
        print(f'Creating GO edge matrix for {goTYPE} terms...')
        edges, edgetable, edgedict = go_edge_matrix(goterms, goDAG, ag, f'Curated GO {goTYPE.capitalize()} IDs', go2parents_isa, go2children_isa)
        df = pd.DataFrame(edgetable, index=edges, columns=edgedict.keys())
        df.to_csv(f'{savepath}_{goTYPE}_edge_matrix.txt', sep='\t')
        print(f'Saved {goTYPE} edge matrix to {savepath}_{goTYPE}_edge_matrix.txt')

if __name__ == "__main__":
    main()