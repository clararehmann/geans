import pandas as pd,  numpy as np
import os, pickle
from goatools import obo_parser

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
        compIDs = np.nan if isinstance(ag.iloc[i]['Computed GO Component IDs'], float) else ag.iloc[i]['Computed GO Component IDs'].split(';')
        funcIDs = np.nan if isinstance(ag.iloc[i]['Computed GO Function IDs'], float) else ag.iloc[i]['Computed GO Function IDs'].split(';')
        procIDs = np.nan if isinstance(ag.iloc[i]['Computed GO Process IDs'], float) else ag.iloc[i]['Computed GO Process IDs'].split(';')

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


def main():
    # Load GO OBO file
    goDAG = obo_parser.GODag(gopath)
    # Load gene statistics data
    ag = pd.read_csv(statpath, sep='\t')
    print(ag.head())
    goDICT = build_go_dict(goDAG, ag)
    #with open(f'{savepath}_GOterms_stats.pkl', 'wb') as f:
    #    pickle.dump(goDICT, f)
    #    print(f'Saved GO terms statistics to {savepath}_GOterms_stats.pkl')
    
    print(goDICT['Component'])  # Print the keys of the GO dictionary to verify structure


    # create array of gene statistics
    for goTYPE in goDICT.keys():
        # Collect all data for this GO type
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
        print(df.head())
        
        # Add level and depth for each GO term
        df['Level'] = df['GO_ID'].apply(lambda x: goDAG[x].level if x in goDAG else None)
        df['Depth'] = df['GO_ID'].apply(lambda x: goDAG[x].depth if x in goDAG else None)
        
        # Save to file
        df.to_csv(f'{savepath}_{goTYPE}_terms_stats.csv', sep='\t', index=False)
        print(f'Saved {goTYPE} terms statistics to {savepath}_{goTYPE}_terms_stats.csv')


if __name__ == "__main__":
    main()