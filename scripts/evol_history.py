import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
import re
import os
import json
from tqdm import tqdm
from networkx.drawing.nx_agraph import graphviz_layout


from outbreak_data import outbreak_data
from outbreak_tools import outbreak_tools

""" Given a cryptic mutation cluster, plot the stepwise evolution of the mutations leading up to it."""

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def parse_query_list(query_list):
    query_list = query_list.strip('[]').split(', ')
    query_list = [x.strip().strip("'").strip('"') for x in query_list]
    return query_list


def get_aa_site(mut):
    if 'DEL' in mut:
        return int(mut.split('DEL')[1].split('/')[0])
    else: 
        return int(mut.split(':')[1][1:-1]
)


def get_duplicate_indices(lst):
    duplicates = {}
    for i, x in enumerate(lst):
        if x in duplicates:
            duplicates[x].append(i)
        else:
            duplicates[x] = [i]
    duplicates = [v for k,v in zip(duplicates.keys(),duplicates.values()) if len(v)>1]
    duplicates_flat = [d for d1 in duplicates for d in d1]
    return duplicates_flat

def main():
    df = pd.read_csv('cryptic_variants_meta.tsv',sep='\t')
    os.system('rm evoplots/*')

    df['Covariants'] = df['query'].apply(parse_query_list)
    df['Covariants'] = df['Covariants'].apply(lambda x: tuple(sorted(x,key=get_aa_site)))

    df['num_muts'] = df['query'].apply(lambda x:len(x))
    df = df[df['num_muts']>=2]

    df.drop_duplicates(subset=['Covariants', 'collection_date', 'location'],keep='first',inplace=True)

    df['coverage_start'] = df['coverage_start'].apply(lambda x: (x - 21563) // 3 ) # TODO replace positions with actual genome positions
    df['coverage_end'] = df['coverage_end'].apply(lambda x: (x - 21563) // 3 )

    # aggregate metadata for each unique cluster
    df_aggregate = df.groupby('Covariants').agg({'count':tuple,'location':tuple,'collection_date':tuple,
                                        'coverage_start':tuple,'coverage_end':tuple,
                                        'num_clinical_detections':'mean','num_muts':'min'}).sort_values(by='num_muts',ascending=False)
    
    df_aggregate['Total_Obs']= df_aggregate['count'].apply(lambda x: len(x))

    df_aggregate.to_csv('df_aggregate.tsv',sep='\t')

    # pull only cryptics (here, defined as 1 or less clinical detections)
    df_cryptics = df_aggregate[df_aggregate['num_clinical_detections']<=1]
    df_cryptics['Covariants'] = df_cryptics.index
    df_cryptics = df_cryptics[df_cryptics['Total_Obs']>=2] # At least 2 detections

    # remove cryptics that are subsets of other cryptics. 
    j=0
    while j<df_cryptics.shape[0]:
        test_clust = set(df_cryptics.index[j])
        df_subset = df_cryptics[df_cryptics['Covariants'].apply(lambda x: set(x).issubset(test_clust) & (set(x)!=test_clust))]
        df_cryptics = df_cryptics.drop(index=df_subset.index)
        j += 1

    ## now generate subtrees leading up to these cryptics. 
    for i0 in tqdm(range(0,df_cryptics.shape[0])):
        test_clust = set(df_cryptics.index[i0])
        # Compare with clusters spanning the same region as test_clust
        df0 = df.copy()
        positions = np.array([get_aa_site(t) for t in test_clust])
        df0 = df0[df0['coverage_start']<=positions.min()]
        df0 = df0[df0['coverage_end']>=positions.max()]

        # now group
        df0 = df0.groupby('Covariants').agg({'count':tuple,'location':tuple,'collection_date':tuple,
                                                'coverage_start':tuple,'coverage_end':tuple,
                                                'num_clinical_detections':'mean',
                                                'num_muts':'min'})
                                                
        df0['Total_Obs']= df0['count'].apply(lambda x: len(x))
        df0.to_csv("df0.tsv",sep='\t')
        df0 = df0[df0['Total_Obs']>=2]
        df0['Covariants'] = df0.index

        df0['num_muts'] = df0['Covariants'].apply(lambda x:len(x))
        df0 = df0.sort_values(by='num_muts',ascending=False)

        df_superset = df0[df0['Covariants'].apply(lambda x: set(x).issubset(test_clust))]

        if df_superset.shape[0]<=1:
            print(f'No detected stepwise evolution for {test_clust}')
            continue

        df_superset.loc[:,'Covariants'] = df_superset['Covariants'].apply(lambda x: tuple(sorted(list(x),key=get_aa_site)))


        parent_list = []
        new_muts_list = []
        for j,c in enumerate(df_superset['Covariants']):
            # get subsets
            subsets = df_superset[df_superset['Covariants'].apply(lambda x: set(x).issubset(c) & (x!=c))]
            if subsets.shape[0]>0:
                subsets = pd.concat((subsets,pd.Series(subsets['Covariants'].apply(lambda x: len(set(c)-set(x))),name='difference')),axis=1)
                minDiff = subsets['difference'].min()
                subsets = subsets[subsets['difference']==minDiff]
                
                if subsets.shape[0] ==1:
                    parent_list.append([subsets['Covariants'].iloc[0]])
                    new_muts_list.append([tuple(set(c)-set(subsets['Covariants'].iloc[0]))])
                else:
                    parent_list.append(list(subsets['Covariants']))
                    new_muts_list.append([tuple(set(c)- set(t0)) for t0 in subsets['Covariants']])
                # print(c,subsets)
            else:
                parent_list.append(None)
                new_muts_list.append(c)


        # for each, define parent, unless smallest cluster. 
        edge_labels = {}
        G = nx.DiGraph()
        # assemble graph
        # make edges to first level
        for j,c in enumerate(df_superset['Covariants']):
            if parent_list[j] is not None:
                for l in range(0,len(parent_list[j])):
                    G.add_edge(parent_list[j][l],c,weight=1)
                    edge_labels[(parent_list[j][l],c)] = ','.join(new_muts_list[j][l])        

        try:
            if len(nx.dag_longest_path(G))<=2:
                continue
        except:
            pass
        fig,ax = plt.subplots()#figsize=(5,1.5))
        pos=graphviz_layout(G, prog='dot',args="-Grankdir='LR' -Goverlap=false",root=0)
        nx.draw_networkx_nodes(G, pos, node_color="cornflowerblue",alpha=0.9,margins=0.1,node_size=80)
        nx.draw_networkx_edges(
            G,
            pos,
            width=1,
            alpha=0.8,
            arrows=True,
            edge_color="grey",
        )
        labels = {c:o for c,o in zip(df_superset['Covariants'],list(df_superset['Total_Obs']))}
        if len(labels)> len(pos):
            labels = {key:labels[key] for key in pos.keys()}
        nx.draw_networkx_labels(G, pos = pos, labels = labels, font_color="white", font_size=6)#,verticalalignment='bottom')
        nx.draw_networkx_edge_labels(G, pos = pos, edge_labels=edge_labels, font_color='black',font_size=5)
        minX = 500
        maxX = 900

        ymin, ymax = ax.get_ylim()
        y_length = ymax - ymin
        
        xmin, xmax = ax.get_xlim()
        x_length = xmax - xmin

        start = df_superset['Covariants'].iloc[-1]
        minXpos = np.min([x[0] for x in pos.values()])

        for pk in pos.keys():
            if G.in_degree(pk)==0:
                ax.text(pos[pk][0]-x_length*0.06,pos[pk][1],'\n'.join(pk),fontsize=5,verticalalignment='center',horizontalalignment='right')
            if G.out_degree(pk)==0:
                ax.text(pos[pk][0]+x_length*0.06,pos[pk][1],'\n'.join(pk),fontsize=5,verticalalignment='center',horizontalalignment='left')
        # add number of clinical detections below ww detection count
        clin_detects = {c:counts for c,counts in zip(df_superset['Covariants'],df_superset['num_clinical_detections'])}
        for pk in pos.keys():
            ax.text(pos[pk][0],pos[pk][1]-y_length*0.03,int(df0.loc[[pk],'num_clinical_detections']),
                    fontsize=5,verticalalignment='center',horizontalalignment='center',color='red')
        # add the sites the variant was detected at. 
        for pk in pos.keys():
            locs = set(list(df0.loc[[pk],'location'][0]))
            for j,l0 in enumerate(['South Bay','Point Loma','Encina']):
                if l0 in locs:
                    ax.scatter(pos[pk][0]+x_length*0.02,pos[pk][1]+ y_length*0.02*(1-j),marker='s',s=5,color='darkblue',linewidth=0.5)
                else:
                    ax.scatter(pos[pk][0]+x_length*0.02,pos[pk][1]+ y_length*0.02*(1-j),marker='s',s=5,edgecolors='darkblue',facecolors='none',linewidth=0.5)
        
        for pk in pos.keys():
            tp = pd.to_datetime(pd.Series(df0.loc[[pk],'collection_date'][0],name='times'))
            if len(tp)==1:
                    ax.text(pos[pk][0],pos[pk][1]+y_length*0.05,
                    f"{tp.iloc[0].strftime('%Y-%m-%d')}",
                    fontsize=4,verticalalignment='center',horizontalalignment='center',color='black')
            else:
                ax.text(pos[pk][0],pos[pk][1]+y_length*0.05,
                        f"{tp.min().strftime('%Y-%m-%d')}\n-{tp.max().strftime('%Y-%m-%d')}",
                        fontsize=4,verticalalignment='center',horizontalalignment='center',color='black')

        plt.gca().set_frame_on(False)

        plt.savefig(f'evoplots/ww_evo_seq{i0}.pdf',transparent=True)
        plt.close('all')

if __name__ == "__main__":
    main()
    