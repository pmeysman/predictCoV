import os
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

def read_txt_files(folder_path):
    temp_frames = {}
    data_frames = {}
        
    # Check if the folder exists
    if not os.path.exists(folder_path):
        print(f"Folder '{folder_path}' not found.")
        return
    
    # Loop through all files in the folder with the .txt extension
    for filename in os.listdir(folder_path): 
        if filename.endswith(".txt"):
            file_path = os.path.join(folder_path, filename)

            # Read the tab-delimited file into a pandas DataFrame
            try:
                temp_frames[filename] = pd.read_csv(file_path, sep='\t')
                
                #Normalize clonecounts
                temp_frames[filename]['cloneFraction'] = temp_frames[filename]['cloneFraction'] / np.sum(temp_frames[filename]['cloneFraction'])

                temp_frames[filename]['TCR'] = temp_frames[filename]['v_call'] + '_' + temp_frames[filename]['junction_aa']


            except Exception as e:
                print(f"Error reading file '{filename}': {e}")

    for filename in temp_frames:
    # Extract sample ID from the filename
        sample_id = filename[:-4]
        print(sample_id)

        data_frames[sample_id] = temp_frames[filename]

    return data_frames

# Specify the folder path
folder_path = '../annotated'

# Read the files and get the dictionary of DataFrames
result_dict = read_txt_files(folder_path)

print(str(len(result_dict)))

database_location="../../../../imdb/"
epirepofile = database_location+"EpiRepo.txt"

epirepo = pd.read_csv(epirepofile,sep="\t")

sarscov2epitopes = list(epirepo.loc[epirepo['origin'] == "SARS-CoV-2"]["epitope"])
#print(sarscov2epitopes)

spike = list(epirepo.loc[(epirepo['origin'] == "SARS-CoV-2") & (epirepo['protein'] == "Spike/surface glycoprotein (S)")]["epitope"])
#print(spike)

cov_results = dict()
for sample,df in result_dict.items():
    cov_results[sample] = df.loc[(df["Epitope"].isin(spike)) & (df["Score"] > 0.23)]

def plot_row_counts(data_frames,filename):
    # Create a bar graph showing the number of rows for each DataFrame
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Extract TCR sequences count and file names
    sample_ids = list(data_frames.keys())
    row_counts = [len(df) for df in data_frames.values()]

    # Sort sample IDs based on CD4 and CD8
    sorted_sample_ids_cd4 = sorted([id for id in sample_ids if id.endswith('CD4')])
    sorted_sample_ids_cd8 = sorted([id for id in sample_ids if id.endswith('CD8')])
    sorted_sample_ids = sorted_sample_ids_cd4 + sorted_sample_ids_cd8

    # Extract all unique TCR chains across all samples
    #all_tcr_chains = set(chain for df in data_frames.values() for chain in df['v_call'].str[:3].unique())

    # Create a color map for TCR chains
    tcr_colors = {'TRA': 'blue', 'TRB': 'green', 'TRG': 'orange', 'TRD': 'red','IGL':'grey'}

    # Initialize bottom values for stacking
    bottom = [0] * len(sorted_sample_ids)

    # Plotting for all samples
    for sample_id in sorted_sample_ids:
        tcr_chains = set(data_frames[sample_id]['v_call'].str[:3].unique())
        for chain in tcr_colors.keys():
            count = len(data_frames[sample_id][data_frames[sample_id]['v_call'].str[:3] == chain])
            ax.bar(sample_id, count, color=tcr_colors[chain], bottom=bottom[sorted_sample_ids.index(sample_id)])
            bottom[sorted_sample_ids.index(sample_id)] += count

    # Set x-tick positions and labels
    ax.set_xticks(range(len(sorted_sample_ids)))
    ax.set_xticklabels(sorted_sample_ids, rotation=45, ha='right')

    ax.set_xlabel('Sample ID')
    ax.set_ylabel('Number of TCR Sequences')
    ax.set_title('Number of SARS-CoV-2 specific TCR Sequences in each sample')

    # Add legend for TCR chains
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=tcr_colors[tc], markersize=10, label=tc) for tc in tcr_colors.keys()]
    ax.legend(handles=legend_elements, title='TCR Chain')

    # Show the plot
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    
samples = set([x[0:5] for x in cov_results.keys()])


def enrich_tcr(sample,data_frame,t1='D0',t2='D2_28',chain='CD8'):
    s1cd = sample+t1+chain
    s2cd = sample+t2+chain
    
    df1 = data_frame[s1cd]
    df2 = data_frame[s2cd]
    
    merged_df = pd.merge(df1, df2, on='TCR', how='outer', suffixes=('_df1', '_df2'))
    
    merged_df = merged_df.fillna(0)
    
    # Define hypergeometric calculation as a function
    def calculate_hypergeom(row, df1, df2):
        M = df1['cloneCount'].sum() + df2['cloneCount'].sum()  # Total number of items in the population
        n = row['cloneCount_df2'] + row['cloneCount_df1'] +1      # Total number of successes in the population (+1 to make it more strict)
        N = df2['cloneCount'].sum()                            # Sample size
        x = row['cloneCount_df2']                              # Number of successes in the sample
        return hypergeom.sf(x - 1, M, n, N)

    # Calculate hypergeometric p-values using list comprehension
    merged_df['p_value'] = [calculate_hypergeom(row, df1, df2) for index, row in merged_df.iterrows()]
    
    reject, corrected_p_values, _, _ = multipletests(merged_df['p_value'], method='fdr_bh')
    merged_df['corrected_p_value'] = corrected_p_values
    
    # Filter significant clonotypes based on a chosen significance level (e.g., 0.05)
    significant_clonotypes = merged_df[merged_df['corrected_p_value'] < 0.05]

    # Print or use significant_clonotypes as needed
    #print(significant_clonotypes)
    
    return(merged_df)

def plot_tcrvenn(first,second,third,filename):
    #fig, ax = plt.subplots(figsize=(10, 8))
    # Extract sets from dataframes
    set_first = set(first['TCR'])
    set_second = set(second['TCR'])
    set_third = set(third['TCR'])

    # Create a Venn diagram
    venn_labels = {'100': len(set_first - set_second - set_third),
                   '010': len(set_second - set_first - set_third),
                   '001': len(set_third - set_first - set_second),
                   '110': len(set_first & set_second - set_third),
                   '101': len(set_first & set_third - set_second),
                   '011': len(set_second & set_third - set_first),
                   '111': len(set_first & set_second & set_third)}

    venn_diagram = venn3(subsets=(len(set_first), len(set_second), len(set_first & set_second),
                                  len(set_third), len(set_first & set_third), len(set_second & set_third),
                                  len(set_first & set_second & set_third)),
                         set_labels=('Primary', '1st Booster', '2nd Booster'))
    
    # Display the intersection numbers
    #for idx, label in enumerate(venn_labels):
    #    venn_diagram.get_label_by_id(idx).set_text(venn_labels[label])

    plt.savefig(filename)
    plt.close()

def plot_clonotype(row,filename):
    #fig, ax = plt.subplots(figsize=(10, 8))
    freq = [row['cloneFraction_df1_d1'],
           row['cloneFraction_df2_d1'],
           row['cloneFraction_df1_d2'],
           row['cloneFraction_df2_d2'],
           row['cloneFraction_df1'],
           row['cloneFraction_df2']]
    ticks = ['D1_0', 'D2_28', 'D3_0', 'D3_28', 'D4_0', 'D4_28']

    # Plot the line
    plt.plot(ticks, freq, marker='o')

    # Labeling the axes and the plot
    plt.xlabel('Sample')
    plt.ylabel('Clonal_frequency')
    plt.title(row['junction_aa']+' ('+str(row['Epitope_df2_d1'])+')')

    # Show the plot
    plt.savefig(filename)
    plt.close()

# Initialize an empty list to store the results
results = []
allfirst = []

for sample in samples:

    print(sample)

    #plot_row_counts({key: value for key, value in cov_results.items() if sample in key},filename='../results/'+sample+'_spike.pdf')

    for chain in ['CD8','CD4']:

        first = enrich_tcr(sample,result_dict,t1='D0',t2='D2_28',chain=chain)

        skip = ['BV137']

        results.append({
            'samplechain':sample+chain,
            'D0D2TCRs':len(first),
            'D0D2enrich':len(first[first['corrected_p_value'] < 0.05]),
            'D0D2cov':len(first[(first["Epitope_df2"].isin(spike)) & (first["Score_df2"] > 0.23)]),
            'D0D2covenrich':len(first[(first['corrected_p_value'] < 0.05) & (first["Epitope_df2"].isin(spike)) & (first["Score_df2"] > 0.23)])
            })

        covfirst = first[(first["Epitope_df2"].isin(spike)) & (first["Score_df2"] > 0.23)].copy()
        covfirst["sample"] = sample
        allfirst.append(covfirst)

        if sample in skip:
            #first.to_csv('../results/'+sample+'_enrichD0.csv')
            next
        else:
            second = enrich_tcr(sample,result_dict,t1='D3_0',t2='D3_28',chain=chain)
            third = enrich_tcr(sample,result_dict,t1='D4_0',t2='D4_28',chain=chain)

            # Venn diagram of the vaccine-expanded TCRs

            #plot_tcrvenn(first[first['corrected_p_value'] < 0.05],
            #         second[second['corrected_p_value'] < 0.05],
            #         third[third['corrected_p_value'] < 0.05],filename='../results/'+sample+'_enrichvenn_'+chain+'.pdf')

            # Same venn diagram of the vaccine-expanded TCRs but restricted to putative Spike-reactive TCRs

            #plot_tcrvenn(first[(first['corrected_p_value'] < 0.05) & (first["Epitope_df2"].isin(spike)) & (first["Score_df2"] > 0.05)],
            #     second[(second['corrected_p_value'] < 0.05) & (second["Epitope_df2"].isin(spike)) & (second["Score_df2"] > 0.05)],
            #     third[(third['corrected_p_value'] < 0.05) & (third["Epitope_df2"].isin(spike)) & (third["Score_df2"] > 0.05)],
            #     filename='../results/'+sample+'_cov_enrichvenn_'+chain+'.pdf')

            # Merging data frames across different vaccinations based on full TCR sequence to allow tracking

            #merged_fs = pd.merge(first[first['corrected_p_value'] < 0.05], second[second['corrected_p_value'] < 0.05], on='TCR', how='inner', suffixes=('_d1', '_d2'))
            #merged_all = pd.merge(merged_fs, third[third['corrected_p_value'] < 0.05], on='TCR', how='inner', suffixes=('_', '_d3'))

            #merged_all.to_csv('../results/sig/'+str(sample)+'_'+chain+'.txt',sep="\t",header=True,index=False)

            #PLot dynamics of each clonotype

            #[plot_clonotype(row,filename='../results/'+str(sample)+'_x'+str(index)+'.pdf') for index, row in merged_all.iterrows()]

            # Only merge spike-specific TCRs

            #merged_cov_fs = pd.merge(first[(first['corrected_p_value'] < 0.05) & (first["Epitope_df2"].isin(spike)) & (first["Score_df2"] > 0.23)], second, on='TCR', how='left', suffixes=('_d1', '_d2')).fillna(0)
            #merged_cov_all = pd.merge(merged_cov_fs, third, on='TCR', how='left', suffixes=('_', '_d3')).fillna(0)

            #merged_cov_all.to_csv('../results/covsig/'+str(sample)+'_'+chain+'.txt',sep="\t",header=True,index=False)


# Convert the results list to a DataFrame
results_df = pd.DataFrame(results)

# Output the results table
results_df.to_csv('../results/stats_per_sample.csv')

pd.concat(allfirst).to_csv('../results/allfirst_cov.csv')
