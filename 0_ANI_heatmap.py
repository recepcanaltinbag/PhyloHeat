import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch

import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt

genome_id_BIOMIG = 'P. sp. BIOMIG1BAC'
genome_id_BIOMIG = '6095431'


"""
distance_matrix = pd.read_csv('matrix5.csv', index_col='Strains')
genome_id_name = pd.read_csv('Closest_Genomes_v2/150_Closest_Genome_Ids/Closest_Genomes_v1.csv')
distance_matrix.to_csv('matrix5names.csv')


"""





def matrix_ANI_change_names(ani_input, ani_output):
    ANI_data = pd.read_csv(ani_input, index_col='Strains')
    distance_matrix = pd.read_csv('matrix5.csv', index_col='Strains')
    genome_id_name = pd.read_csv('Closest_Genomes_v2/150_Closest_Genome_Ids/Closest_Genomes_v1.csv')
    id_dict = dict(zip(genome_id_name['UID'], genome_id_name['Nice Name']))
    print(id_dict)
    print(ANI_data)
    ANI_data = ANI_data.rename(columns={str(key):id_dict[key] for key in id_dict}, index=id_dict)
    distance_matrix = distance_matrix.rename(columns={str(key): id_dict[key] for key in id_dict}, index=id_dict)
    print(distance_matrix)
    distance_matrix.to_csv('matrix5names.csv')
    ANI_data.to_csv(ani_output)
    return 0


def take_ANI_column(distance_matrix, ani_input, a):
    ANI_data = pd.read_csv(ani_input, index_col='Strains')
    print(distance_matrix)
    strains_list = distance_matrix['Strains'].tolist()
    genome_id_BIOMIG = '6095431'
    ANI_out_list = []
    print(strains_list)
    for strain in strains_list:
        print(strain, 'strain_name')
        try:
            print(ANI_data.at[strain, genome_id_BIOMIG], 'data')
            ANI_out_list.append((strain, ANI_data.at[strain, genome_id_BIOMIG]))
        except:
            print('error')
            ANI_out_list.append((strain, 0.95))

    print(ANI_out_list)
    df_ANI = pd.DataFrame(ANI_out_list, columns=['Strains', 'ANI'])
    print(df_ANI['ANI'].max(), df_ANI['ANI'].min())
    max_v = df_ANI['ANI'].max()
    min_v = df_ANI['ANI'].min()
    max_v = 1.00

    if a == 1:
        min_v = 0.85
    else:
        min_v = 0.95

    diff = abs(max_v - min_v)/5.0
    new_ANI_list = []
    for element in ANI_out_list:
        sum = min_v
        counter = 1
        while float(element[1]) > float(sum):
            sum = sum + diff
            counter = counter + 1
        new_ANI_list.append((element[0], counter))


    new_ANI_list = sorted(new_ANI_list,key=lambda l:l[1])
    print(new_ANI_list)
    print('len',len(new_ANI_list), len(ANI_out_list))
    df_ANI = pd.DataFrame(new_ANI_list, columns=['Strains', 'ANI'])
    lut2 = dict(zip(df_ANI['ANI'].unique(), sns.color_palette("coolwarm", 6)))
    #lut2 = dict(zip(df_ANI['ANI'].unique(), "rgb"))
    df_ANI = df_ANI.set_index('Strains')
    print(df_ANI)

    genome_id_name = pd.read_csv('Closest_Genomes_v2/150_Closest_Genome_Ids/Closest_Genomes_v1.csv')
    id_dict = dict(zip(genome_id_name['UID'], genome_id_name['Nice Name']))
    df_ANI = df_ANI.rename(index=id_dict)

    col_colors2_ANI = df_ANI['ANI'].map(lut2)
    return col_colors2_ANI

ani_input1 = 'ANI_Try_v2/output_names_ANI_v3'
ani_input2 = '16S_rRNA_Closest/ANI_values/ANIb_percentage_identity.tab'

ANI_datax = pd.read_csv(ani_input2, delimiter='\t', index_col='Strains')
ANI_datax.to_csv('16S_rRNA_Closest/ANI_values/ANIb_percentage_identity.csv')
ani_input2 = '16S_rRNA_Closest/ANI_values/ANIb_percentage_identity.csv'

ani_output1 = 'ANI_Try_v2/output_names_ANI_v4.csv'
ani_output2 = '16S_rRNA_Closest/ANI_values/ANIb_percentage_identity_names.csv'
#matrix_ANI_change_names(ani_input1, ani_output1)
#matrix_ANI_change_names(ani_input2, ani_output2)
ani_input1 = 'ANI_Try_v2/output_names_ANI_v3'


distance_matrix = pd.read_csv('matrix5.csv')
col_colors2_ANI = take_ANI_column(distance_matrix, ani_input1, 1)
col_colors3_ANI = take_ANI_column(distance_matrix, ani_input2, 2)
distance_matrix = pd.read_csv('matrix5names.csv')
distance_matrix.set_index('Strains', inplace=True)

print(col_colors2_ANI)



print(distance_matrix)


d = sch.distance.pdist(distance_matrix)
L = sch.linkage(d, method='average')
# 0.2 can be modified to retrieve more stringent or relaxed clusters, it is good for ANI
clusters = sch.fcluster(L, 0.2*d.max(), 'distance')

species_list = []
print(d,len(d))
print(0.2*d.max(), 'cluster treshold',d.max())
input()
# clusters indicices correspond to incides of original df
for i,cluster in enumerate(clusters):
    print(distance_matrix.index[i], cluster)
    species_list.append((distance_matrix.index[i], cluster))

species = pd.DataFrame(species_list, columns=['Strains', 'Cluster Groups'])
lut1 = dict(zip(species['Cluster Groups'].unique(), sns.color_palette("Paired")))

print(lut1)
species = species.set_index('Strains')
print(species)

col_colors3 = species.copy()
col_colors3.loc[col_colors3['Cluster Groups'] > 1, "Cluster Groups"] = 2
lut3 = dict(zip(col_colors3['Cluster Groups'].unique(), sns.color_palette("Set2")))

#col_colors3 = col_colors3.rename(columns={'Cluster Groups': 'Group'})

col_colors3 = col_colors3['Cluster Groups'].map(lut3)
print('COLORS 3')
print(col_colors3)



col_colors1 = species['Cluster Groups'].map(lut1)
print('COL COLORS')
print(col_colors1)
print(col_colors2_ANI)

print('ENDDDDDDDDDDDDDD')



#distance_matrix = pd.read_csv('matrix5.csv', index_col='Strains')

print('hey')
print(distance_matrix)
list_genomes = distance_matrix.index.tolist()
print('lol')
col_colors2_ANI = col_colors2_ANI.reindex(col_colors1.index, axis=1)
col_colors3_ANI = col_colors3_ANI.reindex(col_colors1.index, axis=1)

#COMBINING [col_colors1, col_colors2_ANI, col_colors3]
col_colors = pd.merge(col_colors1, col_colors2_ANI, left_index=True, right_index=True)
col_colors = pd.concat([col_colors1, col_colors2_ANI, col_colors3], axis=1)







# for ANI vmin=0.90, vmax=1.0,
g = sns.clustermap(distance_matrix,
                             cbar_pos=(.1, .2, .03, .4),
                             cmap="vlag",
                             figsize=(45, 45), vmin=4000, vmax=6000,
                             cbar_kws={"shrink": .1},
                             metric="correlation",
                             col_colors=[col_colors3, col_colors1, col_colors2_ANI, col_colors3_ANI], col_linkage=L, row_linkage=L
                             )

#annot=True

g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=8)
g.ax_row_dendrogram.set_visible(False)


print(g.ax_heatmap.get_xmajorticklabels())
plt.savefig( 'ANI_Try_v2/deneme_heatmap_V4_ani.pdf', dpi=360)


print(species)
species.to_csv('ANI_Try_v2/species_clusters_name_v4_ani')


input()
info = pd.read_csv('ANI_Try_v2/info_genomesv3', index_col='Strains')
df_outer = pd.merge(species, info, on='Strains', how='outer')
df_outer = df_outer.sort_values(by='Cluster Groups')
df_outer.to_csv('ANI_Try_v2/info_genomes_with_clusters_V3.csv')




