from configparser import ConfigParser
import pandas as pd
import scipy.cluster.hierarchy as sch
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt

from heatmap_functions import take_ANI_column

parser = ConfigParser()
parser.read("config.conf")

#INPUTS (files)
matrix_file = parser.get("inputs", "matrix")
matrix_with_names_file = parser.get("inputs", "matrix_with_names")
genome_info_file = parser.get("inputs", "genome_info")
file_16S_ANI = parser.get("inputs", "16S_ANI")
file_Whole_ANI = parser.get("inputs", "Whole_ANI")

#SETTINGS
ref_genome_name = parser.get("settings", "ref_genome_name")
ref_genome_id = parser.get("settings", "ref_genome_id")
user_min_ANI_Whole = parser.getfloat("settings", "user_min_ANI_Whole")
user_min_ANI_16S = parser.getfloat("settings", "user_min_ANI_16S")
cluster_strictness = parser.getfloat("settings", "cluster_strictness")
vmin = parser.getint("settings", "vmin")
vmax = parser.getint("settings", "vmax")
dpi = parser.getint("settings", "dpi")

#OUTPUTS
pdf_heatmap_out = parser.get("outputs", "pdf")
species_out = parser.get("outputs", "species")
clusters_out = parser.get("outputs", "clusters")


#READS
share_matrix = pd.read_csv(matrix_file)
col_colors_Whole_ANI = take_ANI_column(share_matrix, file_Whole_ANI, user_min_ANI_Whole, ref_genome_id, genome_info_file)
col_colors_16S_ANI = take_ANI_column(share_matrix, file_16S_ANI, user_min_ANI_16S, ref_genome_id, genome_info_file)

share_matrix = pd.read_csv(matrix_with_names_file)
share_matrix.set_index('Strains', inplace=True)


#DISTANCE, CLUSTERING
d = sch.distance.pdist(share_matrix)
L = sch.linkage(d, method='average')
# 0.2 can be modified to retrieve more stringent or relaxed clusters, it is good for ANI
clusters = sch.fcluster(L, cluster_strictness*d.max(), 'distance')


species_list = []
print(d,len(d))
print(0.2*d.max(), 'cluster treshold: ',d.max())

# clusters indicices correspond to incides of original df
for i,cluster in enumerate(clusters):
    species_list.append((share_matrix.index[i], cluster))

species = pd.DataFrame(species_list, columns=['Strains', 'Cluster Groups'])
lut1 = dict(zip(species['Cluster Groups'].unique(), sns.color_palette("Paired")))

species = species.set_index('Strains')

col_colors3 = species.copy()
col_colors3.loc[col_colors3['Cluster Groups'] > 1, "Cluster Groups"] = 2
lut3 = dict(zip(col_colors3['Cluster Groups'].unique(), sns.color_palette("Set2")))
col_colors3 = col_colors3['Cluster Groups'].map(lut3)
col_colors1 = species['Cluster Groups'].map(lut1)


list_genomes = share_matrix.index.tolist()
col_colors2_ANI = col_colors_16S_ANI.reindex(col_colors1.index, axis=1)
col_colors3_ANI = col_colors_Whole_ANI.reindex(col_colors1.index, axis=1)

col_colors = pd.concat([col_colors1, col_colors2_ANI, col_colors3], axis=1)

g = sns.clustermap(share_matrix,
                             cbar_pos=(.1, .2, .03, .4),
                             cmap="vlag",
                             figsize=(45, 45), vmin=vmin, vmax=vmax,
                             cbar_kws={"shrink": .1},
                             metric="correlation",
                             col_colors=[col_colors3, col_colors1, col_colors2_ANI, col_colors3_ANI], col_linkage=L, row_linkage=L
                             )


g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=8)
g.ax_row_dendrogram.set_visible(False)


plt.savefig(pdf_heatmap_out, dpi=dpi)

print(species)
species.to_csv(species_out)

print('END, please check output folder')



