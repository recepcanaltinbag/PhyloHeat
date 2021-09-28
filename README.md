# PhyloHeat
Phylogenomic clustering and visualization for species determination. For now, PhyloHeat can be used to visualize the species and subgroup differentiation of bacterial genomes. There are two main clustering approach as comparing ANI values (with recommended thresholds) and creating a distance matrix based on shared protein families. The two approach is supporting each other most of the time. 

# Installation

Tested on Python 3.8

# The required packages 

- pandas-1.3.3
- seaborn-0.11.2
- scipy-1.7.1
- matplotlib-3.4.3
```
pip install pandas
pip install seaborn
```
(seaborn come with matplolib and scipy, so you do not need to install them again)


# Demo Run

You can try the initially prepared config file and demo folder to see what can be done with PhyloHeat. All you need is 
just to run PhyloHeat.py file.
```
python3 PhyloHeat.py
```
The outputs can be found in 'demo/output' folder. Example output:

![example_output](/example_outputs/screenshot.png)


# Example Run

All you need to do is changing the parameters in the 'config.conf' file.

## Input data formats

All inputs must be in the 'csv' format.

- genomes_info.csv :   UID (Assembly IDs) and Nice Name (You can use desired names) column names must be same. It just needed for final output because using just UIDs or any IDs can be confusing when looking and analyzing the heatmap.
- matrix.csv : A square matrix is needed. Genomes must be represented as UID (Assembly IDs) and their shared protein families must be in intersected cells.
- matrix_with_names.csv: Same matrix with previous file but UIDs are changed with Nice Names.
- Whole_ANI.csv: Whole genome ANI values in a square matrix form. UIDs must be used. 
- 16S_ANI.csv: 16S ANI values in a square matrix form. UIDs must be used. 

> You can use different IDs and names as long as they are compatible in all input files. 

## config.conf file format

### inputs
Inputs were explained previously, you must enter the path to config file for all inputs.

### settings
- ref_genome_name:  Reference genome name, the name must be same with the name in the genomes_info file.
- ref_genome_id: Reference genome id, the id must be same with the id in the other files. 
- user_min_ANI_Whole: Initially 0.85. It is for the defining the lowest ANI color.
- user_min_ANI_16S: Initially 0.95. It is for the defining the lowest ANI color. 
- cluster_strictness: Initially 0.2, can be modified to retrieve more stringent or relaxed clusters
- vmin: Minimum number of shared protein families (ınitially 4000) Just for coloring, you can try different values
- vmax: Maximum number of shared protein families (ınitially 6000) Just for coloring, you can try different values
- dpi: Output pdf dpi

### Outputs
- heatmat.pdf: output heatmap path
- species.csv: clustering results as csv file. Same numbers represents same classes. Different classes can be sign for different species.



