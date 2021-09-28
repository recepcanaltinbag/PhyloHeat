import pandas as pd
import seaborn as sns; sns.set(color_codes=True)

def take_ANI_column(share_matrix, ani_input, min_ANI, ref_genome_id, genome_info):
    print('ANI function')
    ANI_data = pd.read_csv(ani_input, index_col='Strains')
    strains_list = share_matrix['Strains'].tolist()

    ANI_out_list = []
    for strain in strains_list:
        try:
            ANI_out_list.append((strain, ANI_data.at[strain, ref_genome_id]))
        except:
            print('pass,', end='')
            ANI_out_list.append((strain, 0.95))

    df_ANI = pd.DataFrame(ANI_out_list, columns=['Strains', 'ANI'])
    print(df_ANI['ANI'].max(), df_ANI['ANI'].min())
    max_v = df_ANI['ANI'].max()
    max_v = 1.00
    min_v = min_ANI

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
    print('len', len(new_ANI_list), len(ANI_out_list))
    df_ANI = pd.DataFrame(new_ANI_list, columns=['Strains', 'ANI'])
    lut2 = dict(zip(df_ANI['ANI'].unique(), sns.color_palette("coolwarm", 6)))
    df_ANI = df_ANI.set_index('Strains')

    genome_id_name = pd.read_csv(genome_info)
    id_dict = dict(zip(genome_id_name['UID'], genome_id_name['Nice Name']))
    df_ANI = df_ANI.rename(index=id_dict)

    col_colors2_ANI = df_ANI['ANI'].map(lut2)
    print('END of the ANI function')
    return col_colors2_ANI
