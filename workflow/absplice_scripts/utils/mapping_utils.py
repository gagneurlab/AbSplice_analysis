import pandas as pd

def subset_tissues(df, tissue_map, chosen_tissue):
    df['tissue_main'] = df['tissue'].map(tissue_map)
    if chosen_tissue in tissue_map.values():
        tissue_column = 'tissue_main'
    elif chosen_tissue in tissue_map.keys():
        tissue_column = 'tissue'
    else:
        raise KeyError('%s is not in provided tissues', chosen_tissue) 
    return df[df[tissue_column] == chosen_tissue]
