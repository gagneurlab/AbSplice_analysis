from tqdm import tqdm
import re

missing = list()
df = list()
for path in tqdm(snakemake.input['performance']):
    try:
        _df = pd.read_csv(path).assign(
            var_category = re.search('var_category=(.*)\/outlier_category', path).group(1),
            outlier_category = re.search('outlier_category=(.*)\/performance', path).group(1),
        )

        df.append(_df)
    except:
        missing.append(path)
df = pd.concat(df)

for model in tqdm(model_order):
    # set start precision to 0
    df.loc[
        (df['model'] == model)
        & (df['precision'] == 1),
        'precision'
    ] = 0
    # set end precision to 0
    df.loc[
        (df['model'] == model)
        & (df['recall'] == 1),
        'precision'
    ] = 0

df = df.replace(snakemake.params['replace_dict'])

df.to_csv(snakemake.output['performance_plot'], index=False)