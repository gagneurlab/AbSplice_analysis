import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import average_precision_score, precision_recall_curve
from absplice_scripts.visualization.utils import plot_x_eq_y
from tqdm import tqdm


def plot_roc_curve(recall, fpr):
    plot_curve(fpr, recall)
    plt.xlabel('FPR')
    plt.ylabel('Recall')
    plot_x_eq_y()


def plot_curve(recall, other_score, xlim=(0, 1.01), ylim=(0, 1.01), color=None):
    if color == None:
        plt.step(recall, other_score, where='post')
    else:
        plt.step(recall, other_score, where='post', color=color)
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)


def plot_precision_recall_curve(recall, precision, color=None):
    if color == None:
        plot_curve(recall, precision)
    else:
        plot_curve(recall, precision, color=color)
    plt.xlabel('Recall')
    plt.ylabel('Precision')


def plot_multi_precision_recall_curve(benchmark, ndigits=4, colors=None):
    performance = dict()
    legends = list()

    count = 0
    for k, (target, prediction) in benchmark.items():
        precision, recall, threshold = precision_recall_curve(
            target, prediction)
        if colors == None:
            plot_precision_recall_curve(recall, precision)
        else:
            plot_precision_recall_curve(recall, precision, color=colors[count])
            count += 1
        aps = average_precision_score(target, prediction)
        aps = round(aps, ndigits)
        legends.append(f'{k}, AP={aps}')
        performance[k] = {
            'precision': precision,
            'recall': recall,
            'threshold': threshold,
            'aps': aps
        }
    plt.legend(legends)

    return performance


def get_dataframe_for_ggplot(performance):
    df_performance = pd.DataFrame(performance).transpose().reset_index().rename(columns={'index': 'model'})[['model', 'precision', 'recall']]
    df_performance = df_performance.set_index('model').apply(pd.Series.explode).reset_index()
    df_performance['precision'] = df_performance['precision'].astype(float)
    df_performance['recall'] = df_performance['recall'].astype(float)
    
    df_aps = pd.DataFrame(performance).transpose().reset_index().rename(columns={'index': 'model'})[['model', 'aps']]
    df_performance = df_performance.set_index('model').join(df_aps.set_index('model')).reset_index()
    return df_performance


def get_performance(df, outlier_column, model_dict, colors=None):

    try:
        df[outlier_column] = df[outlier_column].fillna(0).astype(int)
    except:
        if 'True' in df[outlier_column].values:
            df[outlier_column] = df[outlier_column] == 'True'
        df[outlier_column] = df[outlier_column].astype(int)
        
    plot_dict = dict()
    for model in tqdm(model_dict.keys()):
        plot_dict[model] = (
            df[outlier_column], 
            np.abs(df[model_dict[model]]).fillna(0)
        )

    _performance = plot_multi_precision_recall_curve(
        plot_dict,
        ndigits=3,
        colors=colors
    )
    
    return get_dataframe_for_ggplot(_performance), _performance


def get_performance_multi(df, model_dict, outlier_splice='outlier', outlier_gene_expr='outlier_gene_expr', target='splice'):

    plot_dict = dict()
    for model in model_dict.keys():
        
        if target == 'splice':
            plot_dict[model] = (
                df[outlier_splice].fillna(0).astype(int), 
                np.abs(df[model_dict[model]]).fillna(0)
            )
        elif target == 'gene_expr':
            plot_dict[model] = (
                df[outlier_gene_expr].fillna(0).astype(int), 
                np.abs(df[model_dict[model]]).fillna(0)
            )
        elif target == 'both':
            plot_dict[model] = (
                (df[outlier_splice].fillna(0).astype(int)) | (df[outlier_gene_expr].fillna(0).astype(int)), 
                np.abs(df[model_dict[model]]).fillna(0)
            )
        else:
            raise NotImplementedError()

    _performance = plot_multi_precision_recall_curve(
        plot_dict,
        ndigits=3
    )
    
    return get_dataframe_for_ggplot(_performance), _performance


def jackknife_performance(df_benchmark, model_dict, outlier_column='outlier'):
    individuals = list(set(df_benchmark.index.get_level_values('sample')))

    model_jack = dict()
    for model, model_name in model_dict.items():
        print(model, model_name)
        aps = list()
        for indiv in tqdm(individuals):
            df, performance = get_performance(
                df = df_benchmark[df_benchmark.index.get_level_values('sample') != indiv], 
                 outlier_column = outlier_column,
                 model_dict = dict({model: model_name})
            )
            assert len(list(set(df['aps']))) == 1
            aps.append(list(set(df['aps']))[0])
        model_jack[model] = aps
    df_jack = pd.DataFrame(model_jack)
    
    df_stats = df_jack\
        .melt(var_name='model', value_name='auPRC')\
        .groupby('model')['auPRC']\
        .agg(['count', 'mean', 'std', 'sem'])\
        .reset_index()
    
    return df_jack, df_stats


def plot_rank_curve(target=None, score=None,
                    ranks=None, target_ranks=None, positive_label=1):
    if target is not None and score is not None:
        ranks = [
            i
            for i, (s, t) in enumerate(sorted(zip(score, target), reverse=True))
            if t == positive_label
        ]
        target_ranks = list(range(len(ranks)))
    elif ranks is None and target_ranks is None:
        raise ValueError('Either target and score or rank,'
                         'target_rank argument need to defined.')
    ax = plot_curve(ranks, target_ranks, xlim=None, ylim=None)
    plt.xlabel('prediction rank')
    plt.ylabel('target rank')
    return ax, ranks, target_ranks



def get_score_df(results, wildcards):
    df_benchmark_gene = pd.concat([
        pd.DataFrame(results, columns=['model', 'Average Precision Score']),
        pd.DataFrame(wildcards, columns=['tissue', 'CAT tissue'])
    ], axis=1)
    
    df_benchmark_gene = df_benchmark_gene[
        ~( 
            (df_benchmark_gene['model'].str.contains('CAT'))
            & (df_benchmark_gene['tissue'] == df_benchmark_gene['CAT tissue'])
        )
    ]
    
    cat_mapper = {
        'Cells_EBV_transformed_lymphocytes': 'Lymphocytes',
        'Cells_Transformed_fibroblasts': 'Fibroblasts',
        'Whole_Blood': 'Whole Blood',
    }

    model_mapper = dict()
    for model in list(set(df_benchmark_gene['model'])):
        model_new = model
        for cat_old, cat_new in cat_mapper.items():
            model_new = model_new.replace(cat_old, cat_new)
        model_mapper[model] = model_new

    df_benchmark_gene['model'] = df_benchmark_gene['model'].apply(lambda x: model_mapper[x])
    
    df_benchmark_gene['CAT tissue'] = df_benchmark_gene['CAT tissue'].replace({
        'Cells_Transformed_fibroblasts': 'Fibroblasts',
        'Cells_EBV_transformed_lymphocytes': 'Lymphocytes',
        'Whole_Blood': 'Whole Blood'
    })
    
    scores_mean = pd.DataFrame(df_benchmark_gene.groupby('model').mean())['Average Precision Score'].to_dict()
    scores_std = pd.DataFrame(df_benchmark_gene.groupby('model').std())['Average Precision Score'].to_dict()
    
    df_benchmark_gene['mean'] = df_benchmark_gene['model'].map(lambda x: scores_mean[x])
    df_benchmark_gene['std'] = df_benchmark_gene['model'].map(lambda x: scores_std[x])
    
    return df_benchmark_gene