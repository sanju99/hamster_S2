import numpy as np
import pandas as pd

import scipy.stats as st
import seaborn as sns
import matplotlib.pyplot as plt

import bokeh.io
import bokeh.plotting
import bokeh.models

import holoviews as hv
hv.extension("bokeh")
from holoviews import opts

import iqplot

import panel as pn
pn.extension()

from bokeh.models import BasicTicker, ColorBar, LinearColorMapper, PrintfTickFormatter
from bokeh.transform import transform
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import LeaveOneOut, cross_val_score


def log_transform_data(df):
    '''
    This function log-transforms the non-log-transformed data in the dataframe
    '''
    
    return pd.concat([df[["Challenge", "Immunization", "Sample"]], 
                      np.log10(df[df.columns[df.columns.str.contains("|".join(["_", "AD"]))]]), 
                      df[df.columns[df.columns.str.contains("pfu")]]], axis=1)



def stripbox(df_plot, kdims_lst, abbrev_dict):
    
    plots = []
    
    for key in abbrev_dict.keys():

        df_small = df_plot.loc[df_plot.variable == key]

       # make box plot
        box = hv.BoxWhisker(
            data=df_small,
            kdims=kdims_lst,
            vdims='value',
        ).opts(
            ylabel="MFI",
            whisker_color='gray',
            box_line_color="gray",
            box_fill_color="white",
            height=500,
            width=600,
            outlier_alpha=0,
            title = abbrev_dict[key],
            fontsize={'labels': 11, 'xticks': 10, 'yticks': 10}
        )

        # extract bokeh object
        p = hv.render(box)
        p.toolbar_location = "above"
                
        plots.append(iqplot.strip(p=p,
                            data=df_small, 
                            q="value", 
                            cats=kdims_lst, 
                            q_axis='y', 
                            jitter=True,
                            marker_kwargs={"size": 8, "line_color": "black"},
                            jitter_kwargs={"width": 0.25},
                            height=450,
                            width=550,
                            y_axis_type="log",
                            x_axis_label="Challenge Variant",
                            y_axis_label="MFI",
                            title=abbrev_dict[key],
                            tooltips=[('Sample', '@Sample')],
                           )
                    )
        
    return bokeh.layouts.gridplot(plots, ncols=2, merge_tools=False)


def polar_area_plot(df_log, col_names):
    '''
    Make a polar area (flower) plot for each hamster
    '''

    df_split = df_log[col_names].sort_values(by="Immunization")
    new_ordering = ["Immunization", "Sample", "S2_IgG1", "S2_IgG2", "S2_IgA", "S2_IgM", 'S2_FcγRIIb', 'S2_FcγRIII', 'S2_FcγRIV', 'ADCD', 'mADCP', 'hADCP', 'hADNP']
    df_split = df_split[new_ordering]
    
    titles = df_split.Sample.values
    del df_split["Immunization"]
    del df_split["Sample"]
    
    variables = [df_split.columns[i].split("_")[-1] if "_" in df_split.columns[i] else df_split.columns[i] for i in range(len(df_split.columns))]
    
    # Compute pie slices
    N = len(df_split.columns)
    
    theta = np.linspace(0, 2 * np.pi, N, endpoint=False)
    width = 2 * np.pi / N

    colors = ["purple" if "AD" in variables[k] else "orange" for k in range(len(df_split.columns))]

    sns.set_style("dark")

    fig, ax = plt.subplots(4, 4, figsize=(18, 18), subplot_kw={'projection': 'polar'}, constrained_layout=True)
    ax = ax.flatten()
    
    for i in range(len(df_split)):

        radii = df_split.iloc[i, :].values
        
        ax[i].set_title(titles[i], pad=50, fontsize=16)
        ax[i].tick_params(axis='x', pad=10)

        ax[i].bar(theta, 
               radii, 
               width=width,
               bottom=0, 
               color=colors,
               edgecolor="black")

        ax[i].set_yticklabels([])
        ax[i].set_ylim([0, int(np.ceil(np.max(np.max(df_split))))])

        ax[i].set_xticks(theta)
        ax[i].set_xticklabels(variables)        

        ax[i].tick_params(axis='both', which='major', labelsize=12)

    plt.close()
    return pn.pane.Matplotlib(fig)


def zscore_heatmap(df):
    
    df_zscore = df.select_dtypes(include=np.number).apply(st.zscore)
    
    # neat ordering of the x axis
    heatmap_order = list(df_zscore.columns[df_zscore.columns.str.contains("IgG1")]) + \
                    list(df_zscore.columns[df_zscore.columns.str.contains("IgG2")]) + \
                    list(df_zscore.columns[df_zscore.columns.str.contains("IgA")]) + \
                    list(df_zscore.columns[df_zscore.columns.str.contains("IgM")]) + \
                    list(df_zscore.columns[df_zscore.columns.str.contains("FcγRIIb")]) + \
                    list(df_zscore.columns[df_zscore.columns.str.contains("FcγRIII")]) + \
                    list(df_zscore.columns[df_zscore.columns.str.contains("FcγRIV")]) + \
                    list(df_zscore.columns[df_zscore.columns.str.contains("pfu")])
    
    df_zscore = pd.concat([df[["Challenge", "Immunization", "Sample"]], df_zscore[heatmap_order], ], axis=1)    
    df_zscore.sort_values(by=["Immunization", "Challenge"], inplace=True)
    
    df_hm = df_zscore.melt(id_vars=["Sample", "Immunization", "Challenge"])
    source = bokeh.models.ColumnDataSource(df_hm)

    color_num_max = np.max([abs(df_hm.value.min()), abs(df_hm.value.max())])

    colors = list(sns.diverging_palette(220, 20).as_hex())
    mapper = LinearColorMapper(palette=colors, 
                               low=-color_num_max,
                               high=color_num_max)

    p = bokeh.plotting.figure(width=1400, 
                              height=600, 
                              # title="Hamster S2 Z-Score Data",
                              x_range=list(df_hm.variable.unique()), 
                              y_range=list(df_hm.Sample.unique())[::-1],
                              toolbar_location=None, x_axis_location="above",
                              tools=[bokeh.models.HoverTool(tooltips=[('Immunization', '@Immunization'), 
                                                                      ('Z-Score', '@value{0.000}'), 
                                                                      ('Feature', '@variable'),
                                                                      ("Challenge", "@Challenge")])])

    p.rect(x="variable", y="Sample", width=1, height=1, source=source,
           line_color="white", 
           fill_color=transform('value', mapper),
          )

    color_bar = ColorBar(color_mapper=mapper,
                         ticker=BasicTicker(desired_num_ticks=len(colors)),
                        )

    p.add_layout(color_bar, 'right')

    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None

    p.xaxis.major_label_text_font_size = "10px"
    p.yaxis.major_label_text_font_size = "12px"

    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = 1.0
    
    return p


def correlation_plot(corr_df, title=""):
    
    corr_plot = hv.Bars(corr_df, 
            kdims=['Feature'], 
            vdims=["Correlation"]
           ).opts(height=450,
                  width=1000,
                  ylim=(0, np.min(corr_df.Correlation)-0.05),
                  xrotation=90,
                  fill_color="#1f77b4",
                  ylabel="Spearman Correlation",
                  title=title,
                  tools=[bokeh.models.HoverTool(tooltips=[('', '@Correlation{0.0000}'), ('', '@Feature')])],
        )

    p = hv.render(corr_plot)
    p.toolbar_location = "above"
    p.title.text_font_size = "16px"

    return p


def pls_regression(df_log, viral_load_corr, index, num_components=3, title=""):
    
    df_model = df_log[df_log.columns[df_log.columns.str.contains("|".join(["_", "AD"]))]]
    
    X = df_model[viral_load_corr.Feature].values
    y = df_log.iloc[:, index].values

    pls2 = PLSRegression(n_components=num_components)
    pls2.fit(X, y)
    Y_pred = pls2.predict(X)
    
    pearson = st.pearsonr(Y_pred.flatten(), df_log.iloc[:, index].values)[0]
    
    df_plot = pd.DataFrame({"pred": Y_pred.flatten(), "actual": df_log.iloc[:, index].values})
        
    fig, ax = plt.subplots(figsize=(8, 7))
    
    sns.set_style("white")
    sns.regplot(x="pred", y="actual", data=df_plot, ax=ax, ci=95)
    ax.set_xlabel("Predicted", fontsize=12)
    ax.set_ylabel("Actual", fontsize=12)
    ax.set_title(f"{title} \n Components = {num_components}, Pearson Correlation = {round(pearson, 4)}", fontsize=16)
    
    plt.close()
    
    return pn.pane.Matplotlib(fig)


def loo_cv(df_log, lung_corr_df, NT_corr_df, full_features=False):
    
    df_model = df_log[df_log.columns[df_log.columns.str.contains("|".join(["_", "AD"]))]]

    if full_features:
        X_lung = df_model.values
        X_NT = df_model.values
    else:
        X_lung = df_model[lung_corr_df.Feature].values
        X_NT = df_model[NT_corr_df.Feature].values
    
    # maximum features for testing here is N / 2 = 8. Realistically, the maximum number of predictors should be more like 4
    max_features = int(len(df_model)/2)

    loo_lung = LeaveOneOut()
    loo_NT = LeaveOneOut()
    
    loo_lung.get_n_splits(X_lung)
    loo_NT.get_n_splits(X_NT)

    scores_lst_lung = []
    scores_lst_NT = []

    y_lung = df_log.iloc[:, -2].values
    y_NT = df_log.iloc[:, -1].values
    
    # With N samples and uncorrelated features, a reasonable number of predictors is N-1
    # With N samples and uncorrelated features, a reasonable number of predictors is sqrt(N)    
    for i in range(1, max_features+1):

        pls2 = PLSRegression(n_components=i)

        scores_lung = cross_val_score(pls2, X_lung, y_lung, scoring='neg_mean_absolute_error',
                             cv=loo_lung, n_jobs=-1)
        scores_NT = cross_val_score(pls2, X_NT, y_NT, scoring='neg_mean_absolute_error',
                             cv=loo_NT, n_jobs=-1)

        # Mean of mean absolute error -- compute mean absolute error, then take mean of the absolute value
        scores_lst_lung.append(np.mean(abs(scores_lung)))
        scores_lst_NT.append(np.mean(abs(scores_NT)))
        
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(np.arange(1, max_features+1), scores_lst_lung, label="Lung")
    ax.plot(np.arange(1, max_features+1), scores_lst_NT, label="NT")
    
    ax.set_title("Mean Absolute Error in Leave-One-Out Cross Validation")
    ax.set_xlabel("Number of Predictors")
    ax.set_ylabel("Mean MAE")

    plt.legend()
    
    plt.close()
    return pn.pane.Matplotlib(fig)