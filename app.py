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

from utils import *


# read in and process the dataframe
#df = pd.read_csv("Halfmann_P1.csv")
df = pd.read_csv("Halfmann_P1_bgSubtracted.csv")

# remove empty columns
df = df[df.columns[~df.columns.str.contains("Unnamed")]]

# make widgets

luminex_cols = df.columns[df.columns.str.contains("_")]

ag_lst = [name.split("_")[0] for name in luminex_cols]
ab_fcr_lst = [name.split("_")[1] for name in luminex_cols]

# make selection widgets
ag_select = pn.widgets.Select(name='Select Antigen', options=list(np.unique(ag_lst)))
ab_fcr_select = pn.widgets.Select(name='Select Ig or FcR', options=list(np.unique(ab_fcr_lst)))
ig_fcr_choose = pn.widgets.RadioButtonGroup(options=['Ig Titer', 'FcR Binding'], button_type='primary')
box_plot_toggle = pn.widgets.Toggle(name="Show Box Plot", button_type="success")

# split dataframe into Luminex and functional data

df_luminex = pd.concat([df.loc[:, df.columns.str.contains("_")], df[["Immunization", "Challenge", "Sample"]]], axis=1)
df_func = pd.concat([df.loc[:, df.columns.str.contains("AD")], df[["Immunization", "Challenge", "Sample"]]], axis=1)
df_challenge = pd.concat([df.loc[:, df.columns.str.contains("pfu")], df[["Immunization", "Challenge", "Sample"]]], axis=1)

# make the data tidy for plotting

df_luminex_plot = df_luminex.melt(id_vars=["Immunization", "Challenge", "Sample"])
df_func_plot = df_func.melt(id_vars=["Immunization", "Challenge", "Sample"])
df_challenge_plot = df_challenge.melt(id_vars=["Immunization", "Challenge", "Sample"])

antigens, igs_fcrs = list(zip(*df_luminex_plot.variable.str.split("_")))
df_luminex_plot["Ag"] = antigens
df_luminex_plot["Ig_FcR"] = igs_fcrs

opts.defaults(
    opts.Scatter(
        size=9,
        line_color="black",
        alpha=0.8,
        jitter=0.3,
        logy=True,
        xlabel="",
        ylabel='MFI',
        cmap=bokeh.palettes.Category10[10],
        tools=[bokeh.models.HoverTool(tooltips=[('Sample', '@Sample'), ('Fluor', '@value{int}')])],
        fontsize={'labels': 11, 'xticks': 10, 'yticks': 10}
    ),
)

assay_abbreviation_dict = {"ADCD": "Complement Deposition",
                           "mADCP": "J774A Monocyte Phagocytosis",
                           "hADCP": "THP-1 Monocyte Phagocytosis",
                           "hADNP": "Human Neutrophil Phagocytosis"}

compartment_abbreviation_dict = {"Lung (log10 pfu/g)": "Lung Viral Load",
                                 "NT (log10 pfu/g)": "Nasal Turbinate Viral Load"}


@pn.depends(ag_select.param.value,
            ig_fcr_choose.param.value,
            box_plot_toggle.param.value)
def ag_strip_plot(antigen=df_luminex_plot.Ag.values[0], 
                  ig_or_fcr="Ig Titer",
                  show_box_plot=False):
    
    if ig_or_fcr == "Ig Titer":
        df_small = df_luminex_plot.loc[(df_luminex_plot.Ag == antigen) & (df_luminex_plot.Ig_FcR.str.contains("Ig")), :]
        title = f"{antigen}-specific Antibodies"
    else:
        df_small = df_luminex_plot.loc[(df_luminex_plot.Ag == antigen) & (df_luminex_plot.Ig_FcR.str.contains("Fc")), :]
        title = f"{antigen}-specific Antibodies to FcRs"

    if show_box_plot:
        
        # make box plot
        box = hv.BoxWhisker(
            data=df_small,
            kdims=["Ig_FcR", "Immunization"],
            vdims='value',
        ).opts(
            ylabel="MFI",
            whisker_color='gray',
            box_line_color="gray",
            box_fill_color="white",
            height=500,
            width=800,
            outlier_alpha=0,
            logy=True,
            title = title,
            fontsize={'labels': 11, 'xticks': 10, 'yticks': 10}
        )

        # extract bokeh object
        p = hv.render(box)
    
        iqplot.strip(p=p,
                    data=df_small, 
                    q="value", 
                    cats=["Ig_FcR", "Immunization"], 
                    q_axis='y', 
                    color_column="Immunization",
                    jitter=True,
                    marker_kwargs={"size": 8, "line_color": "black"},
                    jitter_kwargs={"width": 0.25},
                    tooltips=[('Sample', '@Sample'), ("Fluor", '@value')]
                    )
        
    else:
        
        p = iqplot.strip(data=df_small, 
                    q="value", 
                    cats=["Ig_FcR", "Immunization"], 
                    q_axis='y', 
                    color_column="Immunization",
                    jitter=True,
                    marker_kwargs={"size": 8, "line_color": "black"},
                    jitter_kwargs={"width": 0.25},
                    title=title,
                    height=500,
                    width=800,
                    y_axis_type="log",
                    y_axis_label="MFI",
                    tooltips=[('Sample', '@Sample'), ("Fluor", '@value')]
                )

    p.toolbar_location = "above"
    
    return p


dash1 = pn.Row(pn.Column(pn.layout.VSpacer(), ag_select, pn.Spacer(height=30), ig_fcr_choose, pn.Spacer(height=30), box_plot_toggle, pn.layout.VSpacer()), 
               pn.Spacer(width=40),
               ag_strip_plot)

@pn.depends(ab_fcr_select.param.value)
def ig_fcr_strip_plot(ig_or_fcr=df_luminex_plot.Ig_FcR.values[0]):
    
    df_small = df_luminex_plot.loc[(df_luminex_plot.Ig_FcR == ig_or_fcr), :]
    
    if "Ig" in ig_or_fcr:
        title = f"Antigen-Specific {ig_or_fcr} Titers"
    else:
        title = f"{ig_or_fcr} Binding Across Antigens"
    
    strip = hv.Scatter(
                data=df_small,
                kdims=['Ag'],
                vdims=["value", "Immunization", "Sample"],
            ).opts(
                color='Immunization',
                title=title,
                width=1100,
                height=400,
            )
    
    p = hv.render(strip)
    p.add_layout(p.legend[0], 'right')
    p.toolbar_location = "above"
    p.legend.title = "Immunization"
    
    return p


# Interactive plots for data exploration
dash2 = pn.Column(ab_fcr_select, pn.Spacer(height=30), ig_fcr_strip_plot)
tab1 = pn.Row(pn.layout.HSpacer(), pn.Column(dash1, pn.Spacer(height=20), dash2), pn.layout.HSpacer())

# Functional assays
tab2 = pn.Row(pn.layout.HSpacer(), stripbox(df_func_plot, ["Immunization"], assay_abbreviation_dict), pn.layout.HSpacer())

# Viral Loads
tab3 = pn.Row(pn.layout.HSpacer(), stripbox(df_challenge_plot, ["Immunization", "Challenge"], compartment_abbreviation_dict), pn.layout.HSpacer())

# Z Scores of all the data
df_log = log_transform_data(df)
tab4 = pn.Row(pn.layout.HSpacer(), zscore_heatmap(df_log), pn.layout.HSpacer())

# Correlations to identify features that might be significant predictors of viral load
df_corr_drop = df_log[df_log.columns[~df_log.columns.str.contains("|".join(["pfu", "Ig"]))]]

lung_correlation = df_corr_drop.corrwith(df_log["Lung (log10 pfu/g)"])
NT_correlation = df_corr_drop.corrwith(df_log["NT (log10 pfu/g)"])

lung_corr_df = pd.DataFrame({"Feature": lung_correlation.index, "Correlation": lung_correlation.values}).sort_values(by="Correlation")
NT_corr_df = pd.DataFrame({"Feature": NT_correlation.index, "Correlation": NT_correlation.values}).sort_values(by="Correlation")

# based on LOO cross-validation
num_components = 3

# Linear model of the features computed above
tab5 = pn.Row(pn.layout.HSpacer(),
              pn.Column(pn.Row(
                                pn.Column(pn.Spacer(height=30), correlation_plot(lung_corr_df, "Lung Viral Load Correlates")), 
                                pls_regression(df_log, lung_corr_df, -2, num_components, "Lung Viral Load Model")),
                        pn.Row(
                                pn.Column(pn.Spacer(height=30), correlation_plot(NT_corr_df, "Nasal Turbinate Viral Load Correlates")),
                                pls_regression(df_log, NT_corr_df, -1, num_components, "Nasal Turbinate Viral Load Model"))
                       ),
              pn.layout.HSpacer())

# Flower plots
col_names = ["Immunization", "Sample"] + list(df.columns[df.columns.str.contains("|".join(["S2", "AD"]))])
flower_plots = polar_area_plot(df_log, col_names)
tab6 = pn.Row(pn.layout.HSpacer(), flower_plots, pn.layout.HSpacer())

dashboard = pn.Tabs(('Luminex Plots', tab1), 
                    ("Functional Assays", tab2), 
                    ("Viral Loads", tab3), 
                    ("Z-Score Heatmap", tab4),
                    ("Significant Features", tab5),
                    ("Flower Plots", tab6))

dashboard.servable()