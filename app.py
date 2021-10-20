import numpy as np
import pandas as pd

import scipy.stats as st
import seaborn as sns

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


# read in and process the dataframe
df = pd.read_csv("Halfmann_P1.csv")

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

# split dataframe into Luminex and functional data

df_luminex = pd.concat([df.loc[:, df.columns.str.contains("_")], df[["Treatment", "Challenge", "Sample"]]], axis=1)
df_func = df.loc[:, ~df.columns.str.contains("_")]
df_challenge = pd.concat([df.loc[:, df.columns.str.contains("pfu")], df[["Treatment", "Challenge", "Sample"]]], axis=1)

# make the data tidy for plotting

df_luminex_plot = df_luminex.melt(id_vars=["Treatment", "Challenge", "Sample"])
df_func_plot = df_func.melt(id_vars=["Treatment", "Challenge", "Sample"])
df_challenge_plot = df_func.melt(id_vars=["Treatment", "Challenge", "Sample"])

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


@pn.depends(ag_select.param.value,
            ig_fcr_choose.param.value)
def ag_strip_plot(antigen=df_luminex_plot.Ag.values[0], ig_or_fcr="Ig Titer"):
    
    if ig_or_fcr == "Ig Titer":
        df_small = df_luminex_plot.loc[(df_luminex_plot.Ag == antigen) & (df_luminex_plot.Ig_FcR.str.contains("Ig")), :]
        title = f"{antigen}-specific Antibodies"
    else:
        df_small = df_luminex_plot.loc[(df_luminex_plot.Ag == antigen) & (df_luminex_plot.Ig_FcR.str.contains("Fc")), :]
        title = f"{antigen}-specific Antibodies to FcRs"
        
    strip = hv.Scatter(
                data=df_small,
                kdims=['Ig_FcR'],
                vdims=["value", "Treatment", "Sample"],
            ).opts(
                color='Treatment',
                title=title,
                width=750,
                height=400,
                fontsize={'labels': 11, 'xticks': 10, 'yticks': 10}
            )
    
    p = hv.render(strip)
    p.add_layout(p.legend[0], 'right')
    p.toolbar_location = "above"
    
    return p


dash1 = pn.Row(pn.Column(pn.layout.VSpacer(), ag_select, pn.Spacer(height=50), ig_fcr_choose, pn.layout.VSpacer()), 
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
                vdims=["value", "Treatment", "Sample"],
            ).opts(
                color='Treatment',
                title=title,
                width=1100,
                height=400,
            )
    
    p = hv.render(strip)
    p.add_layout(p.legend[0], 'right')
    p.toolbar_location = "above"
    
    return p

dash2 = pn.Column(ab_fcr_select, pn.Spacer(height=30), ig_fcr_strip_plot)


def func_stripbox(df):
            
    assay_abbreviation_dict = {"ADCD": "Complement Deposition",
                               "mADCP": "J774A Monocyte Phagocytosis",
                               "hADCP": "THP-1 Monocyte Phagocytosis",
                               "hADNP": "Human Neutrophil Phagocytosis"}

    
    plots = []
    
    for func_assay in assay_abbreviation_dict.keys():

        df_small = df.loc[df.variable == func_assay]

        p = iqplot.stripbox(data=df_small, 
                            q="value", 
                            cats="Treatment", 
                            q_axis='y', 
                            jitter=True,
                            marker_kwargs={"size": 8, "line_color": "black"},
                            jitter_kwargs={'width': 0.25},
                            tooltips=[('Sample', '@Sample')],
                            toolbar_location="right",
                            height=400,
                            width=500,
                            y_axis_type="log",
                            x_axis_label="Group",
                            y_axis_label="MFI",
                            title=assay_abbreviation_dict[func_assay],
                           )

        plots.append(p)
        
    return bokeh.layouts.gridplot(plots, ncols=2, merge_tools=False)


abbreviation_dict = {"Lung (log10 pfu/g)": "Lung Viral Load",
                     "NT (log10 pfu/g)": "Nasal Turbinate Viral Load"}

def viral_load_stripbox(df):

    plots = []

    for location in abbreviation_dict.keys():

        df_small = df.loc[df.variable == location]
        
        # make box plot
        box = hv.BoxWhisker(
            data=df_small,
            kdims=['Treatment', 'Challenge'],
            vdims='value',
        ).opts(
            ylabel="Log10 (PFU/g)",
            whisker_color='gray',
            box_line_color="gray",
            box_fill_color="white",
            height=500,
            width=600,
            outlier_alpha=0,
            title = abbreviation_dict[location],
            fontsize={'labels': 11, 'xticks': 10, 'yticks': 10}
        )

        # extract bokeh object
        p = hv.render(box)
        p.toolbar_location = "above"
                
        plots.append(iqplot.strip(p=p,
                            data=df_small, 
                            q="value", 
                            cats=["Treatment", "Challenge"], 
                            q_axis='y', 
                            jitter=True,
                            marker_kwargs={"size": 8, "line_color": "black"},
                            jitter_kwargs={"width": 0.25},
                            toolbar_location="above",
                            height=450,
                            width=550,
                            y_axis_type="log",
                            x_axis_label="Challenge Variant",
                            y_axis_label="MFI",
                            title=abbreviation_dict[location],
                            tooltips=[('Sample', '@Sample')]
                           )
                    )
        
    return bokeh.layouts.gridplot(plots, ncols=2, merge_tools=False)


def zscore_heatmap(df):
    
    df_zscore = df.select_dtypes(include=np.number).apply(st.zscore)
    df_zscore["Treatment"] = df.Treatment
    df_zscore["Sample"] = df.Sample

    df_zscore.sort_values(by="Treatment", inplace=True)

    df_hm = df_zscore.melt(id_vars=["Sample", "Treatment"])
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
                              tools=[bokeh.models.HoverTool(tooltips=[('Treatment', '@Treatment'), ('Z-Score', '@value{0.000}'), ('', '@variable')])])

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


dash3 = zscore_heatmap(df)

tab1 = pn.Row(pn.layout.HSpacer(), pn.Column(dash1, pn.Spacer(height=20), dash2), pn.layout.HSpacer())

tab2 = pn.Row(pn.layout.HSpacer(), func_stripbox(df_func_plot), pn.layout.HSpacer())

tab3 = pn.Row(pn.layout.HSpacer(), viral_load_stripbox(df_challenge_plot), pn.layout.HSpacer())

tab4 = pn.Row(pn.layout.HSpacer(), dash3, pn.layout.HSpacer())

dashboard = pn.Tabs(('Luminex Plots', tab1), ("Functional Assays", tab2), ("Viral Loads", tab3), ("Z-Score Heatmap", tab4))

dashboard.servable()
