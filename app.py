import numpy as np
import pandas as pd

import seaborn as sns
import scipy.stats as st

import bokeh.io
import bokeh.plotting
import bokeh.models

import holoviews as hv
hv.extension("bokeh")

import panel as pn
pn.extension()

from bokeh.models import (BasicTicker, ColorBar, ColumnDataSource,
                          LinearColorMapper, PrintfTickFormatter)
from bokeh.transform import transform


# read in and process the dataframe
df = pd.read_excel("Halfman_Systems SerologyP1_100121 PH.xlsx")

combined_names = []

for i in range(len(df)):
    
    if pd.isnull(df.Challenge.values[i]):
        combined_names.append(df.Treatment.values[i])
    else:
        combined_names.append(df.Treatment.values[i] + "_" + df.Challenge.values[i])
        
df["Treatment"] = combined_names
del df["Challenge"]



# make widgets

ag_lst = [name.split("_")[0] for name in df.columns[3:]]
ab_fcr_lst = [name.split("_")[1] for name in df.columns[3:]]

# make selection widgets
ag_select = pn.widgets.Select(name='Select Antigen', options=list(np.unique(ag_lst)))
ab_fcr_select = pn.widgets.Select(name='Select Ig or FcR', options=list(np.unique(ab_fcr_lst)))
ig_fcr_choose = pn.widgets.RadioButtonGroup(options=['Ig Titer', 'FcR Binding'], button_type='primary')

df_plot = df.melt(id_vars=["Treatment", "Sample"])

antigens, igs_fcrs = list(zip(*df_plot.variable.str.split("_")))
df_plot["Ag"] = antigens
df_plot["Ig_FcR"] = igs_fcrs


@pn.depends(ag_select.param.value,
            ig_fcr_choose.param.value)
def ag_strip_plot(antigen=df_plot.Ag.values[0], ig_or_fcr="Ig Titer"):
    
    if ig_or_fcr == "Ig Titer":
        df_small = df_plot.loc[(df_plot.Ag == antigen) & (df_plot.Ig_FcR.str.contains("Ig")), :]
        title = f"Titers of {antigen}-specific Antibodies"
    else:
        df_small = df_plot.loc[(df_plot.Ag == antigen) & (df_plot.Ig_FcR.str.contains("Fc")), :]
        title = f"Binding of {antigen}-specific Antibodies to FcRs"
        
    strip = hv.Scatter(
                data=df_small,
                kdims=['Ig_FcR'],
                vdims=["value", "Treatment", "Sample"],
            ).opts(
                color='Treatment',
                line_color="black",
                alpha=0.8,
                jitter=0.3,
                size=9,
                xlabel="",
                ylabel='Median Fluorescence',
                title=title,
                logy=True,
                width=700,
                height=450,
                legend_position="bottom",
                tools=[bokeh.models.HoverTool(tooltips=[('Sample', '@Sample'), ('Fluor', '@value{int}')])],
                cmap=['#1f77b4', 'darkorange', 'green'],
                fontsize={'labels': 11, 'xticks': 10, 'yticks': 10}
            )
    
    return strip


dash1 = pn.Row(pn.Column(pn.layout.VSpacer(), ag_select, pn.Spacer(height=50), ig_fcr_choose, pn.layout.VSpacer()), 
               pn.Spacer(width=40),
               ag_strip_plot)

@pn.depends(ab_fcr_select.param.value)
def ig_fcr_strip_plot(ig_or_fcr=df_plot.Ig_FcR.values[0]):
    
    df_small = df_plot.loc[(df_plot.Ig_FcR == ig_or_fcr), :]
    
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
                line_color="black",
                alpha=0.8,
                jitter=0.3,
                size=9,
                xlabel="",
                ylabel='Median Fluorescence',
                title=title,
                logy=True,
                width=900,
                height=450,
                legend_position="bottom",
                tools=[bokeh.models.HoverTool(tooltips=[('Sample', '@Sample'), ('Fluor', '@value{int}')])],
                cmap=['#1f77b4', 'darkorange', 'green'],
                fontsize={'labels': 11, 'xticks': 10, 'yticks': 10}
            )
    
    return strip

dash2 = pn.Column(ab_fcr_select, pn.Spacer(height=30), ig_fcr_strip_plot)


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

    p = bokeh.plotting.figure(width=1350, 
                              height=600, 
                              # title="Hamster S2 Z-Score Data",
                              x_range=list(df_hm.variable.unique()), 
                              y_range=list(df_hm.Sample.unique())[::-1],
                              toolbar_location=None, x_axis_location="above",
                              tools=[bokeh.models.HoverTool(tooltips=[('Group', '@Treatment'), ('Z-Score', '@value{0.000}'), ('', '@variable')])])

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
tab2 = pn.Row(pn.layout.HSpacer(), dash3, pn.layout.HSpacer())

dashboard = pn.Tabs(('Strip Plots', tab1), ("Z-Score Heatmap", tab2))

dashboard.servable()
