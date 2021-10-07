import numpy as np
import pandas as pd

import seaborn as sns
import scipy.stats as st

import bokeh.io
import bokeh.plotting
import bokeh.models

import holoviews as hv
hv.extension("bokeh")
from holoviews import opts

import panel as pn
pn.extension()

from bokeh.models import BasicTicker, ColorBar, LinearColorMapper, PrintfTickFormatter
from bokeh.transform import transform


# read in and process the dataframe
df = pd.read_csv("Halfmann_P1.csv")

combined_names = []

for i in range(len(df)):
    
    if pd.isnull(df.Challenge.values[i]):
        combined_names.append(df.Treatment.values[i])
    else:
        combined_names.append(df.Treatment.values[i] + "_" + df.Challenge.values[i])
        
df["Treatment"] = combined_names
del df["Challenge"]



# make widgets

luminex_cols = df.columns[df.columns.str.contains("_")]

ag_lst = [name.split("_")[0] for name in luminex_cols]
ab_fcr_lst = [name.split("_")[1] for name in luminex_cols]

# make selection widgets
ag_select = pn.widgets.Select(name='Select Antigen', options=list(np.unique(ag_lst)))
ab_fcr_select = pn.widgets.Select(name='Select Ig or FcR', options=list(np.unique(ab_fcr_lst)))
ig_fcr_choose = pn.widgets.RadioButtonGroup(options=['Ig Titer', 'FcR Binding'], button_type='primary')

# split dataframe into Luminex and functional data

df_luminex = pd.concat([df.loc[:, df.columns.str.contains("_")], df[["Treatment", "Sample"]]], axis=1)
df_func = df.loc[:, ~df.columns.str.contains("_")]

# make the data tidy for plotting

df_luminex_plot = df_luminex.melt(id_vars=["Treatment", "Sample"])
df_func_plot = df_func.melt(id_vars=["Treatment", "Sample"])

antigens, igs_fcrs = list(zip(*df_luminex_plot.variable.str.split("_")))
df_luminex_plot["Ag"] = antigens
df_luminex_plot["Ig_FcR"] = igs_fcrs

# widget for the functional assay visualization
assay_select = pn.widgets.Select(name="Select Functional Assay", options=list(df_func_plot.variable.unique()))
groups_toggle = pn.widgets.Toggle(name='Combine Challenged Groups', button_type='success')

treatment_split2 = ["Control" if 'Control' in df_func_plot.Treatment.values[i] else 'S2 Immunized' for i in range(len(df_func_plot))]
df_func_plot["Group"] = treatment_split2


opts.defaults(
    opts.Scatter(
        size=9,
        line_color="black",
        alpha=0.8,
        jitter=0.3,
        logy=True,
        xlabel="",
        ylabel='MFI',
        cmap=['#1f77b4', 'darkorange', 'green'],
        tools=[bokeh.models.HoverTool(tooltips=[('Sample', '@Sample'), ('Fluor', '@value{int}')])],
        legend_position="bottom"
    ),
)


@pn.depends(ag_select.param.value,
            ig_fcr_choose.param.value)
def ag_strip_plot(antigen=df_luminex_plot.Ag.values[0], ig_or_fcr="Ig Titer"):
    
    if ig_or_fcr == "Ig Titer":
        df_small = df_luminex_plot.loc[(df_luminex_plot.Ag == antigen) & (df_luminex_plot.Ig_FcR.str.contains("Ig")), :]
        title = f"Titers of {antigen}-specific Antibodies"
    else:
        df_small = df_luminex_plot.loc[(df_luminex_plot.Ag == antigen) & (df_luminex_plot.Ig_FcR.str.contains("Fc")), :]
        title = f"Binding of {antigen}-specific Antibodies to FcRs"
        
    strip = hv.Scatter(
                data=df_small,
                kdims=['Ig_FcR'],
                vdims=["value", "Treatment", "Sample"],
            ).opts(
                color='Treatment',
                title=title,
                width=700,
                height=500,
                fontsize={'labels': 11, 'xticks': 10, 'yticks': 10}
            )
    
    return strip


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
                width=900,
                height=500,
                fontsize={'labels': 11, 'xticks': 10, 'yticks': 10}
            )
    
    return strip

dash2 = pn.Column(ab_fcr_select, pn.Spacer(height=30), ig_fcr_strip_plot)


@pn.depends(assay_select.param.value, groups_toggle.param.value)
def func_strip_boxplot(func_assay=df_func_plot.variable.unique()[0],
                       combine_imm_groups=False):
    
    if combine_imm_groups == False:
    
        strip = hv.Scatter(
            data=df_func_plot.loc[df_func_plot.variable == func_assay],
            kdims=['Treatment'],
            vdims=['value', 'Sample'],
        ).opts(
            color='Treatment',
            title=f"{func_assay} Challenged Hamsters Separated",
            width=700,
            height=500,
        )

        box = hv.BoxWhisker(
            data=df_func_plot.loc[df_func_plot.variable == func_assay],
            kdims=['Treatment'],
            vdims=['value'],
        ).opts(
            box_fill_color='lightgray',
            outlier_alpha=0,
        )
    
    else:
    
        strip = hv.Scatter(
            data=df_func_plot.loc[df_func_plot.variable == func_assay],
            kdims=['Group'],
            vdims=['value', 'Sample', 'Treatment'],
        ).opts(
            color='Treatment',
            title=f"{func_assay} Challenged Hamsters Combined",
            width=700,
            height=500,
        )

        box = hv.BoxWhisker(
            data=df_func_plot.loc[df_func_plot.variable == func_assay],
            kdims=['Group'],
            vdims=['value'],
        ).opts(
            box_fill_color='lightgray',
            outlier_alpha=0,
        )

    return box * strip


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


tab2 = pn.Row(pn.layout.HSpacer(),
              pn.Column(pn.layout.VSpacer(), assay_select, pn.Spacer(height=50), groups_toggle, pn.layout.VSpacer()), 
              pn.Spacer(width=40),
              func_strip_boxplot,
              pn.layout.HSpacer())

tab3 = pn.Row(pn.layout.HSpacer(), dash3, pn.layout.HSpacer())

dashboard = pn.Tabs(('Luminex Plots', tab1), ("Functional Assays", tab2), ("Z-Score Heatmap", tab3))

dashboard.servable()
