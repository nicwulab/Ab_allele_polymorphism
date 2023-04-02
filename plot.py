import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from sinaplot import sinaplot


def plot_histogram(df, column, graph_title, xlabel, ylabel):
    df['floatddg'] = df[column].astype('float64')

    q = df['floatddg'].quantile(0.99)
    qdf = df[df['floatddg'] < q]
    q_low = df['floatddg'].quantile(0.01)
    q_hi = df['floatddg'].quantile(0.99)

    df_filtered = df[(df["floatddg"] < q_hi) & (df["floatddg"] > q_low)]
    list_filtered = df_filtered["floatddg"].tolist()
    sns.histplot(list_filtered, stat='count')
    plt.title(graph_title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.savefig('ddgs_histogram_0.99.png')
    plt.clf()

    l = df[column].tolist()
    intl = [float(i) for i in l if 'nan' not in i]

    sns.histplot(intl, stat='count')
    plt.title(graph_title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig('ddgs_histogram.png')
    plt.clf()

    # DDG cap = 20 and floor = -10
    l = df[column].tolist()
    intl = [float(i) for i in l if 'nan' not in i]
    for i in range(len(intl)):
        if intl[i] > 20:
            intl[i] = 20
        if intl[i] < -10:
            intl[i] = -10
    sns.histplot(intl, stat='count')
    plt.title(graph_title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig('ddgs_histogram20.png')
    plt.clf()

def plot_ddg_ddgs_antibody_scatter_plot(baseline_freq_ddgs_df, column_name):
    baseline_freq_ddgs_df['DDG'] = baseline_freq_ddgs_df['DDG'].astype('float64')
    baseline_freq_ddgs_df['DDG_antibody'] = baseline_freq_ddgs_df['DDG_antibody'].astype('float64')

    sns.scatterplot(data=baseline_freq_ddgs_df, x="DDG", y="DDG_antibody")
    plt.savefig(column_name+'DDG_DDG_antibody_scatterplot.png')
    plt.clf()

    baseline_freq_ddgs_df.loc[baseline_freq_ddgs_df['DDG_antibody'] < -10, 'DDG_antibody'] = -10
    baseline_freq_ddgs_df.loc[baseline_freq_ddgs_df['DDG_antibody'] > 20, 'DDG_antibody'] = 20

    baseline_freq_ddgs_df.loc[baseline_freq_ddgs_df['DDG'] < -10, 'DDG'] = -10
    baseline_freq_ddgs_df.loc[baseline_freq_ddgs_df['DDG'] > 20, 'DDG'] = 20
    fig, ax = plt.subplots()

    ax.axline((0, 0), slope=1)
    sns.scatterplot(data=baseline_freq_ddgs_df, x="DDG", y="DDG_antibody")
    plt.savefig(column_name+'DDG_DDG_antibody_scatterplot_20_yx.png')
    plt.show()
    plt.clf()

    r = baseline_freq_ddgs_df.corr(method="spearman")
    print(r)

def plot_baseline_freq_ddgs_antibody_scatter_plot(baseline_freq_ddgs_df, column_name):
    baseline_freq_ddgs_df['baseline_freq_original'] = baseline_freq_ddgs_df['baseline_freq_original'].astype('float64')
    baseline_freq_ddgs_df['DDG_antibody'] = baseline_freq_ddgs_df['DDG_antibody'].astype('float64')

    sns.scatterplot(data=baseline_freq_ddgs_df, x="baseline_freq_original", y="DDG_antibody")
    plt.savefig(column_name+'baseline_freq_DDG_antibody_scatterplot.png')
    plt.clf()

    baseline_freq_ddgs_df.loc[baseline_freq_ddgs_df['DDG_antibody'] < -10, 'DDG_antibody'] = -10
    baseline_freq_ddgs_df.loc[baseline_freq_ddgs_df['DDG_antibody'] > 20, 'DDG_antibody'] = 20
    sns.scatterplot(data=baseline_freq_ddgs_df, x="baseline_freq_original", y="DDG_antibody")
    plt.savefig(column_name+'baseline_freq_DDG_antibody_scatterplot_20.png')
    plt.clf()

    r = baseline_freq_ddgs_df.corr(method="spearman")
    print(r)


def plot_baseline_freq_ddgs_scatter_plot(baseline_freq_ddgs_df, column_name):
    baseline_freq_ddgs_df['baseline_freq_original'] = baseline_freq_ddgs_df['baseline_freq_original'].astype('float64')
    baseline_freq_ddgs_df['DDG'] = baseline_freq_ddgs_df['DDG'].astype('float64')

    sns.scatterplot(data=baseline_freq_ddgs_df, x="baseline_freq_original", y="DDG")
    plt.savefig(column_name+'baseline_freq_DDG_scatterplot.png')
    plt.clf()

    baseline_freq_ddgs_df.loc[baseline_freq_ddgs_df['DDG'] < -10, 'DDG'] = -10
    baseline_freq_ddgs_df.loc[baseline_freq_ddgs_df['DDG'] > 20, 'DDG'] = 20
    sns.scatterplot(data=baseline_freq_ddgs_df, x="baseline_freq_original", y="DDG")
    plt.savefig(column_name+'baseline_freq_DDG_scatterplot_20.png')
    plt.clf()

    r = baseline_freq_ddgs_df.corr(method="spearman")
    print(r)

def plot_baseline_freq_species_sinaplot(baseline_freq_species):

    plt.gcf().set_size_inches(25, 25)

    baseline_freq_species['baseline_freq_original'] = baseline_freq_species['baseline_freq_original'].astype('float64')
    ddgs_species = baseline_freq_species.dropna()

    sinaplot(x='baseline_freq_original', y='antigen_species_id', data=ddgs_species)
    plt.savefig('baseline_freq_species_sinaplot.png')
    plt.clf()


def plot_ddgs_species_sinaplot(ddgs_species):

    plt.gcf().set_size_inches(25,25)

    ddgs_species['DDG'] = ddgs_species['DDG'].astype('float64')
    ddgs_species = ddgs_species.dropna()

    list_filtered = ddgs_species['DDG'].tolist()

    sinaplot(x='DDG', y='antigen_species_id', data=ddgs_species)
    plt.savefig('ddgs_species_sinaplot.png')
    plt.clf()

    q_low = ddgs_species['DDG'].quantile(0.01)
    q_hi = ddgs_species['DDG'].quantile(0.99)

    df_filtered = ddgs_species[(ddgs_species["DDG"] < q_hi) & (ddgs_species["DDG"] > q_low)]
    list_filtered = df_filtered['DDG'].tolist()

    sinaplot(x='DDG', y='antigen_species_id', data=df_filtered)
    plt.savefig('ddgs_species_sinaplot0.99.png')
    plt.clf()

    # DDG cap = 20 and floor = -10
    ddgs_species.loc[ddgs_species['DDG'] < -10, 'DDG'] = -10
    ddgs_species.loc[ddgs_species['DDG'] > 20, 'DDG'] = 20

    sinaplot(x='DDG', y='antigen_species_id', data=ddgs_species)
    plt.savefig('ddgs_species_sinaplot20.png')
    plt.clf()

def plot_baseline_freq_ddgs_ddgs_antibody_scatter_plot(baseline_freq_ddgs_df, column_name):
    baseline_freq_ddgs_df['baseline_freq_original'] = baseline_freq_ddgs_df['baseline_freq_original'].astype('float64')
    baseline_freq_ddgs_df['DDG'] = baseline_freq_ddgs_df['DDG'].astype('float64')
    baseline_freq_ddgs_df['DDG_antibody'] = baseline_freq_ddgs_df['DDG_antibody'].astype('float64')
    baseline_freq_ddgs_df['ddg-ddf_antibody'] = baseline_freq_ddgs_df['DDG'] - baseline_freq_ddgs_df['DDG_antibody']
    print(baseline_freq_ddgs_df)

    sns.scatterplot(data=baseline_freq_ddgs_df, x="baseline_freq_original", y="ddg-ddf_antibody")
    plt.savefig(column_name+'baseline_freq_DDG-DDG_antibody_scatterplot.png')
    plt.clf()

    r = baseline_freq_ddgs_df.corr(method="spearman")
    print(r)

def main():

    epitope_identification_with_baseline_variation_df = pd.read_csv(
        'results/epitope_identification_with_antibody_only_ddgs.csv')
    ddgs = epitope_identification_with_baseline_variation_df[['DDG']]
    baseline_freq_ddgs = epitope_identification_with_baseline_variation_df[['baseline_freq_original', 'DDG']]
    baseline_freq_antibody_ddgs = epitope_identification_with_baseline_variation_df[['baseline_freq_original', 'DDG_antibody']]
    baseline_freq_antibody_ddgs_ddgs = epitope_identification_with_baseline_variation_df[
        ['baseline_freq_original', 'DDG_antibody', 'DDG']]
    ddgs_species = epitope_identification_with_baseline_variation_df[['DDG', 'antigen_species_id']]
    baseline_freq_species = epitope_identification_with_baseline_variation_df[['baseline_freq_original', 'antigen_species_id']]
    ddg_antibody_ddgs = epitope_identification_with_baseline_variation_df[
        ['DDG', 'DDG_antibody']]

    plot_histogram(ddgs, 'DDG', 'Distribution of DDG', 'DDG value', 'Count')
    plot_baseline_freq_ddgs_scatter_plot(baseline_freq_ddgs, '')

    plot_baseline_freq_ddgs_antibody_scatter_plot(baseline_freq_antibody_ddgs, '')
    plot_ddg_ddgs_antibody_scatter_plot(ddg_antibody_ddgs,'')
    plot_baseline_freq_ddgs_ddgs_antibody_scatter_plot(baseline_freq_antibody_ddgs_ddgs, '')

    HIV_df = epitope_identification_with_baseline_variation_df[epitope_identification_with_baseline_variation_df['antigen_species_id'] == 'HIV']
    HIV_df = HIV_df[['baseline_freq_original', 'DDG']]
    plot_baseline_freq_ddgs_scatter_plot(HIV_df, 'HIV')

    influenza_df = epitope_identification_with_baseline_variation_df[epitope_identification_with_baseline_variation_df['antigen_species_id'] == 'influenza']
    influenza_df = influenza_df[['baseline_freq_original', 'DDG']]
    plot_baseline_freq_ddgs_scatter_plot(influenza_df, 'influenza')

    plot_ddgs_species_sinaplot(ddgs_species)
    plot_baseline_freq_species_sinaplot(baseline_freq_species)


    dssp_area = pd.read_csv('results/add_dssp_area.csv')
    dssp_area['DDG'] = dssp_area['DDG'].astype('float64')
    dssp_area['DDG_antibody'] = dssp_area['DDG_antibody'].astype('float64')
    sns.set(rc={'font.family': 'Arial', 'font.size': 7})
    sns.set_style("whitegrid")  # Set the style to whitegrid
    sns.set_style("whitegrid", {'axes.grid': False})

    p = sns.color_palette("flare", as_cmap=True)
    fig, ax = plt.subplots()

    a4_dims = (2.5, 2.5)
    fig, ax = plt.subplots(figsize=a4_dims)
    ax.set(xlabel='ΔΔG (Ab-Ag complex)', ylabel='ΔΔG (Ab alone)')
    g = sns.scatterplot(ax=ax,data=dssp_area, x="DDG", y="DDG_antibody", hue="dssp_asa", palette=p)
    ax.axline((0, 0), slope=1, color='lightgrey')
    ax.spines[['right', 'top']].set_visible(False)

    box = g.get_position()
    g.set_position([box.x0, box.y0, box.width * 0.85, box.height])  # resize position
    # Put a legend to the right side
    g.legend(loc='center right', bbox_to_anchor=(1.65, 0.5), ncol=1, title='dssp_asa')
    # plt.tight_layout()
    plt.savefig('DDG_DDG_antibody_scatterplot_yx_hue.png', bbox_inches='tight')
    plt.clf()

    dssp_area.loc[dssp_area['DDG_antibody'] < -10, 'DDG_antibody'] = -10
    dssp_area.loc[dssp_area['DDG_antibody'] > 20, 'DDG_antibody'] = 20

    dssp_area.loc[dssp_area['DDG'] < -10, 'DDG'] = -10
    dssp_area.loc[dssp_area['DDG'] > 20, 'DDG'] = 20

    p = sns.color_palette("flare", as_cmap=True)
    fig, ax = plt.subplots()

    a4_dims = (2.5, 2.5)
    fig, ax = plt.subplots(figsize=a4_dims)
    ax.set(xlabel='ΔΔG (Ab-Ag complex)', ylabel='ΔΔG (Ab alone)')
    g = sns.scatterplot(data=dssp_area, x="DDG", y="DDG_antibody", hue="dssp_asa", palette=p, edgecolor="black", alpha=0.5, s=20)

    ax.axline((0, 0), slope=1, color='#A9A9A9', linestyle='--')
    ax.spines[['right', 'top']].set_visible(False)
    ax.spines['bottom'].set_color('#000000')
    ax.spines['left'].set_color('#000000')
    box = g.get_position()
    g.set_position([box.x0, box.y0, box.width * 0.85, box.height])  # resize position
    # Put a legend to the right side
    ax.get_legend().remove()

    # g.legend(loc='center right', bbox_to_anchor=(1.65, 0.5), ncol=1, title='dssp_asa')
    g.tick_params(left=True)
    g.tick_params(bottom=True)
    plt.xticks([-10, 0, 10, 20])
###
    norm = plt.Normalize(dssp_area['dssp_asa'].min(), dssp_area['dssp_asa'].max())
    sm = plt.cm.ScalarMappable(cmap=p, norm=norm)
    sm.set_array([])
###
    cax = fig.add_axes([ax.get_position().x1 + 0.05, ax.get_position().y0, 0.06, ax.get_position().height * 0.4])
    clb = ax.figure.colorbar(sm, cax=cax)
    clb.ax.set_title('RSA', fontsize=10, loc='left')
    clb.outline.set_edgecolor('#000000')
    clb.outline.set_linewidth(0.5)

    # plt.tight_layout()
    plt.savefig('DDG_DDG_antibody_scatterplot_20_yx_hue.png', bbox_inches='tight', dpi=600)
    plt.clf()

if __name__ == "__main__":
    main()
