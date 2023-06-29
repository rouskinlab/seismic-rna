from .util import *

import plotly.graph_objects as go

from ..draw import util
from itertools import cycle
from sklearn import metrics
from sklearn.linear_model import LogisticRegression
from plotly.subplots import make_subplots
import plotly.express as px
from dms_ci import dms_ci

LIST_COLORS = ['red', 'green', 'blue', 'orange', 'purple', 'black', 'yellow', 'pink', 'brown', 'grey', 'cyan',
               'magenta']


def mutation_fraction(df, show_ci: bool = True) -> dict:
    assert len(df) == 1, "df must have only one row"
    mh = df.iloc[0]
    cmap = {"A": "red", "T": "green", "G": "orange", "C": "blue"}  # Color map

    traces, layouts = [], []
    mh['index_selected'] = [i + 1 for i in
                            range(len(mh['sequence']))]  # TODO[i + 1 for i in mh.index_selected] # index starts at 1
    mh_unrolled = pd.DataFrame(
        {'mut_rate': list(mh.sub_rate), 'base': list(mh.sequence), 'index_reset': list(range(len(mh.index_selected))),
         'index_selected': mh.index_selected, 'paired': list(mh.structure)})

    for bt in set(mh['sequence']):
        df_loc = mh_unrolled[mh_unrolled['base'] == bt]
        if len(df_loc) == 0:
            continue

        hover_attr = pd.DataFrame({'mut_rate': list(df_loc.mut_rate),
                                   'base': list(df_loc.base),
                                   'index': df_loc['index_selected'],
                                   'paired': [{'.': False, '(': True, ')': True}[s] for s in df_loc.paired]})
        traces.append(go.Bar(
            x=np.array(df_loc['index_reset']),
            y=np.array(df_loc['mut_rate']),
            name=bt,
            marker_color=cmap[bt],
            text=hover_attr,
            hovertemplate=''.join(["<b>" + ha + ": %{text[" + str(i) + "]}<br>" for i, ha in enumerate(hover_attr)]),
        ))
        if show_ci:
            traces[-1].update(
                error_y=dict(
                    type='data',
                    symmetric=False,
                ))

    fig = go.Figure(data=traces)

    fig.update_layout(
        title=f"{mh['sample']} - {mh['reference']} - {mh['section']} - {mh['cluster']} - {mh['num_aligned']} reads",
        xaxis=dict(title="Position"),
        yaxis=dict(title="Mutation fraction", range=[0, 0.1]))

    fig.update_yaxes(
        gridcolor='lightgray',
        linewidth=1,
        linecolor='black',
        mirror=True,
        autorange=True
    )
    fig.update_xaxes(
        linewidth=1,
        linecolor='black',
        mirror=True,
        autorange=True
    )

    fig.update_xaxes(
        tickvals=mh_unrolled['index_reset'],
        ticktext=["%s %s" % ({'.': '(U)', '(': '(P)', ')': '(P)'}[x], str(y)) for (x, y) in
                  zip(mh['structure'], mh['index_selected'])],
        tickangle=90,
        autorange=True
    )

    # make the background white 
    fig.update_layout(plot_bgcolor='white', paper_bgcolor='white')

    return {'fig': fig, 'df': mh}


def mutation_fraction_identity(data, show_ci: bool = True) -> dict:
    assert len(data) > 0, "The combination of sample, reference and section does not exist in the dataframe"
    assert len(data) == 1, "The combination of sample, reference and section is not unique in the dataframe"
    data = data.iloc[0]

    df = pd.DataFrame(index=list(data['sequence']))
    fig = go.Figure()
    color_map = {'A': 'red', 'C': 'blue', 'G': 'yellow', 'T': 'green'}

    data['err_min'] = [dms_ci(p, data['num_aligned'])[0] for p in data['sub_rate']]
    data['err_max'] = [dms_ci(p, data['num_aligned'])[1] for p in data['sub_rate']]

    for base in ['A', 'C', 'G', 'T']:
        df[base] = np.array(data['sub_' + base]) / np.array(data['info'])
        fig.add_trace(
            go.Bar(x=np.arange(len(data['sequence'])), y=list(df[base]), marker_color=color_map[base], showlegend=True,
                   name=base))

    # add error bars to stacked_bar[-1]
    if show_ci:
        fig.data[-1].error_y = dict(
            type='data',
            array=[data['err_max'][i] - data['sub_rate'][i] for i in range(len(data['sequence']))],
            arrayminus=[data['sub_rate'][i] - data['err_min'][i] for i in range(len(data['sequence']))],
            visible=True,
            symmetric=False,
            thickness=1.5,
            width=2,
            color='black'
        )

    fig.update_layout(
        title='Mutation fraction identity - {} - {} - {} - {} - {} reads'.format(data['sample'], data['reference'],
                                                                                 data['section'], data['cluster'],
                                                                                 data['num_aligned']),
        xaxis_title='Position',
        yaxis_title='Mutation fraction',
    )

    fig.update_xaxes(tickangle=0,
                     tickvals=np.arange(len(df.index)), ticktext=list(df.index), tickfont={'size': 8})

    fig.update_layout(barmode='stack')
    fig.update_layout(plot_bgcolor='white', paper_bgcolor='white')

    return {'fig': fig, 'df': df}


def mutations_in_barcodes(data):
    fig = go.Figure()

    len_barcode = len(data['sequence'].iloc[2])
    for sample in data['sample'].unique():
        hist = np.sum(np.stack(data[data['sample'] == sample]['sub_hist'].values), axis=0)
        bin_edges = np.arange(0, len_barcode, 1)

        fig.add_trace(
            go.Bar(
                x=bin_edges[:-1],
                y=hist,
                name=sample,
                visible=False,
                hovertemplate='Number of mutations: %{x}<br>Number of reads: %{y}<extra></extra>',
            ))

    fig.data[0].visible = True

    fig.update_layout(barmode='stack', title='Number of mutations in barcodes - {}'.format(data['sample'].unique()[0]))
    fig.update_layout(
        updatemenus=[
            dict(
                active=0,
                buttons=list([
                    dict(label=sample,
                         method="update",
                         args=[{"visible": [sample == s for s in data['sample'].unique()]},
                               {"title": 'Number of mutations in barcodes - {}'.format(sample)}])
                    for sample in data['sample'].unique()
                ]),
                x=0.7,
                xanchor='left',
                y=1.1,
                yanchor='top'
            )
        ])

    fig.update_layout(plot_bgcolor='white', paper_bgcolor='white')

    return {'fig': fig, 'data': data[['sample', 'reference', 'sub_hist']]}


def deltaG_vs_sub_rate(df: pd.DataFrame, models: list[str] = [], savefile=None, auto_open=False, use_iplot=True,
                       title=None) -> dict:
    df_temp = pd.DataFrame()
    for _, row in df.iterrows():
        df_temp = pd.concat([df_temp, pd.DataFrame(
            {'reference': row.reference, 'index': row.index_selected, 'sub_rate': row.sub_rate,
             'num_aligned': row.num_aligned, 'deltaG': row['deltaG'], 'base': list(row.sequence),
             'paired': [s != '.' for s in row.structure]}, index=[i + 1 for i in range(len(row.sequence))])])

    assert len(df_temp) > 0, "No data to plot"
    df = df_temp.reset_index()

    hover_attr = ['num_aligned', 'sub_rate', 'base', 'index', 'reference', 'deltaG']
    tra = {}
    for is_paired, prefix in zip([True, False], ['Paired ', 'Unpaired ']):
        markers = cycle(list(range(153)))
        x = np.array(df[df.paired == is_paired]['deltaG'])
        y = np.array(df[df.paired == is_paired]['sub_rate'])
        tra[prefix] = go.Scatter(
            x=x,
            y=y,
            text=df[df.paired == is_paired][hover_attr],
            marker_size=df[df.paired == is_paired]['num_aligned'] / 200,
            mode='markers',
            name=prefix,
            hovertemplate=''.join(["<b>" + ha + ": %{text[" + str(i) + "]}<br>" for i, ha in enumerate(hover_attr)]),
            line=dict(color='green' if is_paired else 'red'))

        for m in models:
            if len(y) > 0:
                fit = util.Fit()
                x_sorted, pred_y_sorted = fit.predict(x, y, m, prefix)
                tra[fit.get_legend()] = go.Scatter(
                    x=x_sorted,
                    y=pred_y_sorted,
                    mode='lines+markers',
                    name=fit.get_legend(),
                    marker=dict(symbol=next(markers)),
                    line=dict(color='darkseagreen' if is_paired else 'crimson', dash='dash'))

    layout = dict(title='Mutation rates of paired / unpaired residues vs the expected energy of the molecule',
                  xaxis=dict(title='DeltaG', ticklen=5, zeroline=False),
                  yaxis=dict(title='Mutation fraction ', ticklen=5, zeroline=False),
                  )

    fig = go.Figure(data=list(tra.values()), layout=layout)
    fig.update_layout(plot_bgcolor='white', paper_bgcolor='white')

    return {'fig': fig, 'df': df}


def experimental_variable_across_samples(data: pd.DataFrame, experimental_variable: str,
                                         models: List[str] = []) -> dict:
    fig = go.Figure()

    assert len(data) > 0, "No data to plot"
    assert experimental_variable in data.columns, "Experimental variable not found in data"
    assert len(data[
                   'sequence'].unique()) == 1, "More than one sequence found in data. CHeck that reference and section are unique"

    df = pd.DataFrame(
        np.hstack([np.vstack(data['sub_rate'].values), np.vstack(data[experimental_variable].values)]),
        index=data[experimental_variable],
        columns=[c + str(idx + 1) for c, idx in zip(data['sequence'].iloc[0], data['index_selected'].iloc[0])] + [
            experimental_variable]
    ).sort_index()

    for col in df.columns:
        if col == experimental_variable:
            continue
        fig.add_trace(
            go.Scatter(
                x=df.index,
                y=df[col],
                mode='lines+markers',
                name=col,
                text=experimental_variable,
                hovertemplate='<b>Experimental variable: %{x}<br>Mutation fraction: %{y}<extra></extra>',
            ),
        )

    fig.update_layout(
        title='Mutation rates across experimental variable - {}'.format(experimental_variable),
        xaxis_title=experimental_variable,
        yaxis_title='Mutation fraction',
    )
    fig.update_layout(plot_bgcolor='white', paper_bgcolor='white')

    return {'fig': fig, 'data': df}


def auc(df: pd.DataFrame, savefile=None, auto_open=False, use_iplot=True, title=None) -> dict:
    def make_roc_curve(X, y, y_pred, fig, title):
        fpr, tpr, thresholds = metrics.roc_curve(y, y_pred)
        roc_auc = metrics.auc(fpr, tpr)
        fig.add_trace(go.Scatter(x=fpr, y=tpr,
                                 mode='lines',
                                 name=title,
                                 line=dict(width=2)))
        return fig

    fig = go.Figure()
    for row in df.iterrows():
        X = row[1]['sub_rate'].reshape(-1, 1)
        y = np.array([1 if c == '.' else 0 for c in row[1]['structure']]).reshape(-1, 1)
        y_pred = LogisticRegression().fit(X, y.ravel()).predict_proba(X)[:, 1]
        fig = make_roc_curve(X, y, y_pred, fig, row[1]['unique_id'])

    fig.add_trace(
        go.Scatter(x=[0, 1], y=[0, 1], line=dict(color='black', width=2, dash='dash'), showlegend=False)
    )
    fig.update_layout(
        title='ROC Curve',
        xaxis_title='False Positive Rate',
        yaxis_title='True Positive Rate')

    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    fig.update_xaxes(constrain='domain')
    fig.update_layout(plot_bgcolor='white', paper_bgcolor='white')

    return {'fig': fig, 'df': df}


def mutation_fraction_delta(df, savefile=None, auto_open=False, use_iplot=True, title=None) -> dict:
    assert len(df) == 2, "df must have 2 row"
    mp_attr = ['sample', 'reference', 'section', 'cluster']
    df['unique_id'] = df.apply(lambda row: ' - '.join([str(row[attr]) for attr in mp_attr]), axis=1)

    mh = pd.Series(
        {
            'sub_rate': df['sub_rate'].values[0] - df['sub_rate'].values[1],
            'sequence': ''.join(
                [c1 if c1 == c2 else '-' for c1, c2 in zip(df['sequence'].values[0], df['sequence'].values[1])]),
            'title': "{} - {} reads vs {} - {} reads".format(df['unique_id'].values[0], df['num_aligned'].values[0],
                                                             df['unique_id'].values[1], df['num_aligned'].values[1])
        }
    )
    cmap = {"A": "red", "T": "green", "G": "orange", "C": "blue", '-': 'grey'}  # Color map

    traces, layouts = [], []
    mh_unrolled = pd.DataFrame(
        {'mut_rate': list(mh.sub_rate), 'base': list(mh.sequence), 'index_reset': list(range(len(mh.sequence)))})

    for bt in set(mh['sequence']):
        df_loc = mh_unrolled[mh_unrolled['base'] == bt]
        if len(df_loc) == 0:
            continue

        hover_attr = pd.DataFrame({'mut_rate': list(df_loc.mut_rate),
                                   'base': list(df_loc.base)})
        traces.append(go.Bar(
            x=np.array(df_loc['index_reset']),
            y=np.array(df_loc['mut_rate']),
            name=bt,
            marker_color=cmap[bt],
            text=hover_attr,
            hovertemplate=''.join(["<b>" + ha + ": %{text[" + str(i) + "]}<br>" for i, ha in enumerate(hover_attr)]),
        ))

    fig = go.Figure(data=traces,
                    layout=go.Layout(
                        title=go.layout.Title(text=mh['title']),
                        xaxis=dict(title="Position"),
                        yaxis=dict(title="Mutation fraction", range=[0, 0.1])))

    fig.update_yaxes(
        gridcolor='lightgray',
        linewidth=1,
        linecolor='black',
        mirror=True,
        autorange=True
    )
    fig.update_xaxes(
        linewidth=1,
        linecolor='black',
        mirror=True,
        autorange=True
    )
    fig.update_layout(plot_bgcolor='white', paper_bgcolor='white')

    return {'fig': fig, 'df': mh}


def _mutations_per_read_subplot(data):
    hist = np.sum(np.stack(data.values), axis=0)
    if (hist[0] == 0 and len(hist) < 2):
        return go.Bar(x=[], y=[], showlegend=False, marker_color='indianred')
    else:

        bin_edges = np.arange(0, np.max(np.argwhere(hist != 0)), 1)
        return go.Bar(x=bin_edges, y=hist, showlegend=False, marker_color='indianred')


def mutations_per_read_per_sample(data):
    unique_samples = data['sample'].unique()
    fig = make_subplots(rows=len(unique_samples), cols=1, vertical_spacing=0.4 / len(unique_samples),
                        subplot_titles=['Number of mutations per read - {}'.format(sample) for sample in
                                        unique_samples])
    for i_s, sample in enumerate(unique_samples):
        fig.add_trace(_mutations_per_read_subplot(data[data['sample'] == sample]['sub_hist'].reset_index(drop=True)),
                      row=i_s + 1, col=1)
        fig.update_yaxes(title='Count')
        fig.update_xaxes(dtick=10)

    fig.update_layout(autosize=True, height=len(unique_samples) * 500,
                      title='Number of mutation per read across samples')
    fig.update_layout(plot_bgcolor='white', paper_bgcolor='white')

    return {
        'fig': fig,
        'data': data
    }


def num_aligned_reads_per_reference_frequency_distribution(data):
    assert len(samples := data['sample'].unique()) == 1, "data must have 1 sample"
    data = data['num_aligned'].values
    return {
        'fig': go.Figure(
            go.Histogram(
                x=data,
                showlegend=False,
                marker_color='indianred',
                hovertemplate="Number of aligned reads: %{x}<br>Count: %{y}<extra></extra>"
            ),
            layout=go.Layout(
                title=go.layout.Title(text='{} - Reads per reference count'.format(samples[0])),
                xaxis=dict(title="Number of aligned reads"),
                yaxis=dict(title="Count"),
                plot_bgcolor='white',
                paper_bgcolor='white'
            )
        ),
        'data': data
    }


def mutation_per_read_per_reference(data):
    assert len(data) == 1, "data must have 1 row"
    data = data.iloc[0]
    sample, reference = data['sample'], data['reference']
    data = data['sub_hist']

    # normalize by the number of reads
    fig = px.bar(x=np.arange(0, len(data)), y=data)

    fig.update_layout(barmode='stack')
    fig.update_layout(title='Number of mutations per read - {} - {}'.format(sample, reference))
    fig.update_yaxes(title='Count')
    fig.update_xaxes(title='Number of mutations per read')
    fig.update_layout(plot_bgcolor='white', paper_bgcolor='white')

    return {
        'fig': fig,
        'data': data
    }


def base_coverage(data):
    assert_only_one_row(data)
    data = data.iloc[0]

    fig = go.Figure(
        go.Bar(
            x=np.arange(1, 1 + len(data['cov'])),
            y=data['cov'],
            showlegend=False,
            marker_color='indianred',
        ),
        layout=go.Layout(
            title=go.layout.Title(text='Base coverage - {} - {}'.format(data['sample'], data['reference'])),
            xaxis=dict(title="Position"),
            yaxis=dict(title="Count")
        )
    )
    fig.update_layout(plot_bgcolor='white', paper_bgcolor='white')

    return {'fig': fig, 'data': data}
