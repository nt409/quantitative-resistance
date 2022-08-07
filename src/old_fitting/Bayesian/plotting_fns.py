from math import ceil, exp, log
import numpy as np
import itertools
from scipy.stats import norm, lognorm
import pandas as pd

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px


standard_layout = go.Layout(
    font=dict(size=22),
    template="plotly_white",
    width=1600,
    height=800,
    # showlegend = True,
)


def plot_output(output):

    # N = len(output['I0s'])

    I0s = output['I0s']
    betas = output['betas']

    trace = go.Scatter(
        x=I0s,
        y=betas,
        mode='markers',
        marker={'opacity': 0.5},
        marker_color=[
            f'rgb(100,0,{255-i/len(betas)})' for i in range(len(betas)-1)],

        opacity=0.2
    )

    fig = go.Figure(data=trace, layout=standard_layout)
    fig.update_layout(showlegend=False)

    fig.update_xaxes(title='I0')
    fig.update_yaxes(title='beta')

    return fig


def plot_sevs(data, gen_use):

    traces = []

    for yy, colour, loc in zip([data, gen_use], ['red', 'blue'], [0, 1]):

        xx = np.random.normal(loc, 0.1, size=len(yy))

        trace = go.Scatter(
            x=xx,
            y=yy,
            mode='markers',
            marker={'opacity': 0.5},
            marker_color=colour,
            opacity=0.2
        )

        traces.append(trace)

        mean = go.Scatter(
            x=[loc],
            y=[np.mean(yy)],
            mode='markers',
            marker={'opacity': 1},
            marker_color='black',
            marker_size=12,
            marker_symbol='x',
            opacity=1
        )

        traces.append(mean)

    fig = go.Figure(data=traces, layout=standard_layout)
    fig.update_layout(showlegend=False)

    fig.update_xaxes(title='<i>I<sub>0</sub></i>', rangemode='tozero')
    fig.update_yaxes(title='<i>'+u'\u03B2'+'</i>', rangemode='tozero')

    return fig


def contour_plot(z, minI0, maxI0, minBeta, maxBeta):

    XNum = z.shape[0]
    YNum = z.shape[1]

    # give points at the edge of each square
    xx = np.linspace(minI0, maxI0, XNum+1)
    yy = np.linspace(minBeta, maxBeta, YNum+1)

    trace = go.Heatmap(
        z=z,
        x=xx,
        y=yy,
        colorbar=dict(nticks=10, ticks='outside',
                      ticklen=5, tickwidth=1,
                      showticklabels=True,
                      tickangle=0, tickfont_size=12)
    )

    fig = go.Figure(data=trace, layout=standard_layout)

    # range=[0, 8*10**(-3)])
    fig.update_xaxes(title='<i>I<sub>0</sub></i>', rangemode='tozero')
    # range=[0, 1.3*10**(-2)])
    fig.update_yaxes(title='<i>'+u'\u03B2'+'</i>', rangemode='tozero')

    return fig


def likelihood_contour(z, X, Y):

    # plot_log_likelihood = False

    # if plot_log_likelihood:
    #     for i,j in itertools.product(range(z.shape[0]),range(z.shape[1])):
    #         z[i,j] = log(z[i,j])

    trace = go.Contour(
        z=z,
        x=X,
        y=Y,
        #    contours_coloring='heatmap',
        colorbar=dict(nticks=10, ticks='outside',
                      ticklen=5, tickwidth=1,
                      showticklabels=True,
                      tickangle=0, tickfont_size=12)
    )

    fig = go.Figure(data=trace, layout=standard_layout)

    fig.update_xaxes(title='log(<i>I<sub>0</sub></i>)')
    fig.update_yaxes(title='<i>'+u'\u03B2'+'</i>', rangemode='tozero')

    return fig


def histogram(z, min_I, max_I, min_B, max_B):

    fig = make_subplots(rows=1, cols=2)

    # if plot=='I0':

    x = np.linspace(min_I, max_I, z.shape[0])
    y = [sum(z[i, :]) for i in range(z.shape[0])]

    trace = go.Scatter(
        x=x,
        y=y,
        mode='lines',
    )

    fig.add_trace(trace, row=1, col=1)

    x = np.linspace(min_B, max_B, z.shape[1])
    y = [sum(z[:, i]) for i in range(z.shape[1])]

    trace = go.Scatter(
        x=x,
        y=y,
        mode='lines',
    )

    fig.add_trace(trace, row=1, col=2)

    fig.update_layout(standard_layout)

    fig.update_xaxes(title='<i>I<sub>0</sub></i>',
                     rangemode='tozero', row=1, col=1)
    fig.update_xaxes(title='<i>'+u'\u03B2'+'</i>',
                     rangemode='tozero', row=1, col=2)

    fig.update_yaxes(title='Count', row=1, col=1)
    fig.update_yaxes(title='Count', row=1, col=2)

    return fig


def plot_sampled_params(logI0s, betas):

    line = go.Scatter(x=logI0s,
                      y=betas,
                      mode="markers",
                      marker={'opacity': 0.5},

                      )

    traces = [line]
    fig = go.Figure(data=traces, layout=standard_layout)
    fig.update_xaxes(title='log(I0)')
    fig.update_yaxes(title='beta')

    return fig


def prior_hist(priors, plot_I0):

    betas = priors.betas
    L_I0s = priors.logI0s

    df = pd.DataFrame(dict(L_I0s=L_I0s, betas=betas))  # , logI0s=logI0s))

    if plot_I0:
        fig = px.histogram(df, x="L_I0s")
        fig.update_xaxes(title='log(<i>I<sub>0</sub></i>)')
    else:
        fig = px.histogram(df, x="betas")
        fig.update_xaxes(title='<i>'+u'\u03B2'+'</i>')

    # print(fig['layout'])
    # exit()

    traces = []

    for vec, dist_name, bool in zip([L_I0s, betas],
                                    ['L_I0s', 'beta'],
                                    [plot_I0, not plot_I0]):

        if bool:

            if dist_name == 'L_I0s':
                xmin = -20  # min(vec)
                xmax = -5  # max(vec)

                x = np.linspace(xmin, xmax, 100)

                y = priors.L_I_dist.pdf(x, *priors.L_I0_pars)

                count = len(vec)
                Nbars = 3
                height = count/Nbars
                width = (xmax-xmin)
                scale = width*height
                y = np.asarray(y)*scale

            else:
                xmin = 0.002  # min(vec)
                xmax = 0.012  # max(vec)

                x = np.linspace(xmin, xmax, 100)

                y = priors.B_dist.pdf(x, *priors.Bpars)

                count = len(vec)
                Nbars = 5
                height = count/Nbars
                width = (xmax-xmin)

                scale = width*height
                y = np.asarray(y)*scale

            line = go.Scatter(
                x=x,
                y=y,
                mode='lines',
                name='Prior (scaled)',
                line=dict(color='black')
            )

            traces.append(line)

            fig.add_trace(line)

    fig.update_layout(standard_layout)
    fig.update_layout(legend=dict(x=0.05, y=1, xanchor='left', yanchor='top'))

    return fig

    # fig = go.Figure(data=traces, layout=standard_layout)
