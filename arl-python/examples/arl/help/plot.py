import pandas as pd
import numpy as np

import plotly.offline as py
import plotly.graph_objs as go

# plot table and bar chart
def plot_table_bar(time_array, columns, index, **titles):
    fig_title = titles['fig_title']
    table_title = titles['table_title']
    bar_title = titles['bar_title']

    # change numpy array to pandas dataFrame
    df = pd.DataFrame(time_array, columns=columns, index=index)

    # table trace
    trace_table = go.Table(
        domain=dict(x= [0.0, 0.5], y= [0, 1.0]),
        columnwidth = [1.4] + [1] * len(columns),

        header=dict(
            values = [''] + ['<b>'+col for col in columns],
            fill = dict(color='rgb(54,56,128)'),
            align = ['center'],
            font = dict(color='white', size=13),
            height = 35,
        ),

        cells=dict(
            values = [['<b>'+idx for idx in list(index)]]
                    + [df.iloc[:,j] for j in range(df.shape[1])],
            fill = dict(color='#F5F8FF'),
            align = ['right'],
            font = dict(size=13),
            height = 35,
        )
    )

    # bar traces
    trace_b1 = go.Bar(
        x = columns, y = df.loc[index[0]], name = index[0],
        marker=dict(color='rgb(54,141,202)'), width=0.35,
    )
    trace_b2 = go.Bar(
        x = columns, y = df.loc[index[1]], name = index[1],
        marker=dict(color='rgb(104,85,201)'), width=0.35,
    )
    trace_bar = [trace_b1, trace_b2]

    # layout
    axis=dict(showline=False, zeroline=True, mirror=False, 
        ticklen=4, tickfont=dict(size=11)
    )
    title = dict(showarrow=False, font=dict(size=15),
        xref='paper', yref='paper', y=1.01,  
        xanchor='left', yanchor='bottom',
    )
    layout = dict(width=950, height=400, autosize=False, 
        title='<b>' + fig_title,
        margin = dict(t=100,l=0,r=0,b=100),
        xaxis=dict(axis, **dict(domain=[0.62,0.96])),
        yaxis=dict(axis, **dict(domain=[0.3,1], title='Time(s)', titlefont=dict(size=12))),
        annotations=
          [dict(title, **dict(text='<b>' + table_title, x=0.2)), 
           dict(title, **dict(text='<b>' + bar_title, x=0.73))]
    )

    # plot table and figure
    fig = dict(data=[trace_table] + trace_bar, layout=layout)
    py.iplot(fig)
    
    
# plot table and bar chart by some fix parameters
def show(origin, optimized, columns, fig_title):
    speedup = ['x{:.1f}'.format(s) for s in origin / optimized]
    summerize = np.array([origin, optimized, speedup])
    index = ['Origin', 'Optimized', 'Speedup']
        
    table_title = 'Running Time(s)'
    bar_title = 'Compare'

    plot_table_bar(summerize, columns, index, 
               fig_title=fig_title, table_title=table_title, bar_title=bar_title)