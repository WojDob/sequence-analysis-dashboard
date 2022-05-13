from io import StringIO
from dash import Dash, html, dcc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from Bio import SeqIO
import pandas as pd
import plotly.express as px


app = Dash(__name__)


@app.callback(
    [Output("textarea-state-output", "children"), Output("chart", "figure")],
    Input("input-button", "n_clicks"),
    State("textarea-state", "value"),
)
def update_output(n_clicks, value):
    if n_clicks > 0 and value:
        record = string_to_seq(value)
        stats = get_sequence_statistics(record)
        fig = get_sequence_graph_figure(record)
        return stats, fig
    raise PreventUpdate


def string_to_seq(fasta_string):
    fasta_io = StringIO(fasta_string)
    record = next(SeqIO.parse(fasta_io, "fasta"))
    fasta_io.close()
    return record


def get_sequence_statistics(record):
    sequence = record.seq
    output = f"""
        id: {record.id}
        description {record.description}
        length: {len(sequence)}
    """
    if is_nucleotide(sequence):
        try:
            complementary = sequence.reverse_complement()
        except ValueError as e:
            output += f"ERROR - {str(e)}"
        else:
            output += f"complementary sequence: {complementary}"
    return output


def is_nucleotide(sequence):
    return not sequence.upper().strip("ATGCNU")


def get_sequence_graph_figure(record):
    symbols, counts = count_symbols(record.seq)
    df = pd.DataFrame(
        {
            "Symbol": symbols,
            "Count": counts,
        }
    )
    fig = px.bar(df, x="Symbol", y="Count", color="Count")
    return fig


def count_symbols(sequence):
    contents = {}
    for symbol in sequence:
        if symbol in contents:
            contents[symbol] += 1
        else:
            contents[symbol] = 1
    symbols = []
    counts = []
    for key, value in contents.items():
        symbols.append(key)
        counts.append(value)
    return symbols, counts


app.layout = html.Div(
    [
        dcc.Textarea(
            id="textarea-state",
            value=">example_1\nVLSISYSRSESSLETIGQRKPSTFSWSSRAASRSSWERGP",
            style={"width": "100%", "height": 200},
        ),
        html.Button("Submit", id="input-button", n_clicks=0),
        html.Div(id="textarea-state-output", style={"whiteSpace": "pre-line"}),
        html.Div(dcc.Graph(id="chart")),
    ]
)


if __name__ == "__main__":
    app.run_server(debug=True)
