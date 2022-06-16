from io import StringIO
from dash import Dash, html, dcc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp
import pandas as pd
import plotly.express as px


app = Dash(__name__, assets_folder="../assets")


@app.callback(
    [
        Output("textarea-state-output", "children"),
        Output("chart", "figure"),
        Output("pie_chart", "figure"),
    ],
    Input("input-button", "n_clicks"),
    State("textarea-state", "value"),
)
def update_output(n_clicks, value):
    if n_clicks > 0 and value:
        record = string_to_seq(value)
        stats = get_sequence_statistics(record)
        fig = get_sequence_bar_figure(record)
        pie_chart = get_sequence_pie_figure(record)
        return stats, fig, pie_chart
    raise PreventUpdate


@app.callback(
    [
        Output(component_id="element_to_show", component_property="style"),
        Output(component_id="textarea-state-output", component_property="style"),
    ],
    Input("input-button", "n_clicks"),
)
def show_elements(n_clicks):
    if n_clicks > 0:
        return {"display": "flex"}, {"display": "block", "whiteSpace": "pre-line"}
    return {"display": "none"}, {"display": "none"}


def string_to_seq(fasta_string):
    fasta_io = StringIO(fasta_string)
    record = next(SeqIO.parse(fasta_io, "fasta"))
    fasta_io.close()
    return record


def get_sequence_statistics(record):
    sequence = record.seq
    output = f"""
id: {record.id}
description: {record.description}
length: {len(sequence)}
"""
    if is_nucleotide(sequence):
        try:
            complementary = sequence.reverse_complement()
            meltingTemp = MeltingTemp.Tm_Wallace(sequence)
            translation = sequence.translate()
        except ValueError as e:
            output += f"ERROR - {str(e)}"
        else:
            output += f"complementary sequence: {complementary}"
            output += f"\nmelting temperature of a sequence: {meltingTemp}"
            output += f"\ntranslated sequence: {translation}"
    return output


def is_nucleotide(sequence):
    return not sequence.upper().strip("ATGCNU")


def get_sequence_bar_figure(record):
    symbols, counts = count_symbols(record.seq)
    df = pd.DataFrame(
        {
            "Symbol": symbols,
            "Count": counts,
        }
    )
    fig = px.bar(df, x="Symbol", y="Count", color="Count")
    fig.update_layout(
        {
            "plot_bgcolor": "rgba(0, 0, 0, 0)",
            "paper_bgcolor": "rgba(0, 0, 0, 0)",
        }
    )
    fig.update_layout(font_color="white")
    return fig


def get_sequence_pie_figure(record):
    symbols, counts = count_symbols(record.seq)
    df = pd.DataFrame(
        {
            "Symbol": symbols,
            "Count": counts,
        }
    )
    fig = px.pie(df, values="Count", names="Symbol", hole=0.3)
    fig.update_layout(
        {
            "plot_bgcolor": "rgba(0, 0, 0, 0)",
            "paper_bgcolor": "rgba(0, 0, 0, 0)",
        }
    )
    fig.update_layout(font_color="white")
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
            className="form-control",
            value=">1-3582383\nTAGCTTATCAGACTGATGTTGAC",
            style={"width": "100%", "height": 200},
        ),
        html.Button("Submit", className="button-49", id="input-button", n_clicks=0),
        html.Div(
            id="textarea-state-output",
            className="output_field",
            style={"display": "none"},
        ),
        html.Div(
            className="graph_wrapper",
            children=[
                html.Div(dcc.Graph(id="chart")),
                html.Div(dcc.Graph(id="pie_chart")),
            ],
            style={"display": "none"},
            id="element_to_show",
        ),
    ]
)


if __name__ == "__main__":
    print("If you are running this application through docker, go to 'localhost:8050'")
    app.run_server(debug=True, host="0.0.0.0")
