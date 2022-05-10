# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.
from io import StringIO
from dash import Dash, html, dcc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from Bio import SeqIO
import pandas as pd
import plotly.express as px


app = Dash(__name__)


@app.callback(
    [Output("textarea-state-example-output", "children"), Output("chart", "figure")],
    Input("textarea-state-example-button", "n_clicks"),
    State("textarea-state-example", "value"),
)
def update_output(n_clicks, value):
    if n_clicks > 0:
        fasta_io = StringIO(value)
        record = next(SeqIO.parse(fasta_io, "fasta"))
        fasta_io.close()
        contents = {}
        for symbol in record.seq:
            if symbol in contents:
                contents[symbol] += 1
            else:
                contents[symbol] = 1
        symbols = []
        counts = []
        for key in contents:
            symbols.append(key)
            counts.append(contents[key])

        output = f"""
            id: {record.id}
            description {record.description}
            symbols: {symbols}
            counts: {counts}
            length: {len(record.seq)}
        """
        df = pd.DataFrame(
            {
                "Symbols": symbols,
                "Counts": counts,
            }
        )
        fig = px.bar(df, x="Symbols", y="Counts", color="Counts")
        fig.update_layout(barmode="overlay")
        return f"You have entered: \n{output}", fig
    raise PreventUpdate


app.layout = html.Div(
    [
        dcc.Textarea(
            id="textarea-state-example",
            value=">example_1\nVLSISYSRSESSLETIGQRKPSTFSWSSRAASRSSWERGP",
            style={"width": "100%", "height": 200},
        ),
        html.Button("Submit", id="textarea-state-example-button", n_clicks=0),
        html.Div(id="textarea-state-example-output", style={"whiteSpace": "pre-line"}),
        html.Div(dcc.Graph(id="chart")),
    ]
)


if __name__ == "__main__":
    app.run_server(debug=True)
