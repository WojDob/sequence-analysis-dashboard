# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from io import StringIO
from dash import Dash, html, dcc
from dash.dependencies import Input, Output, State
from Bio import SeqIO

app = Dash(__name__)

app.layout = html.Div(
    [
        dcc.Textarea(
            id="textarea-state-example",
            value=">example_1\nVLSISYSRSESSLETIGQRKPSTFSWSSRAASRSSWERGP",
            style={"width": "100%", "height": 200},
        ),
        html.Button("Submit", id="textarea-state-example-button", n_clicks=0),
        html.Div(id="textarea-state-example-output", style={"whiteSpace": "pre-line"}),
    ]
)

@app.callback(
    Output("textarea-state-example-output", "children"),
    Input("textarea-state-example-button", "n_clicks"),
    State("textarea-state-example", "value"),
)
def update_output(n_clicks, value):
    if n_clicks > 0:
        fasta_io = StringIO(value)
        records = SeqIO.parse(fasta_io, "fasta")
        teststring = ""
        for rec in records:
            teststring += rec
            fasta_io.close()
        return f"You have entered: \n{teststring}"


if __name__ == "__main__":
    app.run_server(debug=True)
