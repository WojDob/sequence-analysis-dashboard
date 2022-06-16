import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from plotly.graph_objs._figure import Figure
from plotly.graph_objs._pie import Pie
from plotly.graph_objs._bar import Bar

import app


@pytest.fixture(name="sequence")
def fixture_sequence():
    seq = Seq("GATCGGATCG")
    record = SeqRecord(seq)
    record.id = "test"
    record.description = "test"
    return record


def test_should_return_correct_Seq_object():
    seq = ">test\nATGCATGC"
    result = app.string_to_seq(seq)

    assert isinstance(result, SeqRecord)
    assert result.id == "test"
    assert result.description == "test"
    assert result.seq == Seq("ATGCATGC")


def test_should_return_proper_statistics_string(sequence):

    assert (
        app.get_sequence_statistics(sequence)
        == """
id: test
description: test
length: 10
complementary sequence: CGATCCGATC
melting temperature of a sequence: 32.0
translated sequence: DRI"""
    )


def test_should_return_mixed_RNA_DNA_error_message_when_both_U_and_T_in_sequence():
    seq = Seq("GATCGGAUCG")
    record = SeqRecord(seq)
    record.id = "test"
    record.description = "test"
    assert (
        app.get_sequence_statistics(record)
        == """
id: test
description: test
length: 10
ERROR - Mixed RNA/DNA found"""
    )


@pytest.mark.parametrize(
    "sequence, expected_result",
    [("ATGC", True), ("XD", False), ("XDATGCPL", False), ("ATGCNU", True)],
)
def test_should_check_if_sequence_contains_only_nucleotides(sequence, expected_result):
    assert app.is_nucleotide(sequence) == expected_result


def test_should_count_symbol_occurences(sequence):
    symbols, counts = app.count_symbols(sequence.seq)

    assert len(symbols) == len(counts) == 4
    assert symbols == ["G", "A", "T", "C"]
    assert counts == [4, 2, 2, 2]
    assert sum(counts) == 10


def test_should_return_bar_figure(sequence):
    fig = app.get_sequence_bar_figure(sequence)
    assert isinstance(fig, Figure) is True
    assert isinstance(fig["data"][0], Bar) is True


def test_should_return_pie_figure(sequence):
    fig = app.get_sequence_pie_figure(sequence)
    assert isinstance(fig, Figure) is True
    assert isinstance(fig["data"][0], Pie) is True
