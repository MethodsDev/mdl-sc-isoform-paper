import igraph as ig
import numpy as np
from networkx.drawing.nx_agraph import graphviz_layout


def layout_graph(graph: ig.Graph):
    """
    Given a graph, computes a two-dimensional layout for visualization
    """
    pos = graphviz_layout(graph.to_networkx(), prog="sfdp")
    return np.array([pos[i] for i in range(len(pos))])
