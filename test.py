import networkx as nx
import networkx.algorithms.isomorphism as iso


def isSame(x, y):
    print type(x)
    print x, ",", y, x == y
    return x == y


def main():
    DG1 = nx.DiGraph()
    #DG1.add_nodes_from(["1", "2", "3", "4"])
    for i in range(1, 5):
        DG1.add_node(i, label=i)
    DG1.add_edge(1, 2)
    DG1.add_edge(1, 3)

    DG2 = nx.DiGraph()
    for i in range(1, 5):
        DG2.add_node(i, label=i)
    DG2.add_edge(1, 2)
    DG2.add_edge(1, 3)

    print DG2.nodes

    nm = iso.categorical_node_match("label", "label")

    print nx.is_isomorphic(DG1, DG2, nm)


if __name__=="__main__":
    main()