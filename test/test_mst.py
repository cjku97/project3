# write tests for mst
import pytest
import numpy as np
from mst import Graph
from sklearn.metrics import pairwise_distances


def check_mst(adj_mat: np.ndarray, 
              mst: np.ndarray, 
              expected_weight: int, 
              allowed_error: float = 0.0001):
    """ Helper function to check the correctness of the adjacency matrix encoding an MST.
        Note that because the MST of a graph is not guaranteed to be unique, we cannot 
        simply check for equality against a known MST of a graph. 

        Arguments:
            adj_mat: Adjacency matrix of full graph
            mst: Adjacency matrix of proposed minimum spanning tree
            expected_weight: weight of the minimum spanning tree of the full graph
            allowed_error: Allowed difference between proposed MST weight and `expected_weight`

        TODO: 
            Add additional assertions to ensure the correctness of your MST implementation
        For example, how many edges should a minimum spanning tree have? Are minimum spanning trees
        always connected? What else can you think of?
    """
    def approx_equal(a, b):
        return abs(a - b) < allowed_error
    
    def check_symmetric(a):
    	return np.allclose(a, a.T)

    total = 0
    numEdges = 0
    for i in range(mst.shape[0]):
        for j in range(i+1):
        	if not(mst[i,j] == 0):
        		total += mst[i, j]
        		numEdges += 1
    
    all_total = 0
    all_numEdges = 0
    for i in range(adj_mat.shape[0]):
        for j in range(i+1):
        	if not(adj_mat[i,j] == 0):
        		all_total += adj_mat[i, j]
        		all_numEdges += 1
    
    assert approx_equal(total, expected_weight), 'Proposed MST has incorrect expected weight'
    assert numEdges == len(adj_mat) - 1, 'Proposed MST has unexpected number of edges'
    assert total <= all_total, 'Proposed MST has larger weight than original graph'
    assert numEdges <= all_numEdges, 'Proposed MST has more edges than original graph'
    assert check_symmetric(mst), 'Proposed MST is not symmetric (unweighted)'
    
    


def test_mst_small():
    """ Unit test for the construction of a minimum spanning tree on a small graph """
    file_path = './data/small.csv'
    g = Graph(file_path)
    g.construct_mst()
    check_mst(g.adj_mat, g.mst, 8)


def test_mst_single_cell_data():
    """ Unit test for the construction of a minimum spanning tree using 
    single cell data, taken from the Slingshot R package 
    (https://bioconductor.org/packages/release/bioc/html/slingshot.html)
    """
    file_path = './data/slingshot_example.txt'
    # load coordinates of single cells in low-dimensional subspace
    coords = np.loadtxt(file_path)
    # compute pairwise distances for all 140 cells to form an undirected weighted graph
    dist_mat = pairwise_distances(coords)
    g = Graph(dist_mat)
    g.construct_mst()
    check_mst(g.adj_mat, g.mst, 57.263561605571695)


def test_mst_student():
    """ TODO: Write at least one unit test for MST construction """
    g = Graph('./data/other_small.csv')
    g.construct_mst()
    assert g.mst_weight == 17
    check_mst(g.adj_mat, g.mst, 17)

def test_no_span_exception():
	g = Graph('./data/no_span.csv')
	with pytest.raises(ValueError, match = "Tree spanning all nodes doesn't exist"):
		g.construct_mst()

def test_no_span_exception2():
	g = Graph('./data/no_span2.csv')
	with pytest.raises(ValueError, match = "Tree spanning all nodes doesn't exist"):
		g.construct_mst()