import numpy as np
import heapq as hq
from typing import Union

class Graph:
    def __init__(self, adjacency_mat: Union[np.ndarray, str]):
        """ Unlike project 2, this Graph class takes an adjacency matrix as input. 
        `adjacency_mat` can either be a 2D numpy array of floats or 
        the path to a CSV file containing a 2D numpy array of floats.

        In this project, we will assume `adjacency_mat` corresponds to the 
        adjacency matrix of an undirected graph
        """
        if type(adjacency_mat) == str:
            self.adj_mat = self._load_adjacency_matrix_from_csv(adjacency_mat)
        elif type(adjacency_mat) == np.ndarray:
            self.adj_mat = adjacency_mat
        else: 
            raise TypeError('Input must be a valid path or an adjacency matrix')
        self.mst = None

    def _load_adjacency_matrix_from_csv(self, path: str) -> np.ndarray:
        with open(path) as f:
            return np.loadtxt(f, delimiter=',')

    def construct_mst(self):
        """ Given `self.adj_mat`, the adjacency matrix of a connected undirected graph, 
        implement Prim's algorithm to construct an adjacency matrix encoding the minimum 
        spanning tree of `self.adj_mat`. 
            
        `self.adj_mat` is a 2D numpy array of floats. 
        Note that because we assume our input graph is undirected, `self.adj_mat` is symmetric. 
        Row i and column j represents the edge weight between vertex i and vertex j. 
        An edge weight of zero indicates that no edge exists. 
        
        TODO: 
        This function does not return anything. Instead, store the adjacency matrix 
        representation of the minimum spanning tree of `self.adj_mat` in `self.mst`.
        We highly encourage the use of priority queues in your implementation. See the heapq
        module, particularly the `heapify`, `heappop`, and `heappush` functions.
        """
        n = len(self.adj_mat)
        mst = np.zeros((n,n))
        mstWeight = 0
        inMST = [False] * n
        edgeHeap = []
        # define starting node (default 0 but could be anything)
        u = 0
        # initialize edge priority queue from starting node
        inMST[u] = True
        for v in range(0,n):
        	if not(self.adj_mat[u][v] == 0):
        		hq.heappush(edgeHeap, (self.adj_mat[u][v], u, v))
        # iterate through edges in edgeHeap until all nodes are in mst
        while sum(inMST) < n:
        	# check that edgeHeap isn't empty
        	if edgeHeap == []:
        		raise ValueError("Tree spanning all nodes doesn't exist")
        	# get minimum edge
        	nextEdge = hq.heappop(edgeHeap)
        	u = nextEdge[1]
        	v = nextEdge[2]
        	if not(inMST[v]):
        		# add node to mst
        		inMST[v] = True
        		# add edge to mst
        		mst[u][v] = nextEdge[0]
        		# since the graph is undirected, add symmetrical edge to mst
        		mst[v][u] = nextEdge[0]
        		# add edge weight to total mst weight
        		mstWeight += nextEdge[0]
        		# add edges from v to edgeHeap
        		for w in range(0,n):
        			if not(self.adj_mat[v][w] == 0):
        				hq.heappush(edgeHeap, (self.adj_mat[v][w], v, w))
        self.mst = mst
        self.mst_weight = mstWeight

'''
test = Graph('./data/small.csv')
test.construct_mst()
print("Final MST:")
print(test.mst)
print("Final MST Weight:")
print(test.mst_weight)

test2 = Graph('./data/no_span.csv')
test2.construct_mst()

test3 = Graph('./data/no_span2.csv')
test3.construct_mst()
'''