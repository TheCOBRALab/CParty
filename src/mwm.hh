#ifndef MWM_HEADER
#define MWM_HEADER

#include "mea.hh"
#include "base_types.hh"

#include <algorithm>
/**Weighted maximum matching in general graphs.
 * 
The algorithm is taken from "Efficient Algorithms for Finding Maximum
Matching in Graphs" by Zvi Galil, ACM Computing Surveys, 1986.
It is based on the "blossom" method for finding augmenting paths and
the "primal-dual" method for finding a matching of maximum weight, both
due to Jack Edmonds.
Some ideas came from "Implementation of algorithms for maximum matching
on non-bipartite graphs" by H.J. Gabow, Standford Ph.D. thesis, 1973.

A C program for maximum weight matching by Ed Rothberg was used extensively
to validate this new code.
**/
class maxWeightMatching{
    /**Compute a maximum-weighted matching in the general undirected
    weighted graph given by "edges".  If "maxcardinality" is true,
    only maximum-cardinality matchings are considered as solutions.

    Edges is a sequence of tuples (i, j, wt) describing an undirected
    edge between vertex i and vertex j with weight wt.  There is at most
    one edge between any two vertices; no vertex has an edge to itself.
    Vertices are identified by consecutive, non-negative integers.

    Return a list "mate", such that mate[i] == j if vertex i is
    matched to vertex j, and mate[i] == -1 if vertex i is not matched.

    This function takes time O(n ** 3).
    **/
    cand_pos_t n;
    cand_pos_t nedge;
    std::vector<elem_prob_s> edges;
    pf_t maxweight = 0;
    // mate[v] is the remote endpoint of its matched edge, or -1 if it is single i.e the pairing; intitially all unpaired
    std::vector<cand_pos_t> mate;
    // If b is a top-level blossom, label[b] is 0 if b is unlabeled (free); 1 if b is an S-vertex/blossom; 2 if b is a T-vertex/blossom.
    // The label of a vertex is found by looking at the label of its top-level containing blossom.
    // If v is a vertex inside a T-blossom, label[v] is 2 iff v is reachable from an S-vertex outside the blossom.
    // Labels are assigned during a stage and reset after each augmentation.
    std::vector<cand_pos_t> label;
    // If p is an edge endpoint, endpoint[p] is the vertex to which endpoint p is attached.
    // Not modified by the algorithm.
    std::vector<cand_pos_t> endpoint;
    // If v is a vertex, neighbend[v] is the list of remote endpoints of the edges attached to v.
    // Not modified by the algorithm.
    std::vector< std::vector<cand_pos_t> > neighbend;
    // If b is a labeled top-level blossom, labelend[b] is the remote endpoint of the edge through which b obtained its label, or -1 if b's base vertex is single.
    // If v is a vertex inside a T-blossom and label[v] == 2, labelend[v] is the remote endpoint of the edge through which v is reachable from outside the blossom.
    std::vector<cand_pos_t> labelend;
    // If v is a vertex, inblossom[v] is the top-level blossom to which v belongs.
    // If v is a top-level vertex, v is itself a blossom (a trivial blossom) and inblossom[v] == v.
    // Initially all vertices are top-level trivial blossoms.
    std::vector<cand_pos_t> inblossom;
    // If b is a sub-blossom, blossomparent[b] is its immediate parent (sub-)blossom.
    // If b is a top-level blossom, blossomparent[b] is -1.
    std::vector<cand_pos_t> blossomparent;
    // If b is a non-trivial (sub-)blossom, blossomchilds[b] is an ordered list of its sub-blossoms, starting with
    // the base and going round the blossom.
    std::vector< std::vector<cand_pos_t> > blossomchilds;
    // If b is a (sub-)blossom, blossombase[b] is its base VERTEX (i.e. recursive sub-blossom).
    std::vector<cand_pos_t> blossombase;
    // If b is a non-trivial (sub-)blossom, blossomendps[b] is a list of endpoints on its connecting edges,
    // such that blossomendps[b][i] is the local endpoint of blossomchilds[b][i] on the edge that connects it to blossomchilds[b][wrap(i+1)].
    std::vector< std::vector<cand_pos_t> > blossomendps;
    // If v is a free vertex (or an unreached vertex inside a T-blossom), bestedge[v] is the edge to an S-vertex with least slack, or -1 if there is no such edge.
    // If b is a (possibly trivial) top-level S-blossom, bestedge[b] is the least-slack edge to a different S-blossom, or -1 if there is no such edge.
    // This is used for efficient computation of delta2 and delta3.
    std::vector<cand_pos_t> bestedge;
    // If b is a non-trivial top-level S-blossom, blossombestedges[b] is a list of least-slack edges to neighbouring
    // S-blossoms, or None if no such list has been computed yet. This is used for efficient computation of delta3.
    std::vector<std::vector<cand_pos_t> > blossombestedges;
    // List of currently unused blossom numbers.
    std::vector<cand_pos_t> unusedblossoms;
    // If v is a vertex, dualvar[v] = 2 * u(v) where u(v) is the v's variable in the dual
    // optimization problem (multiplication by two ensures integer values throughout the algorithm if all edge weights are integers).
    // If b is a non-trivial blossom, dualvar[b] = z(b) where z(b) is b's variable in the dual optimization problem.
    std::vector<pf_t> dualvar;
    // If allowedge[k] is true, edge k has zero slack in the optimization
    // problem; if allowedge[k] is false, the edge's slack may or may not be zero.
    std::vector<char> allowedge;
    // Queue of newly discovered S-vertices.
    std::vector<cand_pos_t> queue;
    public:
        maxWeightMatching(std::vector<elem_prob_s> &edges,cand_pos_t n);

    private:
        // Return 2 * slack of edge k (does not work inside blossoms).
        pf_t slack(cand_pos_t k);
        // Generate the leaf vertices of a blossom.
        void blossomLeaves(cand_pos_t b,std::vector<cand_pos_t> &leaves);

        // Assign label t to the top-level blossom containing vertex w
        // and record the fact that w was reached through the edge with
        // remote endpoint p.
        void assignLabel(cand_pos_t w, cand_pos_t t, cand_pos_t p);

        // Trace back from vertices v and w to discover either a new blossom
        // or an augmenting path. Return the base vertex of the new blossom or -1.
        cand_pos_t scanBlossom(cand_pos_t v, cand_pos_t w);

        // Construct a new blossom with given base, containing edge k which
        // connects a pair of S vertices. Label the new blossom as S; set its dual
        // variable to zero; relabel its T-vertices to S and add them to the queue.
        void addBlossom(cand_pos_t base, cand_pos_t k);
        // Expand the given top-level blossom.
        void expandBlossom(cand_pos_t b, bool endstage);
        // Swap matched/unmatched edges over an alternating path through blossom b
        // between vertex v and the base vertex. Keep blossom bookkeeping consistent.
        void augmentBlossom(cand_pos_t b, cand_pos_t v);
        // Swap matched/unmatched edges over an alternating path between two
        // single vertices. The augmenting path runs through edge k, which
        // connects a pair of S vertices.
        void augmentMatching(cand_pos_t k);
    public:
        // Get the structure and MEA
        pf_t maximumWeightedMatching(std::vector<pf_t> &pu,std::string &structure, sparse_tree &tree);
};


#endif