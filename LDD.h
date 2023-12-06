//
// Created by Truong Giang Do on 28/11/2023.
//

#ifndef SSSP_NEW_LDD_H
#define SSSP_NEW_LDD_H

#include "Randd.h"

vector<vector<int>> preLDD(Graph &g, int d, double CALCULATE_SCC_PROB);

bool hasLargeDiameter(Graph &g, int s, int diameter);

/*
OUTPUT: A set of edges e_sep with the following guarantees:
– each SCC of G\e_sep has weak diameter at most D; that is, if u,v are in the
same SCC,
then dist_G(u, v) ≤ D and dist_G(v, u) ≤ D.
– For every e ∈ E, Pr[e ∈ e_sep] = O(w(e)·(logn)^2/D +n^−10). These
probabilities are not
guaranteed to be independent.
Each int[] in the output ArrayList has size two and represents an edge
(int[0], int[1])
*/
vector<vector<int>> LDD(Graph &g, int d);

vector<vector<int>> revEdges(vector<vector<int>> &edges);

double calculateGeoProb(int n, int r);

vector<vector<int>> RandomTrim(Graph &g, Graph &g_rev, int s, int d);

// returns the subgraph of g containing only the vertices in ball
// if setMinus is true, the function returns the subgraph of g containing only
// the vertices outside of the ball
Graph getSubGraph(Graph &g, vector<int> &ball, bool setMinus);

// returns the union of two vertex sets
vector<int> vertexUnion(vector<int> &set1, vector<int> &set2);

void addVerticesToSet(set<int> &set, vector<int> &vertices);

vector<vector<int>> edgeUnion(vector<vector<int>> &set1,
                              vector<vector<int>> &set2,
                              vector<vector<int>> &set3);

void addEdgesToSet(set<vector<int>> &set, vector<vector<int>> &edges);

int diffVertex(vector<int> &set1, vector<int> &set2, int v_max);

// OUTPUT: a pair (Condition,i) where Condition ∈ {1,2,3} and i ≤ D is a
// non-negative integer such that
// – if Condition = 1 then n_G(s,i) > 2n and n_G_rev(s,i) > 2n,
// - if Condition = 2 then n_G(s, i) ≤ 2n and Vol_G(s, i) and Vol_G(s, i −
// ⌈D/(3lgn)⌉) are the same canonical range,
// – if Condition = 3 then n_G_rev(s, i) ≤ 2n and Vol_G_rev(s, i) and
// Vol_G_rev(s, i − ⌈D/(3lgn)⌉) are in the same canonical range.
// If Condition ∈ {2, 3} then i ≥ D/(3lgn).
// Runs LayerRange on G and G_rev in parallel.
vector<int> CoreOrLayerRange(Graph &g, Graph &g_rev, int s, int d);

Graph createGRev(Graph &g);

vector<int> oneIterationLayerRange(Graph &g, priority_queue<Node> &pq, vector<bool> &settled,
                                   int numSettled, vector<vector<int>> &farthestDistancesSeen, double constant, vector<int> &dist, int d);

bool sameCanonicalRange(vector<vector<int>> &farthestDistancesSeen, double constant);

vector<vector<int>> layer(Graph &g, vector<int> &ball);

// returns all the vertices in g within a distance of r from source vertex s
// using Dijkstra's
vector<int> volume(Graph &g, int s, int r);

vector<int> Dijkstra(Graph &g, int s);

// void init(Graph g, priority_queue<Node> pq, vector<int> dist, int s);
void updateNeighbors(Graph &g, int u, vector<bool> &settled, priority_queue<Node> &pq, vector<int> &dist, int d);


#endif //SSSP_NEW_LDD_H
