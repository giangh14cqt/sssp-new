//
// Created by Truong Giang Do on 28/11/2023.
//
#include "LDD.h"

vector<vector<int>> preLDD(Graph &g, int d, double CALCULATE_SCC_PROB) {
    double r = ((double) Random::Get().GenInt() / (RAND_MAX));
    if (r < CALCULATE_SCC_PROB)
        return LDD(g, d);

    vector<vector<int>> SCCs = g.SCC();
    vector<vector<int>> E_sep;

    if (SCCs.size() == 1)
        return LDD(g, d);

    for (vector<int> SCC: SCCs)
        if (SCC.size() > 1) {
            Graph SCCSubgraph(g.v_max, false);
            SCCSubgraph.addVertices(SCC);

            set<int> SCCVerts(SCC.begin(), SCC.end());

            for (int v: SCC) {
                vector<int> outVertices;
                vector<int> weights;

                for (int i = 0; i < g.adjacencyList[v].size(); i++) {
                    if (SCCVerts.find(g.adjacencyList[v][i]) != SCCVerts.end()) {
                        outVertices.push_back(g.adjacencyList[v][i]);
                        weights.push_back(g.weights[v][i]);
                    }
                }

                SCCSubgraph.addEdges(v, outVertices, weights);
            }

            int src = SCC[rand() % SCC.size()];
            if (hasLargeDiameter(SCCSubgraph, src, d)) {
                vector<vector<int>> SCCSubgraphLDD = LDD(SCCSubgraph, d);
                for (const auto &vv: SCCSubgraphLDD) {
                    E_sep.push_back(vv);
                }
            }
        }
    return E_sep;
}

bool hasLargeDiameter(Graph &g, int s, int diameter) {
    vector<bool> settled(g.v_max, false);
    int numSettled = 0;
    priority_queue<Node> pq;
    vector<int> dist(g.v_max, INT_MAX);
    pq.emplace(s, 0);
    dist[s] = 0;

    while (numSettled != g.n) {
        if (pq.empty())
            return false;

        int u = pq.top().node;
        pq.pop();

        if (settled[u])
            continue;

        if (dist[u] > diameter)
            return true;

        settled[u] = true;
        numSettled++;

        for (int i = 0; i < g.adjacencyList[u].size(); i++) {
            int v = g.adjacencyList[u][i];
            if (!settled[v] && dist[u] + g.weights[u][i] < dist[v]) {
                dist[v] = dist[u] + g.weights[u][i];
                pq.emplace(v, dist[v]);
            }
        }
    }
    return false;
}

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
vector<vector<int>> LDD(Graph &g, int d) {
    int LDD_BASE_CASE = 10;
    if (g.n <= max(1, LDD_BASE_CASE))
        return {};

    Graph g_rev = createGRev(g);

    // pick a good, well-connected source vertex
    int s = g.vertices[0];
    bool foundGoodS = false;
    for (int v: g.vertices) {
        if (g.adjacencyList[v].size() >= 1 || !g_rev.adjacencyList[v].empty()) {
            s = v;
            foundGoodS = true;
            break;
        }
    }

    if (!foundGoodS)
        return {};

    vector<int> condAndi_max = CoreOrLayerRange(g, g_rev, s, d);

    if (condAndi_max[0] == 1)
        return RandomTrim(g, g_rev, s, d);

    int r = (int) ceil(d / (3.0 * log(g.n)));
    int i_min = condAndi_max[1] - r;

    random_device rd;
    mt19937 gen(rd());
    geometric_distribution<int> distr(calculateGeoProb(g.n, r));

    int i_tilda = distr(gen);
    int i_rnd = i_min + min(i_tilda, r);

    if (condAndi_max[0] == 2) {
        vector<int> ball = volume(g, s, i_rnd);
        Graph subGraph = getSubGraph(g, ball, false);
        Graph minusSubGraph = getSubGraph(g, ball, true);

        vector<vector<int>> layer_g = layer(g, ball);
        vector<vector<int>> preLDD_subGraph = preLDD(subGraph, d, subGraph.n/1000.0);
        vector<vector<int>> preLDD_minusSubGraph = preLDD(minusSubGraph, d, minusSubGraph.n/1000.0);
        return edgeUnion(layer_g, preLDD_subGraph, preLDD_minusSubGraph);
    }

    if (condAndi_max[0] == 3) {
        vector<int> ball = volume(g_rev, s, i_rnd);
        Graph subGraph = getSubGraph(g_rev, ball, false);
        Graph minusSubGraph = getSubGraph(g_rev, ball, true);

        vector<vector<int>> layer_g_rev = layer(g_rev, ball);
        vector<vector<int>> preLDD_subGraph = preLDD(subGraph, d, subGraph.n/1000.0);
        vector<vector<int>> preLDD_minusSubGraph = preLDD(minusSubGraph, d, minusSubGraph.n/1000.0);
        layer_g_rev = revEdges(layer_g_rev);
        preLDD_subGraph = revEdges(preLDD_subGraph);
        preLDD_minusSubGraph = revEdges(preLDD_minusSubGraph);
        return edgeUnion(layer_g_rev, preLDD_subGraph, preLDD_minusSubGraph);
    }

    throw "Error: condAndi_max[0] is not 1, 2, or 3";
}

vector<vector<int>> revEdges(vector<vector<int>> &edges) {
    vector<vector<int>> revEdgeSet(edges.size());
    for (int i = 0; i < edges.size(); i++) {
        revEdgeSet[i] = vector<int>{edges[i][1], edges[i][0]};
    }
    return revEdgeSet;
}

double calculateGeoProb(int n, int r) {
    double logn = log(n);
    double prob = logn * logn / r;
    return min(1.0, prob);
}

vector<vector<int>> RandomTrim(Graph &g, Graph &g_rev, int s, int d) {
    vector<vector<int>> e_sep;

    if (INT_MAX / 4 <= d)
        return e_sep;

    vector<int> dist = Dijkstra(g, s);
    vector<int> dist_rev = Dijkstra(g_rev, s);

    int numDistOver4D = 0;
    int over4DVert = -1;

    vector<int> v_far;
    for (int v: g.vertices)
        if (max(dist[v], dist_rev[v]) > 2 * d) {
            v_far.push_back(v);
            if (max(dist[v], dist_rev[v]) > 4 * d) {
                numDistOver4D++;
                over4DVert = v;
            }
        }

    // base case
    if (numDistOver4D <= 1) {
        if (numDistOver4D == 1) {
            for (int v: g.adjacencyList[over4DVert]) {
                e_sep.push_back(vector<int>{over4DVert, v});
            }
        }
        return e_sep;
    }

    vector<int> m;
    int i_max = d;
    int r = (int) ceil(d / (3.0 * log(g.n)));
    int i_min = i_max - r;

    int v = diffVertex(v_far, m, g.v_max);

    random_device rd;
    mt19937 gen(rd());
    geometric_distribution<int> distr(calculateGeoProb(g.n, r));

    while (v != -1) {
        int i_rnd = i_min + min(distr(gen), r);

        if (dist_rev[v] > 2 * d) {
            Graph gVminusM = getSubGraph(g, m, true);
            vector<int> ball = volume(gVminusM, v, i_rnd);
            Graph GVMinusMSubGraph = getSubGraph(gVminusM, ball, false);
            vector<vector<int>> layer_gVminusM = layer(gVminusM, ball);
            vector<vector<int>> preLDD_GVMinusMSubGraph = preLDD(GVMinusMSubGraph, d, GVMinusMSubGraph.n/1000.0);
            e_sep = edgeUnion(e_sep, layer_gVminusM, preLDD_GVMinusMSubGraph);
            m = vertexUnion(m, ball);
        } else if (dist[v] > 2 * d) {
            Graph gVminusM_rev = getSubGraph(g_rev, m, true);
            vector<int> ball_rev = volume(gVminusM_rev, v, i_rnd);
            Graph GVMinusMSubGraph_rev = getSubGraph(gVminusM_rev, ball_rev, false);
            vector<vector<int>> layer_gVminusM_rev = layer(gVminusM_rev, ball_rev);
            vector<vector<int>> preLDD_GVMinusMSubGraph_rev = preLDD(GVMinusMSubGraph_rev, d, GVMinusMSubGraph_rev.n/1000.0);
            layer_gVminusM_rev = revEdges(layer_gVminusM_rev);
            preLDD_GVMinusMSubGraph_rev = revEdges(preLDD_GVMinusMSubGraph_rev);
            e_sep = edgeUnion(e_sep, layer_gVminusM_rev, preLDD_GVMinusMSubGraph_rev);
            m = vertexUnion(m, ball_rev);
        } else
            throw "Error: dist[v] and dist_rev[v] are both <= 2 * d";

        v = diffVertex(v_far, m, g.v_max);
    }
    return e_sep;
}

// returns the subgraph of g containing only the vertices in ball
// if setMinus is true, the function returns the subgraph of g containing only
// the vertices outside of the ball
Graph getSubGraph(Graph &g, vector<int> &ball, bool setMinus) {
    vector<bool> contains(g.v_max, false);
    for (int i: ball)
        contains[i] = true;

    vector<int> vert;
    for (int v: g.vertices)
        if (!setMinus && contains[v])
            vert.push_back(v);
        else if (setMinus && !contains[v])
            vert.push_back(v);

    Graph subGraph(g.v_max, false);
    subGraph.addVertices(vert);

    for (int u: vert) {
        vector<int> edges;
        vector<int> weights;

        for (int i = 0; i < g.adjacencyList[u].size(); i++) {
            int v = g.adjacencyList[u][i];

            if (subGraph.containsVertex[v]) {
                edges.push_back(v);
                weights.push_back(g.weights[u][i]);
            }
        }

        subGraph.addEdges(u, edges, weights);
    }
    return subGraph;
}

// returns the union of two vertex sets
vector<int> vertexUnion(vector<int> &set1, vector<int> &set2) {
    set<int> set;
    addVerticesToSet(set, set1);
    addVerticesToSet(set, set2);

    vector<int> output;
    output.reserve(set.size());
    for (int i: set)
        output.push_back(i);

    return output;
}

void addVerticesToSet(set<int> &set, vector<int> &vertices) {
    for (int i: vertices)
        set.insert(i);
}

vector<vector<int>> edgeUnion(vector<vector<int>> &set1,
                              vector<vector<int>> &set2,
                              vector<vector<int>> &set3) {
    set<vector<int>> set(set1.begin(), set1.end());
    set.insert(set2.begin(), set2.end());
    set.insert(set3.begin(), set3.end());
    return {set.begin(), set.end()};
}

int diffVertex(vector<int> &set1, vector<int> &set2, int v_max) {
    vector<bool> contains(v_max, false);
    for (int i : set2)
        contains[i] = true;

    for (int i : set1)
        if (!contains[i])
            return i;

    return -1;
}

// OUTPUT: a pair (Condition,i) where Condition ∈ {1,2,3} and i ≤ D is a
// non-negative integer such that
// – if Condition = 1 then n_G(s,i) > 2n and n_G_rev(s,i) > 2n,
// - if Condition = 2 then n_G(s, i) ≤ 2n and Vol_G(s, i) and Vol_G(s, i − ⌈D/(3lgn)⌉) are the same canonical range,
// – if Condition = 3 then n_G_rev(s, i) ≤ 2n and Vol_G_rev(s, i) and Vol_G_rev(s, i − ⌈D/(3lgn)⌉) are in the same canonical range.
// If Condition ∈ {2, 3} then i ≥ D/(3lgn).
// Runs LayerRange on G and G_rev in parallel.
vector<int> CoreOrLayerRange(Graph &g, Graph &g_rev, int s, int d) {
    vector<vector<int>> farthestDistancesSeen;
    vector<vector<int>> farthestDistancesSeen_rev;
    double constant = d / (3.0 * log(g.n));
    vector<bool> settled(g.v_max, false);
    vector<bool> settled_rev(g.v_max, false);
    int numSettled = 0;
    int numSettled_rev = 0;
    priority_queue<Node> pq;
    priority_queue<Node> pq_rev;
    vector<int> dist(g.v_max, INT_MAX);
    vector<int> dist_rev(g.v_max, INT_MAX);
    pq.push(Node(s, 0));
    pq_rev.push(Node(s, 0));
    dist[s] = 0;
    dist_rev[s] = 0;
    bool finished = false;
    bool finished_rev = false;
    int j = -1;
    int j_rev = -1;

    while (true) {
        if (numSettled == g.n)
            finished = true;

        if (numSettled_rev == g_rev.n)
            finished_rev = true;

        if (finished && finished_rev) {
            // case 1
            if (j == -1 || j_rev == -1)
                return vector<int>{2, (int) ceil(constant)};

            return vector<int>{1, max(j, j_rev)};
        }

        if (!finished) {
            vector<int> result = oneIterationLayerRange(g, pq, settled, numSettled, farthestDistancesSeen, constant,
                                                        dist, d);
            for (vector<int> i: farthestDistancesSeen)
                for (int i: result)
                    if (!result.empty()) {
                        if (result[0] == 1) {
                            j = result[1];
                            finished = true;
                        } else if (result[0] == 2) {
                            // case 2
                            return vector<int>{2, result[1]};
                        }
                    }
            numSettled++;
        }

        if (!finished_rev) {
            vector<int> result_rev = oneIterationLayerRange(g_rev, pq_rev, settled_rev, numSettled_rev,
                                                            farthestDistancesSeen_rev, constant, dist_rev, d);
            for (vector<int> i: farthestDistancesSeen_rev)
                for (int i: result_rev)
                    if (!result_rev.empty()) {
                        if (result_rev[0] == 1) {
                            j_rev = result_rev[1];
                            finished_rev = true;
                        } else if (result_rev[0] == 2) {
                            // case 3
                            return vector<int>{3, result_rev[1]};
                        }
                    }
            numSettled_rev++;
        }
    }
}

// Checked
Graph createGRev(Graph &g) {
    Graph g_rev(g.v_max, false);
    g_rev.addVertices(g.vertices);

    vector<vector<int>> edges(g.v_max);
    vector<vector<int>> weights(g.v_max);

    for (int v: g.vertices) {
        for (int i = 0; i < g.adjacencyList[v].size(); i++) {
            edges[g.adjacencyList[v][i]].push_back(v);
            weights[g.adjacencyList[v][i]].push_back(g.weights[v][i]);
        }
    }

    for (int i = 0; i < g.v_max; i++) {
        g_rev.addEdges(i, edges[i], weights[i]);
    }
    return g_rev;
}

vector<int> oneIterationLayerRange(Graph &g,
                                   priority_queue<Node> &pq,
                                   vector<bool> &settled,
                                   int numSettled,
                                   vector<vector<int>> &farthestDistancesSeen,
                                   double constant,
                                   vector<int> &dist,
                                   int d) {
    if (pq.empty()) {
        /*
         * Nothing left to search.
         * g is disconnected, since n_G(s,i) <= 2n/3.
         * Will never be the case that i ≥ D/(3lgn).
         * Return i_big = i + ceil(D/(3lgn)).
         * Guaranteed that i_big will satisfy i_big >= D/(3lgn) and the canonical
         * ranges.
         */
        int farthestDistanceSeen = farthestDistancesSeen[farthestDistancesSeen.size() - 1][0];
        int i_big = min(d, farthestDistanceSeen + (int) ceil(constant));
        return vector<int>{2, i_big};
    }
    int u = pq.top().node;
    pq.pop();

    if (settled[u])
        return {};

    settled[u] = true;

    if (farthestDistancesSeen.empty() || dist[u] > farthestDistancesSeen[farthestDistancesSeen.size() - 1][0]) {
        farthestDistancesSeen.push_back(vector<int>{dist[u], numSettled + 1});
    }

    int farthesDistanceSeen = farthestDistancesSeen[farthestDistancesSeen.size() - 1][0];

    // case 1
    if (numSettled + 1 > 2.0 * g.n / 3.0)
        return vector<int>{1, farthesDistanceSeen};

    // case 2
    if (farthesDistanceSeen >= constant && sameCanonicalRange(farthestDistancesSeen, constant))
        return vector<int>{2, farthesDistanceSeen};

    updateNeighbors(g, u, settled, pq, dist, d);

    return {};
}

// checked
// Checks whether Vol_G(s, i - ceil[D/(3logn)]) and Vol_G(s, i) are in the same canonical range.
// Two numbers are in the same canonical range if they lie in the same half-open interval
// [2^j, 2^{j+1}), where j is a non-negative integer.
bool sameCanonicalRange(vector<vector<int>> &farthestDistancesSeen, double constant) {
    int i = farthestDistancesSeen[farthestDistancesSeen.size() - 1][0];
    int vol1 = farthestDistancesSeen[farthestDistancesSeen.size() - 1][1];

    for (int j = farthestDistancesSeen.size() - 2; j >= 0; j--)
        if (farthestDistancesSeen[j][0] <= i - ceil(constant)) {
            int vol2 = farthestDistancesSeen[j][1];
            if (floor(log(vol1) / log(2)) == floor(log(vol2) / log(2)))
                return true;
            break;
        }
    return false;
}

// checked
// returns the edges (u, v) where u is in ball and v is not in ball
vector<vector<int>> layer(Graph &g, vector<int> &ball) {
    vector<bool> contains(g.v_max, false);
    for (int i: ball)
        contains[i] = true;
    vector<vector<int>> edges;
    for (int u: ball) {
        for (int v: g.adjacencyList[u]) {
            if (!contains[v])
                edges.push_back(vector<int>{u, v});
        }
    }
    return edges;
}

// checked
// returns all the vertices in g within a distance of r from source vertex s using Dijkstra's
vector<int> volume(Graph &g, int s, int r) {
    vector<int> output;
    vector<bool> settled(g.v_max, false);
    int numSettled = 0;
    priority_queue<Node> pq;
    vector<int> dist(g.v_max, INT_MAX);
    pq.push(Node(s, 0));
    dist[s] = 0;

    while (numSettled != g.n) {
        if (pq.empty())
            return output;

        int u = pq.top().node;
        pq.pop();

        if (settled[u] || dist[u] > r)
            continue;

        output.push_back(u);
        settled[u] = true;
        numSettled++;

        updateNeighbors(g, u, settled, pq, dist, r);
    }
    return output;
}

// checked
vector<int> Dijkstra(Graph &g, int s) {
    vector<bool> settled(g.v_max, false);
    int numSettled = 0;
    priority_queue<Node> pq;
    vector<int> dist(g.v_max, INT_MAX);
    pq.push(Node(s, 0));
    dist[s] = 0;

    while (numSettled != g.n) {
        if (pq.empty())
            return dist;

        int u = pq.top().node;
        pq.pop();

        if (settled[u])
            continue;

        settled[u] = true;
        numSettled++;

        updateNeighbors(g, u, settled, pq, dist, INT_MAX);
    }
    return dist;
}

// checked
void updateNeighbors(Graph &g, int u, vector<bool> &settled, priority_queue<Node> &pq, vector<int> &dist, int d) {
    for (int i = 0; i < g.adjacencyList[u].size(); i++) {
        int v = g.adjacencyList[u][i];
        if (!settled[v]) {
            int newDistance = dist[u] + g.weights[u][i];

            dist[v] = min(dist[v], newDistance);

            // only want to process nodes within a distance of d from the source
            if (dist[v] <= d)
                pq.push(Node(v, dist[v]));
        }
    }
}