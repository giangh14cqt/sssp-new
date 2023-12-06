//
// Created by Truong Giang Do on 01/12/2023.
//
#include "LasVegas.h"

bool MAKE_CONNECTED = false;
int SRC = 1;
bool CHECKS = false;
bool WITH_LDD = false;

Graph readInputTemp(ifstream &inputFile) {
    int g_size;
    int maxWeight;
    inputFile >> g_size >> maxWeight;
    ++g_size;
    // vector<int> phi = createRandomPriceFunction(g_size, maxWeight);

    Graph g(g_size, false);
    vector<vector<int>> edges(g_size);
    vector<vector<int>> weights(g_size);

    vector<vector<bool>> edge_exists(g_size, vector<bool>(g_size, false));
    int u, v, w;
    while (inputFile >> u >> v >> w) {
        if (u > g_size || v > g_size) {
            cout << "Error: vertex out of bounds" << endl;
            exit(1);
        }
        if (!g.containsVertex[u])
            g.addVertex(u);

        if (!g.containsVertex[v])
            g.addVertex(v);

        if (!edge_exists[u][v]) {
            edges[u].push_back(v);
            weights[u].push_back(w);
            edge_exists[u][v] = true;
        }
    }
    for (int i = 0; i < g_size; i++)
        g.addEdges(i, edges[i], weights[i]);
    return g;
}

Graph readInput(ifstream &inputFile) {
    int g_size;
    int maxWeight;
    inputFile >> g_size; // >> maxWeight;
    ++g_size;
//    vector<int> phi = createRandomPriceFunction(g_size, maxWeight);

    Graph g(g_size, false);
    vector<vector<int>> edges(g_size);
    vector<vector<int>> weights(g_size);

    vector<vector<bool>> edge_exists(g_size, vector<bool>(g_size, false));
    int u, v, w;
    while (inputFile >> u >> v >> w) {
        if (u > g_size || v > g_size) {
            exit(1);
        }
        if (!g.containsVertex[u])
            g.addVertex(u);

        if (!g.containsVertex[v])
            g.addVertex(v);

        if (!edge_exists[u][v]) {
            edges[u].push_back(v);
            weights[u].push_back(w);
            edge_exists[u][v] = true;
        }
    }

    if (MAKE_CONNECTED) {
        for (int i: g.vertices) {
            edges[0].push_back(i);
            weights[0].push_back(0);

            edges[i].push_back(0);
            weights[i].push_back(0);
        }
        g.addVertex(0);
    }

    for (int i = 0; i < g_size; ++i)
        g.addEdges(i, edges[i], weights[i]);
    return g;
}

vector<int> createRandomPriceFunction(int g_size, int maxWeight) {
    vector<int> phi(g_size);
    for (int i = 0; i < g_size; ++i)
        phi[i] = Random::Get().GenInt(0, maxWeight);
    return phi;
}

int getWeight(int weight, int u, int v, vector<int> &phi) {
    return weight + phi[u] - phi[v];
}

Graph getConnectedSubgraph(Graph &g) {
    vector<bool> reachable(g.v_max, false);
    findReachable(g, SRC, reachable);

    Graph subGraph = Graph(g.v_max, false);
    for (int v = 0; v < g.v_max; ++v)
        if (reachable[v])
            subGraph.addVertex(v);

    for (int u: subGraph.vertices) {
        vector<int> outVertices;
        vector<int> weights;

        for (int i = 0; i < g.adjacencyList[u].size(); ++i) {
            if (reachable[g.adjacencyList[u][i]]) {
                outVertices.push_back(g.adjacencyList[u][i]);
                weights.push_back(g.weights[u][i]);
            }
        }

        subGraph.addEdges(u, outVertices, weights);
    }

    return subGraph;
}

void findReachable(Graph &g, int s, vector<bool> &reachable) {
    reachable[s] = true;
    for (int v: g.adjacencyList[s])
        if (!reachable[v])
            findReachable(g, v, reachable);
}

vector<int> bitScaling(Graph &g) {
    int LDD_BASE_CASE = 10;
    int CALCULATE_SCC_PROB = 1;
    LDD_BASE_CASE = g.n / (LDD_BASE_CASE * log(g.n));
    CALCULATE_SCC_PROB = g.n / 10000.0;

    int minWeight = INT_MAX;
    for (int u: g.vertices)
        for (int i = 0; i < g.adjacencyList[u].size(); ++i)
            minWeight = min(minWeight, g.weights[u][i]);

    if (minWeight >= 0) {
         cout << "Graph is non-negative" << endl;
        return getShortestPathTree(g, SRC);
    }

    int precision = int(pow(2, int(logBase2(-1 * minWeight))));
    vector<int> phi(g.v_max);

    while (precision >= 1) {
        Graph gScaledS(g.v_max + 1, false);
        gScaledS.addVertices(g.vertices);
        gScaledS.addVertex(g.v_max);

        for (int u: g.vertices) {
            int numOutEdges = g.adjacencyList[u].size();
            vector<int> edges(numOutEdges);
            vector<int> weights(numOutEdges);

            for (int i = 0; i < numOutEdges; ++i) {
                int roundedWeight = phi[u] - phi[g.adjacencyList[u][i]] + ceil(g.weights[u][i] / (double) precision);
                if (roundedWeight < -1)
                    throw_with_nested("Bit scaling produced an edge with weight less than -1");

                edges[i] = g.adjacencyList[u][i];
                weights[i] = roundedWeight;
            }
            gScaledS.addEdges(u, edges, weights);
        }

        vector<int> dummyEdges(g.n);
        vector<int> dummyWeights(g.n);
        for (int i = 0; i < g.n; ++i) {
            dummyEdges[i] = g.vertices[i];
            dummyWeights[i] = 0;
        }
        gScaledS.addEdges(g.v_max, dummyEdges, dummyWeights);

        vector<int> tree = SPMain(gScaledS, g.v_max);
        vector<int> dist(g.v_max + 1);

        vector<bool> vectorTmp(g.v_max + 1, false);

        getDistances(gScaledS, tree, dist, 0, g.v_max, vectorTmp);

        if (CHECKS) {
            verifyTree(gScaledS, tree, dist, g.v_max);

            if (hasNegativeEdges(gScaledS, dist, 0))
                throw_with_nested("Bit scaling produced a graph with negative edges");
        }

        for (int u: g.vertices)
            if (precision == 1)
                phi[u] += dist[u];
            else
                phi[u] = 2 * (phi[u] + dist[u]);

        precision /= 2;
    }

    Graph gFinal(g);

    for (int u: gFinal.vertices) {
        for (int i = 0; i < gFinal.adjacencyList[u].size(); ++i) {
            gFinal.weights[u][i] += phi[u] - phi[gFinal.adjacencyList[u][i]];

            if (gFinal.weights[u][i] < 0) {
                throw_with_nested("Apply phi outputted - existing negative edge");
            }
        }
    }

    vector<int> tree = getShortestPathTree(gFinal, SRC);

    return tree;
}

void verifyTree(Graph &g, vector<int> &tree, vector<int> &dist, int src) {
    vector<vector<bool>> adjList(g.v_max, vector<bool>(g.v_max, false));
    for (int u: g.vertices)
        if (tree[u] != -1) {
            if (!g.containsVertex[tree[u]])
                throw_with_nested("Tree contains a vertex that does not exist");
            adjList[tree[u]][u] = true;
        }
    vector<bool> visited(g.v_max, false);
    if (containsCycles(g, adjList, src, visited))
        throw_with_nested("Tree contains a cycle");
    for (int u: g.vertices) {
        for (int i = 0; i < g.adjacencyList[u].size(); ++i) {
            int v = g.adjacencyList[u][i];
            if (dist[v] > dist[u] + g.weights[u][i])
                throw_with_nested("SPMain returned a tree that is not a shortest paths tree.");
        }
    }
}

bool containsCycles(Graph &g, vector<vector<bool>> &adjList, int src, vector<bool> &visited) {
    if (visited[src])
        return true;
    visited[src] = true;
    for (int v = 0; v < g.v_max; ++v)
        if (adjList[src][v])
            if (containsCycles(g, adjList, v, visited))
                return true;
    return false;
}

void getDistances(Graph &g, vector<int> &tree, vector<int> &dist, int curDis, int curVertex,
                  vector<bool> &visited) {
    if (!visited[curVertex]) {
        dist[curVertex] = curDis;
        visited[curVertex] = true;

        for (int i = 0; i < g.adjacencyList[curVertex].size(); i++) {
            int v = g.adjacencyList[curVertex][i];
            if (tree[v] == curVertex) {
                getDistances(g, tree, dist, curDis + g.weights[curVertex][i], v, visited);
            }
        }
    }
}

vector<int> SPMain(Graph &g_in, int s) {
    cout << "SPMain start at: " << Timer::getDuration() << endl;
    int scaleFactor = 2 * g_in.n;
    Graph g = getScaledGraph(g_in, scaleFactor);
    int B = roundPower2(scaleFactor);
    vector<int> phi(g.v_max);
    set<vector<int>> nullErem;
    for (int i = 1; i <= logBase2(B); i++) {
        nullErem.clear();
        Graph g_phi = createModifiedGB(g, 0, false, nullErem, phi);
        vector<int> phi_i = ScaleDown(g_phi, g.n, B / (int) pow(2, i));

        if (CHECKS && hasNegativeEdges(g_phi, phi_i, B / (int) pow(2, i))) {
            throw_with_nested("ScaleDown failed.");
        }

        phi = addPhi(phi, phi_i);
    }

    // create G^*
    for (int u: g.vertices) {
        for (int i = 0; i < g.adjacencyList[u].size(); i++) {
            g.weights[u][i] += phi[u] - phi[g.adjacencyList[u][i]] + 1;

            if (CHECKS && g.weights[u][i] < 0) {
                throw_with_nested(
                        "After applying the phi outputted from SPMain, there exists an edge with negative weight.");
            }
        }
    }

    vector<int> tree = getShortestPathTree(g, s);

    if (CHECKS && invalidTree(g, s, tree)) {
        throw_with_nested("SPMain get shortest path tree failed.");
    }
    cout << "SPMain end at: " << Timer::getDuration() << endl;

    return tree;
}

Graph getScaledGraph(Graph &g_in, int scaleFactor) {
    Graph g = Graph(g_in.v_max, false);
    g.addVertices(g_in.vertices);

    for (int u: g_in.vertices) {
        vector<int> edges(g_in.adjacencyList[u].size());
        vector<int> weights(g_in.adjacencyList[u].size());

        for (int i = 0; i < g_in.adjacencyList[u].size(); ++i) {
            edges[i] = g_in.adjacencyList[u][i];
            weights[i] = g_in.weights[u][i] * scaleFactor;
        }
        g.addEdges(u, edges, weights);
    }

    return g;
}

bool invalidTree(Graph &g, int s, vector<int> &tree) {
    for (int u = 0; u < tree.size(); ++u) {
        if (g.containsVertex[u]) {
            if (u != s && tree[u] == -1)
                return true;
            if (u == s && tree[u] != -1)
                return true;
        }
    }
    return false;
}

int roundPower2(int n) {
    return int(pow(2, ceil(logBase2(n))));
}

double logBase2(int n) {
    return log(n) / log(2);
}

/*
 * 1. INPUT REQUIREMENTS:
 * (a) B is positive integer, w is integral, and w(e) ≥ −2B for all e ∈ E
 * (b) If the graph G does not contain a negative-weight cycle then the input
 * must satisfy η(GB) ≤ ∆; that is, for every v ∈ V there is a shortest
 * sv-path in GBs with at most ∆ negative edges
 * (c) All vertices in G have constant out-degree
 * 2. OUTPUT: If it terminates, the algorithm returns an integral price function
 * φ
 * such that wφ(e) ≥ −B for all e ∈ E
 * 3. RUNNING TIME: If G does not contain a negative-weight cycle, then the
 * algorithm has epected runtime O(m log3(n) log(∆)).
 * Remark: If G contains a negative-weight cycle, there is no guarantee
 * on the runtime, and the algorithm might not even terminate; but if the
 * algorithm does terminate, it always produces a correct output.
 */
vector<int> ScaleDown(Graph &g, int delta, int B) {
    vector<int> phi_2(g.v_max);
    vector<int> emptyPhi;
    set<vector<int>> emptyRem;

    if (delta > 2) {
        double d = delta / 2.0;
        Graph g_b_nneg = createModifiedGB(g, B, true, emptyRem, emptyPhi);

        // phase 0
        vector<vector<int>> E_sep;
        if (WITH_LDD)
            E_sep = SPmainLDD(g_b_nneg, int(4 * d * B));

        set<vector<int>> E_sep_hash(E_sep.begin(), E_sep.end());
        Graph g_B_Esep = createModifiedGB(g, B, false, E_sep_hash, emptyPhi);
        vector<vector<int>> SCCs = g_B_Esep.SCC();

        // phase 1
        vector<int> vertexToSCCMap = getVertexToSCCMap(SCCs, g.v_max);
        set<vector<int>> edgesBetweenSCCs = getEdgesBetweenSCCs(g, vertexToSCCMap);
        Graph H = createModifiedGB(g, 0, false, edgesBetweenSCCs, emptyPhi);
        vector<int> phi_1 = ScaleDown(H, delta / 2, B);

        // phase 2
        Graph g_B_E_sep_phi1 = createModifiedGB(g, B, false, E_sep_hash, phi_1);
        vector<int> phi = FixDAGEdges(g_B_E_sep_phi1, SCCs, vertexToSCCMap, edgesBetweenSCCs);
        phi_2 = addPhi(phi_1, phi);

        if (CHECKS && hasNegativeEdges(g_B_Esep, phi_2, 0))
            throw_with_nested("FixDAGEdges failed.");
    }

    // phase 3
    Graph g_B_phi2 = createModifiedGB(g, B, false, emptyRem, phi_2);
    vector<int> phi_prime = ElimNeg(g_B_phi2);
    vector<int> phi_3 = addPhi(phi_2, phi_prime);

    if (CHECKS && hasNegativeEdges(g_B_phi2, phi_prime, 0))
        throw_with_nested("ElimNeg failed.");

    return phi_3;
}

vector<vector<int>> SPmainLDD(Graph &g, int diameter) {
    vector<vector<int>> E_sep;

    // first remove all the large edges in G

    Graph largeEdgesRemoved(g.v_max, false);
    largeEdgesRemoved.addVertices(g.vertices);

    for (int v: g.vertices) {
        vector<int> outVertices;
        vector<int> weights;

        for (int i = 0; i < g.adjacencyList[v].size(); i++) {
            if (g.weights[v][i] > diameter) {
                // edge is too big, can add to E_sep
                vector<int> edge = {v, g.adjacencyList[v][i]};
                E_sep.push_back(edge);
            } else {
                outVertices.push_back(g.adjacencyList[v][i]);
                weights.push_back(g.weights[v][i]);
            }
        }

        largeEdgesRemoved.addEdges(v, outVertices, weights);
    }
    vector<vector<int>> LDD = preLDD(largeEdgesRemoved, diameter);
    E_sep.insert(E_sep.end(), LDD.begin(), LDD.end());

    return E_sep;
}

bool hasNegativeEdges(Graph &g, vector<int> &phi, int B) {
    for (int u: g.vertices) {
        for (int i = 0; i < g.adjacencyList[u].size(); i++) {
            if (g.weights[u][i] + phi[u] - phi[g.adjacencyList[u][i]] < -1 * B) {
                return true;
            }
        }
    }

    return false;
}

/*
 * Creates G^B_phi = (V, E, w^B_phi), where
 * w^B_phi(e) = w(e) + phi(u) - phi(v) if w(e) >= 0,
 * and w(e) + B + phi(u) - phi(v) if w(e) < 0.
 * If nneg == true, w(e) = max{0, w^B(e)}.
 * Removes all the edges in remEdges.
 */
Graph createModifiedGB(Graph &g, int B, bool nneg, set<vector<int>> &remEdges, vector<int> &phi) {
    Graph modG(g.v_max, false);
    modG.addVertices(g.vertices);
    vector<int> edges, weights;

    for (int u: g.vertices) {
        edges = vector<int>(g.adjacencyList[u].size());
        weights = vector<int>(g.adjacencyList[u].size());

        for (int i = 0; i < g.adjacencyList[u].size(); i++) {
            int v = g.adjacencyList[u][i];

            vector<int> edge = {u, v};
            if ((remEdges.empty()) || (remEdges.find(edge) != remEdges.end())) {
                int weight = g.weights[u][i];

                if (weight < 0) {
                    weight += B;
                }
                if (nneg) {
                    weight = max(0, weight);
                }
                if (!phi.empty()) {
                    weight += phi[u] - phi[v];
                }

                edges[i] = v;
                weights[i] = weight;
            }
        }

        modG.addEdges(u, edges, weights);
    }

    return modG;
}

set<vector<int>> getEdgesBetweenSCCs(Graph &g, vector<int> &vertexToSCCMap) {
    set<vector<int>> edgesBetweenSCCs;
    for (int u: g.vertices) {
        for (int v: g.adjacencyList[u]) {
            if (vertexToSCCMap[u] != vertexToSCCMap[v]) {
                vector<int> edge = {u, v};
                edgesBetweenSCCs.insert(edge);
            }
        }
    }

    return edgesBetweenSCCs;
}

vector<int> getVertexToSCCMap(vector<vector<int>> &SCCs, int numVertices) {
    vector<int> vertexToSCCMap(numVertices, 0);
    for (int i = 0; i < SCCs.size(); i++) {
        for (int v: SCCs[i]) {
            vertexToSCCMap[v] = i;
        }
    }

    return vertexToSCCMap;
}

vector<int> addPhi(vector<int> &phi_1, vector<int> &phi_2) {
    if (phi_1.size() != phi_2.size())
        throw_with_nested("addPhi: phi_1 and phi_2 must have the same length.");

    vector<int> phi(phi_1.size());
    for (int i = 0; i < phi_1.size(); i++) {
        phi[i] = phi_1[i] + phi_2[i];
    }
    return phi;
}

vector<int> FixDAGEdges(Graph &g,
                        vector<vector<int>> &SCCs,
                        vector<int> &vertexToSCCMap,
                        set<vector<int>> &edgesBetweenSCCs) {
    int n = SCCs.size();
    vector<vector<int>> SCCAdjList = createSCCAdjList(SCCs, vertexToSCCMap, edgesBetweenSCCs);
    vector<int> topOrdering = topSort(n, SCCAdjList);

    vector<int> mu(n, 0); // indices are in topological order (e.g., index 0 corresponds to the first SCC
    // in topological order)
    for (int u: g.vertices) {
        for (int i = 0; i < g.adjacencyList[u].size(); i++) {
            int v = g.adjacencyList[u][i];

            int SCCu = vertexToSCCMap[u];
            int SCCv = vertexToSCCMap[v];
            int edgeWeight = g.weights[u][i];

            if ((SCCu != SCCv) && edgeWeight < mu[topOrdering[SCCv]]) {
                mu[topOrdering[SCCv]] = edgeWeight;
            }
        }
    }

    vector<int> phi(g.v_max);
    int m = 0;
    for (int j = 1; j < n; j++) {
        m += mu[j];
        for (int v: SCCs[topOrdering[j + n]]) {
            phi[v] = m;
        }
    }

    return phi;
}

// returns the adjacency list for the DAG where every SCC is viewed as a single vertex
vector<vector<int>> createSCCAdjList(vector<vector<int>> &SCCs,
                                     vector<int> &vertexToSCCMap,
                                     set<vector<int>> &edgesBetweenSCCs) {
    vector<vector<int>> SCCAdjList(SCCs.size());
    vector<vector<bool>> containsEdge(SCCs.size(), vector<bool>(SCCs.size(), false));

    for (vector<int> edge: edgesBetweenSCCs) {
        int u = vertexToSCCMap[edge[0]];
        int v = vertexToSCCMap[edge[1]];
        if (!containsEdge[u][v]) {
            containsEdge[u][v] = true;
            SCCAdjList[u].push_back(v);
        }
    }

    return SCCAdjList;
}

/*
 * Input: DAG
 * Calculates a topological ordering of the n vertices in the DAG
 * Returns a map of vertex v to its index i in the ordering
 * The index in the order i is also mapped to the vertex v (bijective mapping);
 * however, to prevent duplicate keys, instead of using the key i we use the
 * key i + n.
 */
vector<int> topSort(int n, vector<vector<int>> &adjList) {
    vector<int> topOrdering(2 * n);
    stack<int> stack;
    vector<bool> visited(n, false);

    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            topSortUtil(i, visited, stack, adjList);
        }
    }

    int i = 0;
    while (!stack.empty()) {
        int v = stack.top();
        stack.pop();
        topOrdering[v] = i;
        topOrdering[i + n] = v;
        i++;
    }

    return topOrdering;
}

void topSortUtil(int u, vector<bool> &visited, stack<int> &stack, vector<vector<int>> &adjList) {
    visited[u] = true;

    for (int v: adjList[u]) {
        if (!visited[v]) {
            topSortUtil(v, visited, stack, adjList);
        }
    }
    stack.push(u);
}

/*
 * ElimNeg takes as input a graph G = (V,E,w) in which all vertices have
 * constant out-degree.
 * The algorithm outputs a price function φ such that w_φ(e) ≥ 0 for all e ∈ E
 * and has running time O(log(n)*(n + sum_v∈V ηG(v)));
 * note that if G contains a negative-weight cycle then v∈V ηG(v) = ∞
 * so the algorithm will never terminate and hence not produce any output.
 */
vector<int> ElimNeg(Graph &g) {
    Graph Gs = createGs(g);
    vector<int> dist(Gs.v_max);
    int s = Gs.v_max - 1;

    dist[s] = 0;
    for (int v = 0; v < s; v++) {
        dist[v] = INT_MAX;
    }

    custom_priority_queue<Node> pq;
    vector<bool> inPQ(Gs.v_max);
    pq.emplace(s, dist[s]);
    inPQ[s] = true;

    set<int> marked;
//    map<int, bool> tobeRemoved;

    while (!pq.empty()) {
        // Dijkstra Phase
        while (!pq.empty()) {
            int v = pq.top().node;
            int w = pq.top().cost;
            pq.pop();
            if (w > dist[v]) {
                continue;
            }


            inPQ[v] = false;

            marked.emplace(v);

            for (int i = 0; i < Gs.adjacencyList[v].size(); i++) {
                int x = Gs.adjacencyList[v][i];
                int edgeWeight = Gs.weights[v][i];

                if (edgeWeight >= 0 && (dist[v] + edgeWeight < dist[x])) {
                    marked.emplace(x);

                    dist[x] = dist[v] + edgeWeight;
                    pq.push(Node(x, dist[x]));
                    inPQ[x] = true;
                }
            }
        }

        // Bellman-Ford Phase
        for (int v: marked) {
            for (int i = 0; i < Gs.adjacencyList[v].size(); i++) {
                int x = Gs.adjacencyList[v][i];
                int edgeWeight = Gs.weights[v][i];

                if (edgeWeight < 0 && (dist[v] + edgeWeight < dist[x])) {
                    dist[x] = dist[v] + edgeWeight;
                    pq.push(Node(x, dist[x]));
                    inPQ[x] = true;
                }
            }
        }
        marked.clear();
    }

    vector<int> phi(g.v_max);
    for (int v: g.vertices) {
        phi[v] = dist[v];
    }
    return phi;
}

/*
 * Returns the graph that is g with an added dummy vertex s
 * and edges of weight 0 connecting s to every vertex in G.
 */
Graph createGs(Graph &g) {
    int s = g.v_max;
    Graph Gs(g.v_max + 1, false);
    Gs.addVertices(g.vertices);
    Gs.addVertex(s);

    for (int u: g.vertices) {
        Gs.addEdges(u, g.adjacencyList[u], g.weights[u]);
    }

    vector<int> edges(g.v_max);
    vector<int> weights(g.v_max);
    for (int i = 0; i < g.v_max; i++) {
        edges[i] = i;
        weights[i] = 0;
    }
    Gs.addEdges(s, edges, weights);

    return Gs;
}

vector<int> getShortestPathTree(Graph &g, int s) {
    set<int> settled;
    priority_queue<Node> pq;
    vector<int> dist(g.v_max);
    vector<int> tree(g.v_max);

    for (int i = 0; i < g.v_max; i++) {
        dist[i] = INT_MAX;
        tree[i] = -1;
    }

    pq.emplace(s, 0);
    dist[s] = 0;

    while (settled.size() != g.n) {
        if (pq.empty()) {
            return tree;
        }
        int u = pq.top().node;
        pq.pop();

        if (settled.find(u) != settled.end()) {
            continue;
        }

        settled.emplace(u);
        updateTreeNeighbors(g, u, tree, settled, pq, dist);
    }

    return tree;
}

void updateTreeNeighbors(Graph &g, int u, vector<int> &tree,
                         set<int> &settled, priority_queue<Node> &pq,
                         vector<int> &dist) {
    for (int i = 0; i < g.adjacencyList[u].size(); i++) {
        int v = g.adjacencyList[u][i];
        if (settled.find(v) == settled.end()) {
            int weight = g.weights[u][i];
            int newDistance = dist[u] + weight;
            if (newDistance < dist[v]) {
                dist[v] = newDistance;
                tree[v] = u;
            }

            pq.emplace(v, dist[v]);
        }
    }
}

vector<int> bellmanFord(Graph &g) {
    vector<int> dist(g.v_max);
    for (int i = 0; i < g.v_max; i++) {
        dist[i] = INT_MAX;
    }
    dist[SRC] = 0;
    for (int i = 1; i < g.v_max; i++) {
        for (int u: g.vertices) {
            for (int j = 0; j < g.adjacencyList[u].size(); j++) {
                int v = g.adjacencyList[u][j];
                int weight = g.weights[u][j];

                if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                }
            }
        }
    }

    for (int u: g.vertices) {
        for (int j = 0; j < g.adjacencyList[u].size(); j++) {
            int v = g.adjacencyList[u][j];
            int weight = g.weights[u][j];

            if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
                throw_with_nested("Negative cycle detected");
            }
        }
    }

    return dist;
}

vector<int> getDistFromTree(Graph &g, vector<int> &tree) {
    vector<int> dist(g.v_max);
    for (int i = 0; i < g.v_max; i++) {
        dist[i] = INT_MAX;
    }
    dist[SRC] = 0;

    updateDistFromTree(g, tree, dist, SRC);

    return dist;
}

void updateDistFromTree(Graph &g, vector<int> &tree, vector<int> &dist, int u) {
    for (int i = 0; i < g.adjacencyList[u].size(); i++) {
        int v = g.adjacencyList[u][i];
        int weight = g.weights[u][i];

        if (tree[v] == u) {
            dist[v] = dist[u] + weight;
            updateDistFromTree(g, tree, dist, v);
        }
    }
}