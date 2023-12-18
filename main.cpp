#include "LasVegas.h"

int main(int argc, char *argv[]) {
    string filename = "../graph/5000_25.txt";
    WITH_LDD = true;
//    string filename = argv[1];
//    if (argc > 2)
//        SRC = stoi(argv[3]);
//    if (argc == 4)
//        WITH_LDD = stoi(argv[2]);
    Random::Get().Seed();
    ifstream inputFile(filename);
    Graph g = readInput(inputFile);
    vector<int> BellmanFord = bellmanFord(g);

//    vector<int> LasVegas = bitScaling(g);
    vector<int> LasVegasReal = lasVegas(g);
    for (unsigned long i = 0; i < BellmanFord.size(); i++) {
//        cout << BellmanFord[i] << " " << LasVegasReal[i] << endl;
        if (LasVegasReal[i] != BellmanFord[i]) {
            cout << "Failed: Bellman-Ford and Las Vegas are not equal" << endl;
            return 0;
        }
    }
    return 0;
}
