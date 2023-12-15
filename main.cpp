#include "LasVegas.h"

int main(int argc, char *argv[]) {
    string filename = argv[1];
    WITH_LDD = stoi(argv[2]);
    if (argc == 4) {
        SRC = stoi(argv[3]);
    } else {
        SRC = 0;
    }
//    string filename = "../graph/5000_25.txt";
//    WITH_LDD = 1;
//    SRC = 0;
    Random::Get().Seed();
    ifstream inputFile(filename);
    Graph g = readInput(inputFile);
    vector<int> BellManFord;
    try {
        BellManFord = bellmanFord(g);
    } catch (const char *msg) {
        cout << "Failed: " << msg << endl;
        return 0;
    }
    vector<int> LasVegas;
    try {
        LasVegas = bitScaling(g);
    }
    catch (const char *msg) {
        cout << "Failed: "  << msg << endl;
        return 0;
    }

    for (unsigned long i = 0; i < BellManFord.size(); i++) {
        if (BellManFord[i] != LasVegas[i]) {
            cout << "Failed: Bellman-Ford and Las Vegas are not equal" << endl;
            return 0;
        }
    }
    cout << "Succeeded: Bellman-Ford and Las Vegas are equal" << endl;
    return 0;
}
