#include "LasVegas.h"

int main() {
    Random::Get().Seed();
    ifstream inputFile("../graph/2500_1.txt");
    Graph g = readInput(inputFile);
    vector<int> BellManFord;
    try {
        BellManFord = bellmanFord(g);
    } catch (const char *msg) {
        cout << msg << endl;
        return 0;
    }
    vector<int> LasVegas;
    try {
        LasVegas = bitScaling(g);
    }
    catch (const char *msg) {
        cout << msg << endl;
        return 0;
    }

    for (int i = 0; i < BellManFord.size(); i++) {
        if (BellManFord[i] != LasVegas[i]) {
            cout << "Bellman-Ford and Las Vegas are not equal" << endl;
            return 0;
        }
    }
    cout << "Bellman-Ford and Las Vegas are equal" << endl;
    return 0;
}
