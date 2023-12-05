#include "LasVegas.h"

int main() {
    Random::Get().Seed();
    ifstream inputFile("../graph/002.txt");
    Graph g = readInput(inputFile);
    vector<int> BellManFord;
    try {
        BellManFord = bellmanFord(g);
    } catch (const char *msg) {
        cout << msg << endl;
        return 0;
    }
    vector<int> LasVegas;
    Timer::startTimer();
    try {
        LasVegas = bitScaling(g);
    }
    catch (const char *msg) {
        cout << msg << endl;
        return 0;
    }
    for (int LasVega : LasVegas)
        cout << LasVega << " ";
    cout << endl;
    vector<int> LasVegasDist = getDistFromTree(g, LasVegas);

    for (int i = 0; i < BellManFord.size(); i++) {
        if (BellManFord[i] != LasVegasDist[i]) {
            cout << "Bellman-Ford and Las Vegas are not equal" << endl;
            return 0;
        }
    }
    cout << "Bellman-Ford and Las Vegas are equal" << endl;
    return 0;
}
