#include "LasVegas.h"

int main() {
    string filename;
    cout << "Enter filename (eg: \"../graph/500_2.txt\"): ";
    cin >> filename;
    cout << "Enter source vertex: ";
    cin >> SRC;
    cout << "Enter 1 for LDD, 0 for no LDD: ";
    cin >> WITH_LDD;
    Random::Get().Seed();
    ifstream inputFile(filename);
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
            cout << i << ' ' << BellManFord[i] << " " << LasVegas[i] << endl;
            return 0;
        }
    }
    cout << "Bellman-Ford and Las Vegas are equal" << endl;
    return 0;
}
