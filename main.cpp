#include "LasVegas.h"

int main() {
    Random::Get().Seed();
    ifstream inputFile("../graph/2500_1.txt");
    Graph g = readInput(inputFile);

    Timer::startTimer();
    preLDD(g, 50);
    cout << "LDD: " << Timer::getDuration() << " seconds" << endl;

//    vector<int> BellManFord;
//    try {
//        Timer::startTimer();
//        BellManFord = bellmanFord(g);
//        cout << "Bellman-Ford took " << Timer::getDuration() << " seconds" << endl;
//    } catch (const char *msg) {
//        cout << msg << endl;
//        return 0;
//    }
//    vector<int> LasVegas;
//    try {
//        Timer::startTimer();
//        LasVegas = bitScaling(g);
//        cout << "Las Vegas took " << Timer::getDuration() << " seconds" << endl;
//    }
//    catch (const char *msg) {
//        cout << msg << endl;
//        return 0;
//    }
//    vector<int> LasVegasDist = getDistFromTree(g, LasVegas);
//
//    for (int i = 0; i < BellManFord.size(); i++) {
//        if (BellManFord[i] != LasVegasDist[i]) {
//            cout << "Bellman-Ford and Las Vegas are not equal" << endl;
//            return 0;
//        }
//    }
//    cout << "Bellman-Ford and Las Vegas are equal" << endl;
    return 0;
}
