#include<iostream>
#include<fstream>
#include<vector>


using namespace std;

typedef vector<vector<int>> VV;
typedef pair<int, int> P;

class BiMatrix {
    int n;
    VV main_matrix;
    VV det_matrix;

public:
    BiMatrix(const VV &main_matrix) {
        this->n = main_matrix.size();
        this->main_matrix = main_matrix;
        this->det_matrix = VV(n - 1, vector<int>(n - 1, 1));
    }

    double calc_minor_det(P c1, unsigned int n) {

    }
};





void print_matrix(VV vv) {
    for (int i = 0; i < vv.size(); i++) {
        for (int j = 0; j < vv.size(); j++) {
            cout << vv[i][j] << ' ';
        }
        cout << '\n';
    }
}


int main() {
    ifstream in("1.txt");
    int n = 0;
    in >> n;
    VV vv(n, vector<int>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            in >> vv[i][j];
        }
    }
    print_matrix(vv);
}
