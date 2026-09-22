#include<iostream>
#include<fstream>
#include<vector>


using namespace std;

typedef vector<vector<int>> VV;
typedef pair<int, int> P;


void print_matrix(VV vv) {
    for (int i = 0; i < vv.size(); i++) {
        for (int j = 0; j < vv.size(); j++) {
            cout << vv[i][j] << ' ';
        }
        cout << '\n';
    }
}


class BiMatrix {
    int n;
    VV origin_matrix;
    VV main_matrix;
    VV det_matrix;

public:
    BiMatrix(const VV &origin_matrix): origin_matrix(origin_matrix) {
        this->n = this->main_matrix.size();
        this->main_matrix = this->origin_matrix;
        this->det_matrix = VV(n - 1, vector<int>(n - 1, 1));
    }

    // 0-index
    int calc_2minor_det(P c1) {
        int a = main_matrix[c1.first][c1.second];
        int b = main_matrix[c1.first][c1.second + 1];
        int c = main_matrix[c1.first + 1][c1.second];
        int d = main_matrix[c1.first + 1][c1.second + 1];
        return a * d - b * c;
    }

    double calc_Dodgson_det() {

    }
};


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
