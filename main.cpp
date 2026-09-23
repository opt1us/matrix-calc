#include<iostream>
#include<fstream>
#include<vector>
#include<iomanip>


using namespace std;

typedef vector<vector<int>> VV;


void print_matrix(VV vv) {
    for (int i = 0; i < vv.size(); i++) {
        for (int j = 0; j < vv.size(); j++) {
            cout << setw(3) << vv[i][j] << ' ';
        }
        cout << '\n';
    }
    cout << '\n';
}


class BiMatrix {
    int n;
    VV origin_matrix;
    VV main_matrix;
    VV det_matrix;

public:
    BiMatrix(const VV &origin_matrix): origin_matrix(origin_matrix) {
        this->n = this->origin_matrix.size();
        this->main_matrix = this->origin_matrix;
        this->det_matrix = VV(n - 1, vector<int>(n - 1, 1));
    }

    VV getDetMatrix(){
        return this->det_matrix;
    }

    VV getMainMatrix(){
        return this->main_matrix;
    }

    // 0-index
    int calc_2minor_det(int x, int y) {
        int a = main_matrix[x][y];
        int b = main_matrix[x][y + 1];
        int c = main_matrix[x + 1][y];
        int d = main_matrix[x + 1][y + 1];
        return a * d - b * c;
    }

    int calc_Dodgson_det() {
        while(det_matrix.size() > 0){
            print_BiMatrix();
            VV tmp_main(n-1, vector<int>(n-1, 0));
            for (int i = 0; i < n-1; i++) {
                for (int j = 0; j < n-1; j++) {
                    tmp_main[i][j] = (calc_2minor_det(i,j) / det_matrix[i][j]);
                }
            }

            VV tmp_det(n-2, vector<int>(n-2, 0));
            for (int i = 1; i < n-1; i++) {
                for (int j = 1; j < n-1; j++) {
                    tmp_det[i-1][j-1] = main_matrix[i][j];
                }
            }
            main_matrix = tmp_main;
            det_matrix = tmp_det;
            n--;
        }
        return main_matrix[0][0];
    }

    void print_BiMatrix() {
        cout << "bi\n";
        for (int i = 0; i < main_matrix.size() - 1; i++) {
            for (int j = 0; j < main_matrix.size(); j++) {
                cout << setw(5) << main_matrix[i][j] << ' ';
            }
            cout << "\033[35m";
            cout << '\n' << "   ";
            for (int j = 0; j < det_matrix.size(); j++) {
                cout << setw(5) << det_matrix[i][j] << ' ';
            }
            cout << "\033[37m";
            cout << '\n';
        }
        for (int j = 0; j < main_matrix.size(); j++) {
            cout << setw(5) << main_matrix[main_matrix.size() - 1][j] << ' ';
        }
        cout << "\n\n";
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
    BiMatrix bm = BiMatrix(vv);
    cout << bm.calc_Dodgson_det() << '\n';

    cout << "------------------\n";
    print_matrix(bm.getMainMatrix());

}
