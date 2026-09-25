#include<iostream>
#include<fstream>
#include<vector>
#include<iomanip>
#include<math.h>


using namespace std;

template<typename T>
using VVt = vector<vector<T>>;

using VV = VVt<double>;

class Matrix {
    VV origin_matrix;
    VV matrix;

public:
    Matrix(const VV &matrix): origin_matrix(matrix), matrix(matrix) {}

    VV getMatrix() {
        return this->matrix;
    }

    int getSize() {
        return this->matrix.size();
    }

    // 0-index
    auto calc_2minor_det(int row, int col) {
        auto a = matrix[row][col];
        auto b = matrix[row][col + 1];
        auto c = matrix[row + 1][col];
        auto d = matrix[row + 1][col + 1];
        return a * d - b * c;
    }

    // 0-index
    auto calc_2minor_det(int r1, int r2, int c1, int c2) {
        auto a = matrix[r1][c1];
        auto b = matrix[r1][c2];
        auto c = matrix[r2][c1];
        auto d = matrix[r2][c2];
        return a * d - b * c;
    }

    // 0-index
    vector<double> operator[](int index) const {
        return this->matrix[index];
    }

    void print_matrix() {
        for (int i = 0; i < matrix.size(); i++) {
            for (int j = 0; j < matrix.size(); j++) {
                cout << setw(3) << matrix[i][j] << ' ';
            }
            cout << '\n';
        }
        cout << '\n';
    }

};


class BiMatrix {
    VV origin_matrix;
    VV main_matrix;
    VV det_matrix;

public:
    BiMatrix(const VV &origin_matrix) {
        this->origin_matrix = origin_matrix;
        this->main_matrix = origin_matrix;
        this->det_matrix = VV(this->main_matrix.size() - 1, vector<double>(this->main_matrix.size() - 1, 1));
    }

    BiMatrix(const VV &main_matrix, const VV &det_matrix):
        origin_matrix(main_matrix), main_matrix(main_matrix), det_matrix(det_matrix) {}

    VV getDetMatrix(){
        return this->det_matrix;
    }

    VV getMainMatrix(){
        return this->main_matrix;
    }

    int getSize() {
        return this->main_matrix.size();
    }

    void restore() {
        this->main_matrix = origin_matrix;
        this->det_matrix = VV(main_matrix.size() - 1, vector<double>(main_matrix.size() - 1, 1));
    }

    // 0-index
    auto calc_2minor_det(int row, int col) {
        auto a = main_matrix[row][col];
        auto b = main_matrix[row][col + 1];
        auto c = main_matrix[row + 1][col];
        auto d = main_matrix[row + 1][col + 1];
        return a * d - b * c;
    }

    // 0-index
    auto calc_2minor_det(int r1, int r2, int c1, int c2) {
        auto a = main_matrix[r1][c1];
        auto b = main_matrix[r1][c2];
        auto c = main_matrix[r2][c1];
        auto d = main_matrix[r2][c2];
        return a * d - b * c;
    }

    auto calc_Dodgson_det() {
        while(det_matrix.size() > 0){
            int n = main_matrix.size();
            print_BiMatrix();
            VV tmp_main(n-1, vector<double>(n-1, 0));
            for (int i = 0; i < n-1; i++) {
                for (int j = 0; j < n-1; j++) {
                    tmp_main[i][j] = (calc_2minor_det(i,j) / det_matrix[i][j]);
                }
            }

            VV tmp_det(n-2, vector<double>(n-2, 0));
            for (int i = 1; i < n-1; i++) {
                for (int j = 1; j < n-1; j++) {
                    tmp_det[i-1][j-1] = main_matrix[i][j];
                }
            }
            main_matrix = tmp_main;
            det_matrix = tmp_det;
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

auto calc_Dodgson_det(BiMatrix biMatrix) {
    VV main_matrix = biMatrix.getMainMatrix();
    VV det_matrix = biMatrix.getDetMatrix();
    while(det_matrix.size() > 0){
        int n = main_matrix.size();
        biMatrix.print_BiMatrix();
        VV tmp_main(n-1, vector<double>(n-1, 0));
        for (int i = 0; i < n-1; i++) {
            for (int j = 0; j < n-1; j++) {
                tmp_main[i][j] = (biMatrix.calc_2minor_det(i,j) / det_matrix[i][j]);
            }
        }

        VV tmp_det(n-2, vector<double>(n-2, 0));
        for (int i = 1; i < n-1; i++) {
            for (int j = 1; j < n-1; j++) {
                tmp_det[i-1][j-1] = main_matrix[i][j];
            }
        }
        main_matrix = tmp_main;
        det_matrix = tmp_det;
        biMatrix = BiMatrix(main_matrix, det_matrix);
    }
    return main_matrix[0][0];
}

auto calc_Chio_det(Matrix matrix) {
    while (matrix.getSize() > 1) {
        int nn = matrix.getSize();
        double a = matrix[0][0];
        cout << "\t\033[32m" << a << "    " << "\033[37m" << '\n';
        matrix.print_matrix();
        VV tmp_main(nn - 1, vector<double>(nn - 1, 0));
        for (int i = 0; i < nn - 1; i++) {
            for (int j = 0; j < nn - 1; j++) {
                if (i == 0) {
                    tmp_main[i][j] = matrix.calc_2minor_det(0, i + 1, 0, j + 1) / pow(a, nn-2);
                }
                else {
                    tmp_main[i][j] = matrix.calc_2minor_det(0, i + 1, 0, j + 1);
                }

            }
        }
        matrix = Matrix(tmp_main);
    }
    return matrix[0][0];
}

int main() {
    ifstream in("1.txt");
    int n = 0;
    in >> n;
    VV vv(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            in >> vv[i][j];
        }
    }
    BiMatrix bm(vv);
    // cout << bm.calc_Dodgson_det();
    cout << calc_Dodgson_det(bm);
    cout << '\n';
    cout << "-------------CHIO's METHOD-------------\n";
    Matrix m(vv);
    cout << calc_Chio_det(m) << '\n';
}
