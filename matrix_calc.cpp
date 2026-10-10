#include<iostream>
#include<fstream>
#include<vector>
#include<iomanip>
#include<math.h>
#include<string>

class Matrix;
class BiMatrix;
double calc_Chio_det(Matrix matrix);
double calc_Dodgson_det(Matrix matrix);



using namespace std;



template<typename T>
using VVt = vector<vector<T>>;

using VV = VVt<double>;

template<typename T>
void print_vec(vector<T> v) {
    for (T e: v) cout << e+1 << ' ';
    cout << '\n';
}


class Matrix {
    // 0-index

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

    double calc_2minor_det(int row, int col) {
        auto a = matrix[row][col];
        auto b = matrix[row][col + 1];
        auto c = matrix[row + 1][col];
        auto d = matrix[row + 1][col + 1];
        return a * d - b * c;
    }


    double calc_2minor_det(int r1, int r2, int c1, int c2) {
        auto a = matrix[r1][c1];
        auto b = matrix[r1][c2];
        auto c = matrix[r2][c1];
        auto d = matrix[r2][c2];
        return a * d - b * c;
    }

    double calc_minor_det(vector<int> r, vector<int> c) {
        if (r.size() != c.size()) {
            throw new exception;
        }
        VV minor(r.size(), vector<double>(c.size(), 0));
        for (int i = 0; i < r.size(); i++) {
            for (int j = 0; j < c.size(); j++) {
                minor[i][j] = this->matrix[r[i]][c[j]];
            }
        }
        Matrix m(minor);
        double chio = calc_Chio_det(m);
        return chio;
    }

    vector<double> operator[](int index) const {
        return this->matrix[index];
    }

    void print_matrix(string row_prefix = "    ") {
        for (int i = 0; i < matrix.size(); i++) {
            for (int j = 0; j < matrix.size(); j++) {
                cout << row_prefix << setw(3) << matrix[i][j] << ' ';
            }
            cout << '\n';
        }
        cout << '\n';
    }

};



class BiMatrix {
    // 0-index

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
        origin_matrix(main_matrix), main_matrix(main_matrix), det_matrix(det_matrix) {
            if (this->main_matrix.size() - this->det_matrix.size() != 1) throw new exception;
        }

    BiMatrix(Matrix matrix) {
        this->origin_matrix = matrix.getMatrix();
        this->main_matrix = matrix.getMatrix();
        this->det_matrix = VV(this->main_matrix.size() - 1, vector<double>(this->main_matrix.size() - 1, 1));
    }

    VV getDetMatrix(){
        return this->det_matrix;
    }

    VV getMainMatrix(){
        return this->main_matrix;
    }

    int getSize() {
        return this->main_matrix.size();
    }

    auto calc_2minor_det(int row, int col) {
        auto a = main_matrix[row][col];
        auto b = main_matrix[row][col + 1];
        auto c = main_matrix[row + 1][col];
        auto d = main_matrix[row + 1][col + 1];
        return a * d - b * c;
    }

    auto calc_2minor_det(int r1, int r2, int c1, int c2) {
        auto a = main_matrix[r1][c1];
        auto b = main_matrix[r1][c2];
        auto c = main_matrix[r2][c1];
        auto d = main_matrix[r2][c2];
        return a * d - b * c;
    }


    // all det_matrix can't have zero elements
    // centre element (a[1][1]) can't be zero
    auto calc_3minor_det(vector<int> r, vector<int> c) {
        cout << "-------------------->>>\n";
        cout << "BIMATRIX DET\n";
        if (r.size() != c.size()) {
            throw new exception;
        }
        VV main_minor(r.size(), vector<double>(c.size(), 0));
        for (int i = 0; i < r.size(); i++) {
            for (int j = 0; j < c.size(); j++) {
                main_minor[i][j] = this->main_matrix[r[i]][c[j]];
            }
        }
        VV det_minor(r.size() - 1, vector<double>(c.size() - 1, 0));
        for (int i = 0; i < r.size() - 1; i++) {
            for (int j = 0; j < c.size() - 1; j++) {
                det_minor[i][j] = this->det_matrix[r[i]][c[j]];
            }
        }
        Matrix a(main_minor);
        Matrix b(det_minor);
        BiMatrix bm_minor(main_minor, det_minor);
        cout << "BIMINOR:\n";
        bm_minor.print_BiMatrix();
        if (a.getSize() == 1) {
            cout << "RESULT: " << a[0][0] << "\n";
            cout << "<<<--------------------\n";
            return a[0][0];
        }
        if (a.getSize() == 2) {
            double res = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / (b[0][0]);
            cout << "RESULT: " << res << "\n";
            cout << "<<<--------------------\n";
            return res;
        }
        cout
            << "(\n"
            << "    (\n"
            << "        " << a[0][1] << " * " << a[1][0] << " * "
                            << a[1][2] << " * " << a[2][1] << "\n"
            << "        *\n"
            << "        (" << b[0][0] << " * " << b[1][1]
                            << " - " << b[0][1] << " * " << b[1][0]
                            << ") / (-" << a[1][1] << ")\n"
            << "    )\n"
            << "    -\n"
            << "    (\n"
            << "        + " << a[1][0] << " * (\n"
            << "            " << a[0][1] << " * " << a[2][2]
                            << " * " << b[0][1] << " * " << b[1][0] << "\n"
            << "            - " << a[0][2] << " * " << a[2][1]
                            << " * " << b[0][0] << " * " << b[1][1] << "\n"
            << "        )\n"
            << "        - " << a[1][1] << " * (\n"
            << "            " << a[0][0] << " * " << a[2][2]
                            << " * " << b[0][1] << " * " << b[1][0] << "\n"
            << "            - " << a[0][2] << " * " << a[2][0]
                            << " * " << b[0][0] << " * " << b[1][1] << "\n"
            << "        )\n"
            << "        + " << a[1][2] << " * (\n"
            << "            " << a[0][0] << " * " << a[2][1]
                            << " * " << b[0][1] << " * " << b[1][0] << "\n"
            << "            - " << a[0][1] << " * " << a[2][0]
                            << " * " << b[0][0] << " * " << b[1][1] << "\n"
            << "        )\n"
            << "    )\n"
            << ")\n"
            << "/\n"
            << "(" << b[0][0] << " * " << b[0][1] << " * "
            << b[1][0] << " * " << b[1][1] << ");\n";
        double res =
            (
                 (
                    a[0][1] * a[1][0] * a[1][2] * a[2][1]
                    *
                    (b[0][0] * b[1][1] - b[0][1] * b[1][0]) / (-a[1][1])
                 )
                 -
                 (
                    + a[1][0] * (
                        a[0][1] * a[2][2] * b[0][1] * b[1][0]
                        - a[0][2] * a[2][1] * b[0][0] * b[1][1]
                    )
                    - a[1][1] * (
                        a[0][0] * a[2][2] * b[0][1] * b[1][0]
                        - a[0][2] * a[2][0] * b[0][0] * b[1][1]
                    )
                    + a[1][2] * (
                        a[0][0] * a[2][1] * b[0][1] * b[1][0]
                        - a[0][1] * a[2][0] * b[0][0] * b[1][1]
                    )
                 )
             )
            /
            (b[0][0] * b[0][1] * b[1][0] * b[1][1]);
        cout << "RESULT: " << res << "\n";
        cout << "<<<--------------------\n";
        return res;
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

double calc_Dodgson_det(BiMatrix biMatrix) {
    cout << "-------------------->>>\n";
    cout << "DODJSON METHOD\n";
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
    cout << "RESULT: " << main_matrix[0][0] << "\n";
    cout << "<<<--------------------\n";
    return main_matrix[0][0];
}

double calc_Chio_det(Matrix matrix) {
    cout << "-------------------->>>\n";
    cout << "CHIO METHOD\n";
    while (matrix.getSize() > 1) {
        int nn = matrix.getSize();
        double a = matrix[0][0];
        cout << "\033[32m" << a << "    " << "\033[37m" << '\n';
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
    cout << "RESULT: " << matrix[0][0] << "\n";
    cout << "<<<--------------------\n";
    return matrix[0][0];
}

double calc_common_Chio_det(Matrix matrix, int k) {
    cout << "MATRIX:\n";
    matrix.print_matrix();
    cout << "k = " << k << '\n';
    int n = matrix.getSize();
    vector<int> r(k, -1);
    vector<int> c(k, -1);
    cout << "===1-index===\n";
    cout << "Enter minor " << k << " rows:\n";
    for (int i = 0; i < k; i++) {
        cin >> r[i];
        --r[i];
    }
    cout << "Enter minor " << k << " cols:\n";
    for (int i = 0; i < k; i++) {
        cin >> c[i];
        --c[i];
    }
    auto mdet = matrix.calc_minor_det(r, c);
    int power = n - k - 1;
    // lexicographic order of r and c
    VVt<int> border_minor_rows;
    int previous = -1;
    for (int i = 0; i <= k; i++) {
        int current = i < k ? r[i] : n;
        for (int j = previous + 1; j < current; j++) {
            vector<int> rows;
            rows.insert(rows.end(), r.begin(), r.begin() + i);
            rows.push_back(j);
            rows.insert(rows.end(), r.begin() + i, r.end());
            border_minor_rows.push_back(rows);
        }
        previous = current;
    }
    VVt<int> border_minor_cols;
    previous = -1;
    for (int i = 0; i <= k; i++) {
        int current = i < k ? c[i] : n;
        for (int j = previous + 1; j < current; j++) {
            vector<int> cols;
            cols.insert(cols.end(), c.begin(), c.begin() + i);
            cols.push_back(j);
            cols.insert(cols.end(), c.begin() + i, c.end());
            border_minor_cols.push_back(cols);
        }
        previous = current;
    }
    // calc
    VV less_matrix(n - k, vector<double>(n - k, 0));
    for (int i = 0; i < border_minor_rows.size(); i++) {
        vector<int> rows = border_minor_rows[i];
        for (int j = 0; j < border_minor_cols.size(); j++) {
            vector<int> cols = border_minor_cols[j];
            cout << "\t  ";
            print_vec<int>(cols);
            cout << "\tA\n";
            cout << "\t  ";
            print_vec<int>(rows);
            less_matrix[i][j] = matrix.calc_minor_det(rows, cols);
        }
    }
    Matrix lm(less_matrix);
    cout << "1 / " << "(" << mdet << "^" << power << ") = " << "1 / " << pow(mdet, power) << '\n';
    lm.print_matrix();
    cout << '\n';
}

double calc_common_Dodgson_det(BiMatrix bm, int k) {
    int n = bm.getSize();
    bm.print_BiMatrix();
    Matrix bmain(bm.getMainMatrix());
    Matrix bdet(bm.getDetMatrix());
    VVt<int> rc_main_minor;
    for (int i = 0; i < n - k; i++) {
        vector<int> rc;
        for (int j = i; j < i + k + 1; j++) {
            rc.push_back(j);
        }
        rc_main_minor.push_back(rc);
    }
    VVt<int> rc_det_minor;
    for (int i = 1; i < n - k; i++) {
        vector<int> rc;
        for (int j = i; j < i + k; j++) {
            rc.push_back(j);
        }
        rc_det_minor.push_back(rc);
    }
    VV new_main(n - k, vector<double>(n - k));
    for (int i = 0; i < rc_main_minor.size(); i++) {
        auto rows = rc_main_minor[i];
        for (int j = 0; j < rc_main_minor.size(); j++) {
            auto cols = rc_main_minor[j];
            cout << "\t  ";
            print_vec<int>(cols);
            cout << "\tA\n";
            cout << "\t  ";
            print_vec<int>(rows);
            // bimatrix minor
            new_main[i][j] = bm.calc_3minor_det(rows, cols);
        }
    }
    VV new_det(n - k - 1, vector<double>(n - k - 1));
    cout << "rc_det_minor.size()" << rc_det_minor.size() << '\n';
    for (int i = 0; i < rc_det_minor.size(); i++) {
        auto rows = rc_det_minor[i];
        for (int j = 0; j < rc_det_minor.size(); j++) {
            auto cols = rc_det_minor[j];
            cout << "\t  ";
            print_vec<int>(cols);
            cout << "\tB\n";
            cout << "\t  ";
            print_vec<int>(rows);
            new_det[i][j] = bm.calc_3minor_det(rows, cols);
        }
    }
    BiMatrix new_bm(new_main, new_det);
    new_bm.print_BiMatrix();
}

Matrix read_matrix() {
    ifstream in("1.txt");
    int n = 0;
    in >> n;
    VV vv(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            in >> vv[i][j];
        }
    }
    return Matrix(vv);
}




BiMatrix read_biMatrix() {
    ifstream in("1.txt");
    int n = 0;
    in >> n;
    VV main(n, vector<double>(n, 0));
    VV det(n - 1, vector<double>(n - 1, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            in >> main[i][j];
        }
    }
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - 1; j++) {
            in >> det[i][j];
        }
    }
    return BiMatrix(main, det);
}

int main() {
    Matrix m1 = read_matrix();
    BiMatrix bm1 = read_biMatrix();
    cout << "-------------DODGSON's METHOD-------------\n";
    calc_Dodgson_det(BiMatrix(m1));
    cout << '\n';
//    cout << "-------------THIRD ORDER BIMATRIX DET-------------\n";

    cout << "-------------CHIO's METHOD-------------\n";
    calc_Chio_det(m1);

    cout << "-------------DODGSON'S COMMON METHOD-------------\n";
    calc_common_Dodgson_det(m1, 2);

    cout << "-------------CHIO'S COMMON METHOD-------------\n";
    calc_common_Chio_det(m1, 2);

}
