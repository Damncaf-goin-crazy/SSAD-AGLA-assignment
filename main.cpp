#include <iostream>
#include <vector>
#include <bits/stdc++.h>

using namespace std;

//multiplication of matrices
vector<vector<double>> multiply_matrices(vector<vector<double>> A, vector<vector<double>> B) {
    int rows_A = A.size();
    int cols_A = A[0].size();
    int cols_B = B[0].size();
    vector<vector<double>> C(rows_A, vector<double>(cols_B, 0.0));
    for (int i = 0; i < rows_A; i++) {
        for (int j = 0; j < cols_B; j++) {
            for (int k = 0; k < cols_A; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

vector<double> multiply_matrix_vector(vector<vector<double>> A, vector<double> B) {
    int rows_A = A.size();
    int cols_A = A[0].size();
    int cols_B = B.size();
    vector<double> C(rows_A);
    for (int i = 0; i < rows_A; i++) {
        for (int k = 0; k < cols_A; k++) {
            C[i] += A[i][k] * B[k];
        }
    }
    return C;
}

//transpose matrix function
vector<vector<double>> transpose_matrix(vector<vector<double>> A) {
    int rows = A.size();
    int cols = A[0].size();
    vector<vector<double>> A_T(cols, vector<double>(rows, 0.0));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            A_T[j][i] = A[i][j];
        }
    }
    return A_T;
}

//inversion of matrices
vector<vector<double>> invert_matrix(vector<vector<double>> A) {
    int n = A.size();
    vector<vector<double>> B(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        B[i][i] = 1.0;
    }
    for (int k = 0; k < n; k++) {
        double max_el = abs(A[k][k]);
        int max_row = k;
        for (int i = k + 1; i < n; i++) {
            if (abs(A[i][k]) > max_el) {
                max_el = abs(A[i][k]);
                max_row = i;
            }
        }

        swap(A[k], A[max_row]);
        swap(B[k], B[max_row]);
        for (int i = k + 1; i < n; i++) {
            double c = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] -= c * A[k][j];
            }
            for (int j = 0; j < n; j++) {
                B[i][j] -= c * B[k][j];
            }
        }
    }
    for (int k = n - 1; k >= 0; k--) {
        for (int i = 0; i < k; i++) {
            double c = A[i][k] / A[k][k];
            for (int j = 0; j < n; j++) {
                B[i][j] -= c * B[k][j];
            }
        }
        double d = 1.0 / A[k][k];
        for (int j = 0; j < n; j++) {
            B[k][j] *= d;
        }
    }
    return B;
}

vector<double> least_square_approximation(vector<pair<double, double>> data, int n) {

    int m = data.size();
    vector<vector<double>> A(m, vector<double>(n + 1, 1.0));
    vector<vector<double>> A_T;
    vector<vector<double>> A_T_A;
    vector<vector<double>> A_T_A_inv;
    vector<double> A_T_b;
    vector<double> b(m, 0.0);
    vector<double> x(n + 1, 0.0);
    for (int i = 0; i < m; i++) {
        for (int j = 1; j < n + 1; j++) {
            A[i][j] = pow(data[i].first, j);
        }
        b[i] = data[i].second;
    }
    A_T = transpose_matrix(A);
    A_T_A = multiply_matrices(A_T, A);
    A_T_A_inv = invert_matrix(A_T_A);
    A_T_b = multiply_matrix_vector(A_T, b);
    x = multiply_matrix_vector(A_T_A_inv, A_T_b);
    cout << "A:" << endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n + 1; j++) {
            cout << fixed << setprecision(4) << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << "A_T*A:" << endl;
    for (int i = 0; i < n + 1; i++) {
        for (int j = 0; j < n + 1; j++) {
            cout << fixed << setprecision(4) << A_T_A[i][j] << " ";
        }
        cout << endl;
    }
    cout << "(A_T*A)^-1:" << endl;
    for (int i = 0; i < n + 1; i++) {
        for (int j = 0; j < n + 1; j++) {
            cout << fixed << setprecision(4) << A_T_A_inv[i][j] << " ";
        }
        cout << endl;
    }
    cout << "A_T*b:" << endl;
    for (int i = 0; i < n + 1; i++) {
        cout << fixed << setprecision(4) << A_T_b[i] << endl;
    }
      cout << "x~:" << endl;
     for (int i = 0; i < n + 1; i++)
     {
         cout << fixed << setprecision(4) << x[i] << endl;
     }
    return x;
}

int main() {
    int m, n;
    vector<pair<double, double>> data;
    cin >> m;
    for (int i = 0; i < m; i++) {
        double t, b;
        cin >> t >> b;
        data.emplace_back(t, b);
    }
    cin >> n;
    least_square_approximation(data, n);
    return 0;
}
