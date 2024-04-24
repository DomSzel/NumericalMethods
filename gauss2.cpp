#include <iostream>
#include <fstream>

using namespace std;

const int MAX_SIZE = 10;

// Function to print the matrix
void printMatrix(double matrix[MAX_SIZE][MAX_SIZE + 1], int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= n; ++j) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

// Function to perform partial pivot Gauss-Crout elimination
void partialPivotGaussCrout(double matrix[MAX_SIZE][MAX_SIZE + 1], int n, int columnOrder[]) {
    for (int i = 0; i < n; ++i) {
        columnOrder[i] = i;
        int maxCol = i;
        for (int j = i + 1; j < n; ++j) {
            if (abs(matrix[j][i]) > abs(matrix[j][maxCol])) {
                maxCol = j;
            }
        }
        swap(columnOrder[i], columnOrder[maxCol]);
        for (int j = 0; j <= n; ++j) {
            swap(matrix[i][j], matrix[maxCol][j]);
        }
        for (int j = i + 1; j < n; ++j) {
            double ratio = matrix[j][i] / matrix[i][i];
            for (int k = i; k <= n; ++k) {
                matrix[j][k] -= ratio * matrix[i][k];
            }
        }
    }
}

// Function to perform partial pivot Gauss elimination
void partialPivotGauss(double matrix[MAX_SIZE][MAX_SIZE + 1], int n) {
    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int j = i + 1; j < n; ++j) {
            if (abs(matrix[j][i]) > abs(matrix[maxRow][i])) {
                maxRow = j;
            }
        }
        for (int k = i; k <= n; ++k) {
            swap(matrix[i][k], matrix[maxRow][k]);
        }
        for (int j = i + 1; j < n; ++j) {
            double ratio = matrix[j][i] / matrix[i][i];
            for (int k = i; k <= n; ++k) {
                matrix[j][k] -= ratio * matrix[i][k];
            }
        }
    }
}

// Function to solve equations using the resulting upper triangular matrix
void solveEquations(double matrix[MAX_SIZE][MAX_SIZE + 1], int n, double solutions[MAX_SIZE]) {
    for (int i = n - 1; i >= 0; --i) {
        solutions[i] = matrix[i][n];
        for (int j = i + 1; j < n; ++j) {
            solutions[i] -= matrix[i][j] * solutions[j];
        }
        solutions[i] /= matrix[i][i];
    }
}

int main() {
    ifstream inputFile("RURL_data2.txt");
    if (!inputFile) {
        cerr << "Error opening file!";
        return 1;
    }

    int n;
    inputFile >> n;

    double matrix[MAX_SIZE][MAX_SIZE + 1];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= n; ++j) {
            inputFile >> matrix[i][j];
        }
    }
    inputFile.close();

    cout << "Matrix before computations:\n";
    printMatrix(matrix, n);
    cout << endl;

    int columnOrder[MAX_SIZE];

    partialPivotGaussCrout(matrix, n, columnOrder);

    cout << "Matrix after the first step of computations - forward process:\n";
    printMatrix(matrix, n);
    cout << endl;

    cout << "Column order after swaps:\n";
    for (int i = 0; i < n; ++i) {
        cout << columnOrder[i] << " ";
    }
    cout << endl << endl;

    double solutionsGaussCrout[MAX_SIZE];
    solveEquations(matrix, n, solutionsGaussCrout);

    cout << "Solution of equations using Gauss-Crout method:\n";
    for (int i = 0; i < n; ++i) {
        cout << "x" << i << " = " << solutionsGaussCrout[i] << endl;
    }
    cout << endl;

    inputFile.open("RURL_data2.txt");
    if (!inputFile) {
        cerr << "Error opening file!";
        return 1;
    }

    inputFile >> n;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= n; ++j) {
            inputFile >> matrix[i][j];
        }
    }
    inputFile.close();

    cout << "Matrix before computations using partial pivot Gauss method:\n";
    printMatrix(matrix, n);
    cout << endl;

    partialPivotGauss(matrix, n);

    cout << "Matrix after the first step of computations - forward process using partial pivot Gauss method:\n";
    printMatrix(matrix, n);
    cout << endl;

    double solutionsPartialPivotGauss[MAX_SIZE];
    solveEquations(matrix, n, solutionsPartialPivotGauss);

    cout << "Solution of equations using partial pivot Gauss method:\n";
    for (int i = 0; i < n; ++i) {
        cout << "x" << i << " = " << solutionsPartialPivotGauss[i] << endl;
    }

    return 0;
}
