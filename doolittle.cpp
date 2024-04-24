#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip> // std::setprecision

void printMatrix(const std::vector<std::vector<double>>& matrix) {
    // Iterate over each row of the matrix
    for (const auto& row : matrix) {
        // For each element in the row
        for (double elem : row) {
            // Set the width of the displayed field to 8 characters with 4 decimal places
            std::cout << std::setw(8) << std::setprecision(4) << elem << " ";
        }
        std::cout << std::endl;
    }
}

// Decompose the matrix into lower triangular matrix L and upper triangular matrix U
// Return true if successful, false otherwise
bool luDecomposition(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& L,
                     std::vector<std::vector<double>>& U, int n) {
    // Iterate over all rows of the matrix
    for (int i = 0; i < n; ++i) {
        // Compute elements of matrix U
        for (int j = i; j < n; ++j) {
            double sum = 0.0;
            // Sum the products of calculated elements of matrix L and U in the current row
            for (int k = 0; k < i; ++k) {
                sum += L[i][k] * U[k][j];
            }
            // Compute element of matrix U using formula (12)
            U[i][j] = A[i][j] - sum;
        }

        // Compute elements of matrix L
        for (int j = i; j < n; ++j) {
            // If we are on the diagonal, set the element to 1
            if (i == j)
                L[i][i] = 1;
            else {
                double sum = 0.0;
                // Sum the products of calculated elements of matrix L and U in the current row
                for (int k = 0; k < i; ++k) {
                    sum += L[j][k] * U[k][i];
                }
                // If the element of matrix U on the diagonal is zero, return false
                if (U[i][i] == 0.0) {
                    std::cout << "Division by zero!" << std::endl;
                    return false;
                }
                // Compute element of matrix L using formula (13)
                L[j][i] = (A[j][i] - sum) / U[i][i];
            }
        }
    }
    return true;
}

// Solve the LY = B equation system
void solveLY(const std::vector<std::vector<double>>& L, std::vector<double>& Y, const std::vector<double>& B, int n) {
    // Iterate over all rows
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        // Sum the products of elements of matrix L and vector Y
        for (int j = 0; j < i; ++j) {
            sum += L[i][j] * Y[j];
        }
        // Compute element Y[i] using formula (8)
        Y[i] = B[i] - sum;
    }
}

// Solve the UX = Y equation system
void solveUX(const std::vector<std::vector<double>>& U, std::vector<double>& X, const std::vector<double>& Y, int n) {
    // Iterate from the last row to the first
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        // Sum the products of elements of matrix U and vector X
        for (int j = i + 1; j < n; ++j) {
            sum += U[i][j] * X[j];
        }
        // If the element of matrix U on the diagonal is zero, return error
        if (U[i][i] == 0.0) {
            std::cout << "Division by zero!" << std::endl;
            return;
        }
        // Compute element X[i] using formula (11)
        X[i] = (Y[i] - sum) / U[i][i];
    }
}

int main() {
    std::ifstream file("RURL_data2.txt");
    if (!file.is_open()) {
        std::cout << "Cannot open the file." << std::endl;
        return 1;
    }

    int n;
    file >> n;

    std::vector<std::vector<double>> A(n, std::vector<double>(n)), L(n, std::vector<double>(n, 0.0)),
            U(n, std::vector<double>(n, 0.0));
    std::vector<double> B(n), Y(n), X(n);

    // Read matrix A and vector B from the file
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file >> A[i][j];
        }
        file >> B[i];
    }
    file.close();

    // Print the read matrix A
    std::cout << "Matrix A:" << std::endl;
    printMatrix(A);
    std::cout << "Vector B:" << std::endl;
    for (double b : B) {
        std::cout << std::setw(8) << std::setprecision(4) << b << std::endl;
    }

    // Decompose matrix A into matrices L and U
    if (!luDecomposition(A, L, U, n)) {
        return 1;
    }

    // Print matrices L and U
    std::cout << "Matrix L:" << std::endl;
    printMatrix(L);
    std::cout << "Matrix U:" << std::endl;
    printMatrix(U);

    // Solve the LY = B equation system
    solveLY(L, Y, B, n);
    std::cout << "Vector Y:" << std::endl;
    for (double y : Y) {
        std::cout << std::setw(8) << std::setprecision(4) << y << std::endl;
    }

    // Solve the UX = Y equation system
    solveUX(U, X, Y, n);
    std::cout << "Solution X:" << std::endl;
    for (double x : X) {
        std::cout << std::setw(8) << std::setprecision(4) << x << std::endl;
    }

    return 0;
}
