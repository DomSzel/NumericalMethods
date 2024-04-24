#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

void printMatrix(const std::vector<std::vector<double>>& matrix, const std::string& name) {
    std::cout << name << ":" << std::endl;
    // Print each element of the matrix
    for (const auto& row : matrix) {
        for (double elem : row) {
            // Set width and precision for better formatting
            std::cout << std::setw(10) << std::setprecision(4) << elem << " ";
        }
        std::cout << std::endl;
    }
}

void printVector(const std::vector<double>& vector, const std::string& name) {
    // Print the name of the vector
    std::cout << name << ":" << std::endl;
    // Print each element of the vector
    for (double elem : vector) {
        // Set precision for better formatting
        std::cout << std::setw(10) << std::setprecision(4) << elem << std::endl;
    }
}

bool isDiagonallyDominant(const std::vector<std::vector<double>>& A) {
    int n = A.size();
    bool at_least_one_strict = false;

    // Iterate over each row of the matrix
    for (int i = 0; i < n; i++) {
        double sum = 0;
        // Iterate over each element of the row
        for (int j = 0; j < n; j++) {
            if (i != j) {
                // Add the absolute value of each off-diagonal element
                sum += std::abs(A[i][j]);
            }
        }
        // Check if the diagonal element is strictly greater than the sum of off-diagonal elements
        if (std::abs(A[i][i]) < sum) return false;
        // Check if the diagonal element is strictly greater than the sum of off-diagonal elements for at least one row
        if (std::abs(A[i][i]) > sum) at_least_one_strict = true;
    }
    return at_least_one_strict;
}

std::vector<std::vector<double>> decomposeToLDU(const std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& L, std::vector<std::vector<double>>& U, std::vector<double>& D_inv) {
    int n = A.size();
    // Initialize matrices L and U and vector D_inv with zeros
    L = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
    U = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
    D_inv = std::vector<double>(n, 0.0);

    // Decompose matrix A into matrices L, U, and D_inv
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i > j) {
                L[i][j] = A[i][j]; // Lower triangular matrix
            }
            else if (i < j) {
                U[i][j] = A[i][j]; // Upper triangular matrix
            }
        }
        D_inv[i] = 1 / A[i][i]; // Inverse of diagonal elements
    }
    return L; // Return matrix L
}

void jacobiMethod(const std::vector<std::vector<double>>& A, const std::vector<double>& b, std::vector<double>& x, int max_iter) {
    int n = A.size();
    std::vector<std::vector<double>> L, U;
    std::vector<double> D_inv;
    // Decompose matrix A into matrices L, U, and D_inv
    std::vector<std::vector<double>> L_plus_U = decomposeToLDU(A, L, U, D_inv);

    std::vector<double> temp_x(n, 0.0);
    // Perform Jacobi iteration
    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            // Compute the sum of products of elements in L, U and current solution vector x
            for (int j = 0; j < n; j++) {
                sum += (L[i][j] + U[i][j]) * x[j];
            }
            // Compute the new value of the solution vector x for the current equation
            temp_x[i] = -D_inv[i] * sum + D_inv[i] * b[i];
        }
        // Update the solution vector x for the next iteration
        x = temp_x;
    }

    // Print matrices L+U, vector D^-1, and solution vector x
    printMatrix(L_plus_U, "Matrix L+U");
    printVector(D_inv, "Matrix D^-1");
    printVector(x, "Solution");
}

int main() {
    // Open the file containing input data
    std::ifstream file("RURL_Jacobi_data1.txt");
    if (!file.is_open()) {
        std::cout << "Cannot open the file." << std::endl;
        return 1;
    }

    int n;
    file >> n; // Read the size of matrix A

    // Initialize matrix A, vector b, and solution vector x
    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    std::vector<double> b(n), x(n, 0.0);

    // Read matrix A and vector b from the file
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file >> A[i][j];
        }
        file >> b[i];
    }
    file.close(); // Close the file after reading data

    // Print matrix A and vector b
    printMatrix(A, "Matrix A");
    printVector(b, "Vector b");

    // Check if matrix A is diagonally dominant
    if (!isDiagonallyDominant(A)) {
        std::cout << "The matrix is not diagonally dominant." << std::endl;
        return 1;
    }
    else {
        std::cout << "The matrix is diagonally dominant." << std::endl;
    }

    int max_iter;
    // Prompt the user to enter the maximum number of iterations
    std::cout << "Enter the number of iterations: ";
    std::cin >> max_iter;

    // Apply Jacobi method to solve the system of linear equations
    jacobiMethod(A, b, x, max_iter);

    return 0;
}
