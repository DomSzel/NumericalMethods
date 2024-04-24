#include <iostream>
#include <fstream>
#include <cmath>

// Maximum size for the matrix
const int max_size = 10;

// Function to print the matrix
void printMatrix(double matrix[max_size][max_size + 1], int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= n; ++j) {
            std::cout << matrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

// Function to perform forward elimination
void eliminateForward(double matrix[max_size][max_size + 1], int n) {
    for (int i = 0; i < n; ++i) {
        if (matrix[i][i] == 0.0) {
            std::cout << "Zero on the diagonal! Cannot continue calculations." << std::endl;
            return;
        }
        for (int k = i + 1; k < n; ++k) {
            double factor = matrix[k][i] / matrix[i][i];
            for (int j = i; j <= n; ++j) {
                matrix[k][j] -= factor * matrix[i][j];
            }
        }
    }
}

// Function to solve the system of equations backwards
void solveBackwards(double matrix[max_size][max_size + 1], int n) {
    double solutions[max_size];
    for (int i = n - 1; i >= 0; --i) {
        solutions[i] = matrix[i][n];
        for (int j = i + 1; j < n; ++j) {
            solutions[i] -= matrix[i][j] * solutions[j];
        }
        solutions[i] /= matrix[i][i];
    }
    std::cout << "Solution of the system of equations:" << std::endl;
    for (int i = 0; i < n; ++i) {
        std::cout << "x" << i << " = " << solutions[i] << std::endl;
    }
}

int main() {
    // Open the file
    std::fstream file("RURL_dane2.txt");
    // Check if the file is opened successfully
    if (!file) {
        std::cout << "There was an error while opening the file" << std::endl;
        return 1;
    }

    // Read the size of the matrix from the file
    int n;
    file >> n;

    // Create a matrix to store the augmented matrix
    double matrix[max_size][max_size + 1];

    // Read the matrix from the file
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= n; ++j) {
            file >> matrix[i][j];
        }
    }

    file.close();

    // Display the augmented matrix before calculations
    std::cout << "Augmented matrix (before calculations):" << std::endl;
    printMatrix(matrix, n);

    // Perform forward elimination
    eliminateForward(matrix, n);

    // Display the augmented matrix after forward elimination
    std::cout << "Augmented matrix (after forward elimination):" << std::endl;
    printMatrix(matrix, n);

    // Solve the system of equations backwards
    solveBackwards(matrix, n);

    return 0;
}
