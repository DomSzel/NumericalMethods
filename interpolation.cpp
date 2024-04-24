#include <iostream>
#include <fstream>

int main() {
    // Open the file
    std::fstream file("MN.txt");
    // Check if the file is opened successfully
    if (!file.is_open()) {
        // Display error message if file cannot be opened
        std::cout << "ERROR: problem with opening the file \"MN.txt\"!" << std::endl;
        return 1;
    }
    int d;
    // Read the number of nodes from the file
    file >> d;
    // Create dynamic arrays to store node values
    double *x = new double[d]; // x values of nodes
    double *y = new double[d]; // y values of nodes

    // Read node values from the file
    for (int i = 0; i < d; i++) {
        file >> x[i] >> y[i];
    }

    // Close the file
    file.close();

    // Variable to store the point of evaluation
    double w;
    // Ask the user to input the point of evaluation
    std::cout << "Enter the point where you want to evaluate the polynomial: ";
    std::cin >> w;

    // Variable to store the result of polynomial evaluation
    double result = 0;
    // Calculate the polynomial value using Lagrange interpolation
    for (int i = 0; i < d; i++) {
        // Initialize wi to 1 for Lagrange interpolation
        double wi = 1;
        // Calculate the Lagrange polynomial for the given point
        for (int j = 0; j < d; j++) {
            if (i != j) {
                wi *= (w - x[j]) / (x[i] - x[j]);
            }
        }
        // Accumulate the polynomial value for each node
        result += wi * y[i];
    }

    // Display the number of nodes
    std::cout << "Number of nodes: " << d << std::endl;
    // Display node values
    std::cout << "Data:" << std::endl;
    for (int i = 0; i < d; i++) {
        std::cout << "  x[" << i << "] = " << x[i] << ", y[" << i << "] = " << y[i] << std::endl;
    }
    // Display the point of evaluation
    std::cout << "Point of evaluation: " << w << std::endl;
    // Display the value of the polynomial for the chosen point
    std::cout << "Value of the polynomial for the chosen value: " << result << std::endl;

    // Deallocate memory for dynamic arrays
    delete[] x;
    delete[] y;

    return 0;
}
