#include <iostream>
#include <fstream>

struct Node {
    double value; // Value of the function at a given point
    double divided_difference; // Divided difference for a node
};

int main() {
    // Open the file with data
    std::ifstream file("MNp1.txt");
    if (!file.is_open()) {
        std::cout << "ERROR: problem with opening the file \"MNp1.txt\"!" << std::endl;
        return 1;
    }

    int d;
    file >> d; // Read the number of interpolation nodes

    // Initialize arrays for interpolation nodes, function values, and divided differences
    double *x = new double[d]; // Array for interpolation nodes
    double *y = new double[d]; // Array for function values at nodes
    Node *divided_differences = new Node[d]; // Array for divided differences

    // Read data from the file
    for (int i = 0; i < d; i++) {
        file >> x[i] >> y[i]; // Read the next interpolation node and its corresponding function value
        divided_differences[i].value = y[i]; // Assign the function value to the 'value' field
        divided_differences[i].divided_difference = y[i]; // Initial divided differences are equal to function values
    }

    file.close(); // Close the file after reading data

    // Calculate divided differences
    for (int i = 1; i < d; i++) { // Iterate over subsequent divided differences
        for (int j = d - 1; j >= i; j--) { // Iterate over nodes from the end to the current index
            // Calculate the next divided difference
            divided_differences[j].divided_difference = (divided_differences[j].divided_difference - divided_differences[j - 1].divided_difference) / (x[j] - x[j - i]);
        }
    }

    double w;
    std::cout << "Enter the point where you want to evaluate the polynomial: ";
    std::cin >> w; // Read the point at which we evaluate the polynomial

    // Calculate the value of the polynomial at point w
    double result = 0.0;
    double term = 1.0; // Terms contributing to the polynomial value

    for (int i = 0; i < d; i++) {
        result += term * divided_differences[i].divided_difference;
        term *= (w - x[i]);
    }

    // Display the results
    std::cout << "Number of nodes: " << d << std::endl; // Display the number of interpolation nodes
    std::cout << "Data:" << std::endl; // Display the data (interpolation nodes and their corresponding function values)
    for (int i = 0; i < d; i++) {
        std::cout << "  x[" << i << "] = " << x[i] << ", y[" << i << "] = " << y[i] << std::endl;
    }
    std::cout << "Your value: " << w << std::endl; // Display the point at which we evaluate the polynomial
    std::cout << "Value of the Newton's polynomial for the chosen value: " << result << std::endl; // Display the value of the polynomial at the given point

    // Display the coefficients of the Newton's polynomial
    std::cout << "Newton's polynomial coefficients bk: ";
    for (int i = 0; i < d; i++) {
        std::cout << divided_differences[i].divided_difference;
        if (i < d - 1) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;

    // Free allocated memory
    delete[] x;
    delete[] y;
    delete[] divided_differences;

    return 0;
}
