
#include <iostream>
#include <vector>

class Matrix {
private:
    int rows;
    int cols;
    std::vector<std::vector<double>>data;
public:
    Matrix(int rows, int cols):rows(rows),cols(cols),data(rows,std::vector<double>(cols,0.0)){}
    void set(int row, int col, double value) {
        data[row][col] = value;
    }
    double get(int row, int col) const {
        return data[row][col];
    }
    Matrix operator+(const Matrix& other)const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.set(i, j, data[i][j] + other.get(i, j));
            }
        }
        return result;
    }
    Matrix operator-(const Matrix& other)const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.set(i, j, data[i][j] - other.get(i, j));
            }
        }
        return result;
    }
    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            // Matrix dimensions are not compatible for multiplication
            
            throw std::runtime_error("Matrix dimensions are not compatible for multiplication");
        }

        Matrix result(rows, other.cols); // Resulting matrix with appropriate dimensions

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < other.cols; j++) {
                double sum = 0.0;
                for (int k = 0; k < cols; k++) {
                    sum += data[i][k] * other.data[k][j];
                }
                result.set(i, j, sum);
            }
        }

        return result;
    }
    double determinant() const {
        if (rows != cols) {
            // Matrix must be square for determinant calculation
            throw std::runtime_error("Matrix must be square for determinant calculation");
        }

        if (rows == 1) {
            // Base case for 1x1 matrix
            return data[0][0];
        }

        double det = 0.0;
        for (int col = 0; col < cols; col++) {
            double submatrix_det = getSubmatrix(0, col).determinant(); // Recursive call for submatrix determinant
            double cofactor = (col % 2 == 0) ? 1.0 : -1.0; // Alternating sign for cofactor
            det += cofactor * data[0][col] * submatrix_det;
        }

        return det;
    }
    Matrix getSubmatrix(int row, int col) const {
        Matrix submatrix(rows - 1, cols - 1); // Create a submatrix with one row and column less than the current matrix
        int submatrix_row = 0; // Row index for the submatrix
        for (int i = 0; i < rows; i++) {
            if (i == row) {
                // Skip the specified row
                continue;
            }
            int submatrix_col = 0; // Column index for the submatrix
            for (int j = 0; j < cols; j++) {
                if (j == col) {
                    // Skip the specified column
                    continue;
                }
                // Copy the element to the submatrix
                submatrix.set(submatrix_row, submatrix_col, data[i][j]);
                submatrix_col++;
            }
            submatrix_row++;
        }
        return submatrix;
    }
    void print() const {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                std::cout << data[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
};
int main()
{
    Matrix matrix1(2, 2);
    matrix1.set(0, 0, 1.0);
    matrix1.set(0, 1, 2.0);
    matrix1.set(1, 0, 3.0);
    matrix1.set(1, 1, 4.0);

    Matrix matrix2(2, 2);
    matrix2.set(0, 0, 5.0);
    matrix2.set(0, 1, 6.0);
    matrix2.set(1, 0, 7.0);
    matrix2.set(1, 1, 8.0);

    // Perform matrix addition
    Matrix matrixSum = matrix1 - matrix2;

    // Print the result
    std::cout << "Matrix 1:" << std::endl;
    matrix1.print();
    std::cout << "Matrix 2:" << std::endl;
    matrix2.print();
    std::cout << "Matrix Mult:" << std::endl;
    matrixSum.print();
  
    // Calculate determinant
    double det = matrix1.determinant();

    // Print the result
    std::cout << "Matrix 1:" << std::endl;
    matrix1.print();
    std::cout << "Determinant: " << det << std::endl;
    return 0;
}
