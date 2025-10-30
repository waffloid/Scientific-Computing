#include <iostream>
#include <vector>

void calcNextRow(const int* prevRow, int* nextRow, int rowNo) {
	nextRow[0] = 1;
	nextRow[rowNo] = 1;

	for (int j = 1; j<rowNo; j++) {
		nextRow[j] = prevRow[j-1] + prevRow[j];
	}
	
	for (int k = 0; k < rowNo; k++) {
		std::cout << prevRow[k] << " ";
	}
	std::cout << std::endl;
}

void calculateRows(int noRows) {
	int* prevRow = nullptr;
	int* currRow = new int[1];

	currRow[0] = 1;

	for (int i = 1; i<noRows; i++) {
		if (prevRow) {
			delete[] prevRow;
		}

		prevRow = currRow;
		currRow = new int[i+1];

		calcNextRow(prevRow, currRow, i);
	}
}

/*
Write separate functions that do the following:
1. Allocate space for an N by N matrix and return a double pointer that can be used to access it. 2. Delete the memory allocated by the first function.
1
The functions should be usable in the following way:
  int** matrix = allocateMatrix(3);
  int m22 = matrix[2][2];
  freeMatrix(matrix);

 Hint: You will need the construct int** a = new int*[N] to allocate an array of pointers to each row. Write a function that finds the product of two matrices declared in this way. Allow a user to input a size
N and then two N × N matrices and use this function to find their product.
Note that there are (at least) two ways to allocate the matrix memory, either as a single contiguous block,
or as separate rows. For practice, you should try both approaches.
 */



int** allocateMatrix(int N) {
	int** matrix = nullptr;
	matrix[0] = new int[N * N];

	for (int i = 0; i < N; ++i) {
		matrix[i] = *matrix + i * N;
	}

	return matrix;
}

void deleteMatrix(int** matrix) {
	delete matrix;
}

int** multiplyMatrices(int** A, int** B, int N) {
	int** output = allocateMatrix(N);

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			for (int k = 0; k < N; k++) {
				output[i][j] = A[i][k] * B[k][j];
			}
		}
	}

	return output;
}

int main() {
	int noRows;

	std::cin >> noRows;

	calculateRows(noRows);



	return 0;
}



