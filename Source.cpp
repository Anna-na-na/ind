#include <stdio.h>
#include <iostream>
#include "mpi.h"
#include <time.h> 

using namespace std;

/// <summary>
/// ���������������� �������
/// </summary>
/// <param name="A">�������</param>
/// <param name="n">����</param>
void Flip(double*& A, int n)
{
	double* temp;
	temp = new double[n * n];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			*(temp + i * n + j) = *(A + i * n + j);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			A[i * n + j] = temp[j * n + i];
}

/// <summary>
/// ���������������� ������������ ������
/// </summary>
/// <param name="A">������ ���������</param>
/// <param name="B">������ ���������</param>
/// <param name="C">������������</param>
/// <param name="n">����</param>
void MatrixMultiplication(double*& A, double*& B, double*& C, int& n) {
	int i, j, k;
	Flip(B, n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			C[i * n + j] = 0;
			for (k = 0; k < n; k++)
				C[i * n + j] += A[i * n + k] * B[j * n + k];
		}
	}
}

int main(int* argc, char** argv)
{
	int size, rank;

	MPI_Init(argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// ������������� ������
	double* A, * B, * C;
	int n = 1000;
	A = new double[n * n];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			*(A + i * n + j) = rand() % 10;
		}
	B = new double[n * n];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			*(B + i * n + j) = rand() % 10;
		}
	C = new double[n * n];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			*(C + i * n + j) = 0;
		}

	/*cout << "Matrix A:\n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) {
			cout << A[i * n + j] << " ";
		}
		cout << endl;
	}
	cout << "\nMatrix B:\n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) {
			cout << B[i * n + j] << " ";
		}
		cout << endl;
	}
	cout << "\n\n";*/

	// ���������������� ������������
	if (rank == 0) {
		double start = MPI_Wtime();
		MatrixMultiplication(A, B, C, n);
		double end = MPI_Wtime();
		double seconds = (end - start) ;
		printf("Sequential time: %f seconds\n\n", seconds);
		Flip(B, n);

		/*cout << "Matrix C:\n";
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++) {
				cout << C[i * n + j] << " ";
			}
			cout << endl;
		}
		cout << "\n\n";*/
	}

	// ������������ ������������
	MPI_Status Status;
	double temp;

	int partedSize = n / size;
	int newsize = partedSize * n;
	double* tempA = new double[newsize];
	double* tempB = new double[newsize];
	double* tempC = new double[newsize];
	if (rank == 0) {
		Flip(B, n);
		// ������� ����� �������
		cout << "size = " << size << endl;
	}
	
	if (size > 1) {

	}
	double start = MPI_Wtime();
	MPI_Scatter(A, newsize, MPI_DOUBLE, tempA, newsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(B, newsize, MPI_DOUBLE, tempB, newsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int i = 0; i < partedSize; i++) {
		for (int j = 0; j < partedSize; j++) {
			temp = 0;
			for (int k = 0; k < n; k++) {
				temp += tempA[i * n + k] * tempB[j * n + k];
			}				
			tempC[i * n + j + partedSize * rank] = temp;
		}
	}

	int r, next, prev;
	for (int p = 1; p < size; p++) {
		// ���� ����� ������� ������, ��� 1, ������� �� ������
		if (p == 1) cout << rank << " proc" << endl;
		next = rank + 1;
		if (rank == size - 1)
			next = 0;
		prev = rank - 1;
		if (rank == 0)
			prev = size - 1;
		MPI_Sendrecv_replace(tempB, newsize, MPI_DOUBLE, next, 0, prev, 0, MPI_COMM_WORLD, &Status);
		for (int i = 0; i < partedSize; i++) {
			for (int j = 0; j < partedSize; j++) {
				temp = 0;
				for (int k = 0; k < n; k++) {
					temp += tempA[i * n + k] * tempB[j * n + k];
				}
				if (rank - p >= 0)
					r = rank - p;
				else r = (rank - p + size);
				tempC[i * n + j + partedSize * r] = temp;
			}
		}
	}

	MPI_Gather(tempC, newsize, MPI_DOUBLE, C, newsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	double end = MPI_Wtime();
	double seconds = (end - start) ;
	if (rank == 0)
	{
		printf("Parallel time: %f seconds\n", seconds);

		/*cout << "Parallel matrix C:\n";
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++) {
				cout << C[i * n + j] << " ";
			}
			cout << endl;
		}
		cout << "\n\n";*/
	}

	delete[]tempA;
	delete[]tempB;
	delete[]tempC;

	MPI_Finalize();
}
