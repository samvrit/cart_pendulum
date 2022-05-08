#include "matrix_operations.h"

void matrix_vector_multiply(const float A[N_STATES][N_STATES], const float x[N_STATES], float output[N_STATES])
{
    for(int i = 0; i < N_STATES; i++)
	{
	    float sum = 0.0f;
	    for(int j = 0; j < N_STATES; j++)
	    {
	        sum += A[i][j] * x[j];
	    }
	    output[i] = sum;
	}
}

void matrix_matrix_multiply(const float A[N_STATES][N_STATES], const float B[N_STATES][N_STATES], float output[N_STATES][N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			output[i][j] = 0.0f;
			for (int k = 0; k < N_STATES; k++)
			{
				output[i][j] += (A[i][k] * B[k][j]);
			}
		}
	}
}

float dot_product(const float K[N_STATES], const float x[N_STATES])
{
    float sum = 0.0f;
    for(int i = 0; i < N_STATES; i++)
    {
        sum += K[i] * x[i];
    }
    return sum;
}

void vector_sum(const float x1[N_STATES], const float x2[N_STATES], float output[N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		output[i] = x1[i] + x2[i];
	}
}

void matrix_sum(const float A[N_STATES][N_STATES], const float B[N_STATES][N_STATES], float output[N_STATES][N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			output[i][j] = A[i][j] + B[i][j];
		}
	}
}

void matrix_diff(const float A[N_STATES][N_STATES], const float B[N_STATES][N_STATES], float output[N_STATES][N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			output[i][j] = A[i][j] - B[i][j];
		}
	}
}

void vector_diff(const float x1[N_STATES], const float x2[N_STATES], float output[N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		output[i] = x1[i] - x2[i];
	}
}

void vector_scale(const float x[N_STATES], const float a, float output[N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		output[i] = x[i] * a;
	}
}

void matrix_scale(const float A[N_STATES][N_STATES], const float a, float output[N_STATES][N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			output[i][j] = A[i][j] * a;
		}
	}
}

void matrix_transpose(const float A[N_STATES][N_STATES], float output[N_STATES][N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			output[i][j] = A[j][i];
		}
	}
}

void matrix_assign(const float input[N_STATES][N_STATES], float output[N_STATES][N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			output[i][j] = input[i][j];
		}
	}
}

void vector_initialize(float vector[N_STATES], const float init_val)
{
	for (int i = 0; i < N_STATES; i++)
	{
		vector[i] = init_val;
	}
}

void matrix_initialize(float matrix[N_STATES][N_STATES], const float init_val)
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			matrix[i][j] = init_val;
		}
	}
}

void identity(float matrix[N_STATES][N_STATES])
{
	for (int i = 0; i < N_STATES; i++)
	{
		for (int j = 0; j < N_STATES; j++)
		{
			matrix[i][j] = (i == j) ? 1.0f : 0.0f;
		}
	}
}

/*-------------------------- Inverse functions -----------------------------*/

// Adapted from geeks for geeks: https://www.geeksforgeeks.org/adjoint-inverse-matrix/

// Function to get cofactor of A[p][q] in temp[][]. n is current
// dimension of A[][]
static void get_cofactor(float A[N_STATES][N_STATES], float temp[N_STATES][N_STATES], int p, int q, int n)
{
    int i = 0, j = 0;
 
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q)
            {
                temp[i][j++] = A[row][col];
 
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

static float determinant(const float A[N_STATES][N_STATES], int n)
{
    float D = 0.0f;

    if (n == 1)
        return A[0][0];

    float temp[N_STATES][N_STATES];

    float sign = 1.0f;  // To store sign multiplier

    for (int f = 0; f < n; f++)
    {
        get_cofactor(A, temp, 0, f, n);
        D += sign * A[0][f] * determinant(temp, n - 1);

        sign = -sign;
    }

    return D;
}

static void adjoint(const float A[N_STATES][N_STATES], float adj[N_STATES][N_STATES])
{
    if (N_STATES == 1)
    {
        adj[0][0] = 1.0f;
        return;
    }

    float sign = 1.0f, temp[N_STATES][N_STATES] = {{0.0f}};
 
    for (int i = 0; i < N_STATES; i++)
    {
        for (int j = 0; j < N_STATES; j++)
        {
            // Get cofactor of A[i][j]
            get_cofactor(A, temp, i, j, N_STATES);
 
            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ( (i + j) % 2 == 0 )? 1.0f: -1.0f;
 
            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign) * (determinant(temp, N_STATES - 1));
        }
    }
}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool matrix_inverse(const float A[N_STATES][N_STATES], float inverse[N_STATES][N_STATES])
{
    // Find determinant of A[][]
    float det = determinant(A, N_STATES);
    if (det == 0)
    {
        return false;
    }
 
    // Find adjoint
    float adj[N_STATES][N_STATES] = {{0.0f}};
    adjoint(A, adj);
 
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i = 0; i < N_STATES; i++)
	{
        for (int j = 0; j < N_STATES; j++)
		{
            inverse[i][j] = adj[i][j] / det;
		}
	}
 
    return true;
}

