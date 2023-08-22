#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// ran1 parameters
#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1 + (IM - 1) / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)
#define SEED -2023

// Function declarations
int **read_ldpc_G_from_file(const char *, int, int);
int **read_node_connection_matrix_from_file(const char *, int, int);
int **createGeneratorMatrix(int **, int, int);
void freeMatrix(int **, int);
void printMatrix(int **, int, int);
//=========================================
int *generated_sequence(int);
int **create_dynamic_matrix_int(int, int);
double **create_dynamic_matrix_double(int, int);
int *Encoder(int **, int *, int, int);
double *Channel(double *, long *, double, int);
double *Modulator(int *, int);
//==========================================
int *decoder(double **, double **, int **, double *, double, int, int, int, int);
double *CHK(double **, int **, int, int);
double VAR(double **, int **, int, int, int);
double SPA(double, double);
//==========================================
double noise_variance(double);
long *rand_init(long *, long);
double ran1(long *);
void normal(long *, double *, double);

int main()
{
    double SNR = 3.0;
    int block_num = 10000000, iteration_max = 100;

    // Start the timer to measure the execution time of the simulation
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    // Define block length, dimension, and row weight
    int block_length = 1023, dimension = 781, row_weight = 32;

    // Read LDPC matrix from files
    int **ldpc_G = read_ldpc_G_from_file("ldpc_G_1023.txt", dimension, block_length);
    int **node_connection_matrix = read_node_connection_matrix_from_file("ldpc_H_1023.txt", 2 * block_length, row_weight);

    // Initialize the random number generator
    long *idum = (long *)calloc(1, sizeof(long));
    rand_init(idum, SEED);

    // Compute the standard deviation of the noise
    double noise_variance_val = noise_variance(SNR);
    double sigma = sqrt(noise_variance_val);

    // Create dynamic matrices for edges
    double **check_to_variable_msgs = create_dynamic_matrix_double(block_length, block_length);
    double **variable_to_check_msgs = create_dynamic_matrix_double(block_length, block_length);

    // Generate a sequence of data to be transmitted
    int *u_tx = generated_sequence(dimension * 63);

    // Variables for error counting
    int err_num = 0, err_block_num = 0;

    // Loop through each block
    for (int i = 0; i < block_num; i++)
    {
        // Terminate if too many error blocks
        if (err_block_num > 100)
        {
            break;
        }

        // Extract a segment from the transmitted data
        int *u = (int *)calloc(dimension, sizeof(int));
        for (int j = 0; j < dimension; j++)
        {
            u[j] = u_tx[(i % 63) * dimension + j];
        }

        // Encode, Modulate, and Transmit through Channel
        int *c = Encoder(ldpc_G, u, dimension, block_length);
        double *x = Modulator(c, block_length);
        double *y = Channel(x, idum, sigma, block_length);

        // Decode the received signal
        int *u_estimated = decoder(check_to_variable_msgs, variable_to_check_msgs, node_connection_matrix,
                                   y, noise_variance_val, block_length, dimension, row_weight, iteration_max);

        // Check for errors in the decoded block
        int flag = 0;
        for (int j = 0; j < dimension; j++)
        {
            if (u[j] != u_estimated[j])
            {
                flag = 1;
                err_num++;
            }
        }
        if (flag == 1)
        {
            err_block_num++;
        }

        // Print the current status
        printf("block no.%d,  err_block_num = %d,  err_num = %d \n", i + 1, err_block_num, err_num);

        // Free allocated memory
        free(u);
        free(c);
        free(x);
        free(y);
    }
    // Free allocated memory
    free(u_tx);
    freeMatrix(ldpc_G, dimension);
    freeMatrix(node_connection_matrix, 2 * block_length);

    // Calculate and print the total execution time
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("running time: %f sec\n", cpu_time_used);
    system("pause");
    return 0;
}

// Function to read the Generator Matrix from a file and store it in a matrix
int **read_ldpc_G_from_file(const char *filename, int rows, int cols)
{
    // Open the file for reading
    FILE *file = fopen(filename, "r");

    // If the file could not be opened, print an error message and return NULL
    if (file == NULL)
    {
        printf("Unable to open the file.\n");
        return NULL;
    }
    // Create a dynamic 2D array (matrix) to store the data
    int **matrix = create_dynamic_matrix_int(rows, cols);

    char ch;

    // Loop through each row
    for (int i = 0; i < rows; i++)
    {
        // Loop through each column
        for (int j = 0; j < cols; j++)
        {
            // Read a character from the file
            if ((ch = fgetc(file)) != EOF)
            {
                if (ch == '\n')
                {
                    // If the character is a newline, skip this iteration
                    j--;
                    continue;
                }

                // Convert the character to an integer and store it in the matrix
                matrix[i][j] = ch - '0';
            }
            else
            {
                // If the end of file is reached before the matrix is filled, print an error message and return NULL
                printf("Not enough data in file\n");
                return NULL;
            }
        }
    }

    fclose(file);

    return matrix;
}

// Function to read data from file and store in a matrix
int **read_node_connection_matrix_from_file(const char *filename, int rows, int cols)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        printf("Unable to open the file.\n");
        return NULL;
    }

    // Create a dynamic 2D array (matrix) to store the data
    int **matrix = create_dynamic_matrix_int(rows, cols);

    int temp;

    // Loop through each row
    for (int i = 0; i < rows; i++)
    {
        // Loop through each column
        for (int j = 0; j < cols; j++)
        {
            // Read an integer from the file and store it in 'temp'
            if (!fscanf(file, "%d", &temp))
            {
                // If the end of file is reached before the matrix is filled, print an error message and return NULL
                printf("Not enough data in file\n");
                return NULL;
            }
            // Subtract 1 from the read value and store it in the matrix
            matrix[i][j] = temp - 1;
        }
    }

    fclose(file);

    return matrix;
}

// Function to free allocated memory
void freeMatrix(int **matrix, int rows)
{
    for (int i = 0; i < rows; i++)
        free(matrix[i]);
    free(matrix);
}

// Function to print the contents of the matrix
void printMatrix(int **matrix, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

// Function to create a dynamically allocated 2D array (matrix) of integers
int **create_dynamic_matrix_int(int rows, int cols)
{
    int **matrix = (int **)malloc(rows * sizeof(int *));
    for (int i = 0; i < rows; i++)
        matrix[i] = (int *)calloc(cols, sizeof(int));

    return matrix;
}

// Function to create a dynamically allocated 2D array (matrix) of doubles
double **create_dynamic_matrix_double(int rows, int cols)
{
    double **matrix = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++)
        matrix[i] = (double *)calloc(cols, sizeof(double));

    return matrix;
}

// Function to encode a binary input sequence u using a generator matrix ldpc_G
int *Encoder(int **ldpc_G, int *u, int rows, int cols)
{

    int *c = (int *)calloc(cols, sizeof(int));
    for (int i = 0; i < cols; i++)
    {
        for (int j = 0; j < rows; j++)
        {
            c[i] += u[j] * ldpc_G[j][i];
        }
        c[i] = c[i] % 2;
    }

    return c;
}

// Function to perform Binary Phase Shift Keying (BPSK) modulation on an encoded binary sequence c
double *Modulator(int *c, int length)
{
    double *x = (double *)calloc(length, sizeof(double));
    for (int i = 0; i < length; i++)
    {
        x[i] = (c[i] == 0) ? 1.0 : -1.0;
    }

    return x;
}

// Function to simulate the transmission of a modulated sequence x through an additive white Gaussian noise (AWGN) channel
double *Channel(double *x, long *idum, double sigma, int block_length)
{
    //
    int exec_times = block_length / 2;
    double *y = (double *)calloc(block_length, sizeof(double));

    // Loop to add Gaussian noise to each symbol in the modulated sequence
    for (int i = 0; i < exec_times; i++)
    {
        double z[2] = {0.0};
        normal(idum, z, sigma);
        y[2 * i] = x[2 * i] + z[0];
        y[2 * i + 1] = x[2 * i + 1] + z[1];
    }

    // If the block length is odd, add noise to the last symbol
    if (block_length % 2 != 0)
    {
        double z[2] = {0.0};
        normal(idum, z, sigma);
        y[block_length - 1] = x[block_length - 1] + z[0];
    }

    return y;
}

/**
 * Function: decoder
 * Decodes the received sequence `y` using LDPC (Low-Density Parity-Check) with Sum-Product Algorithm.
 *
 * Parameters:
 *  - check_to_variable_msgs: Messages from check nodes to variable nodes.
 *  - variable_to_check_msgs: Messages from variable nodes to check nodes.
 *  - node_connection_matrix: Matrix indicating node connections.
 *  - y: Received sequence.
 *  - noise_variance: Channel's noise variance.
 *  - block_length: Codeword's block length.
 *  - dimension: Codeword's dimension.
 *  - row_weight: LDPC matrix's row weight.
 *  - iteration_max: Maximum iterations for decoding.
 *
 * Returns:
 *  - Decoded codeword.
 */
int *decoder(double **check_to_variable_msgs, double **variable_to_check_msgs, int **node_connection_matrix, double *y, double noise_variance,
             int block_length, int dimension, int row_weight, int iteration_max)
{
    // Allocate memory for the estimated codeword
    int *estimate_codeword = (int *)calloc(dimension, sizeof(int));

    // Initialize L values based on the received symbols and noise variance
    double L[1023] = {0.0};
    for (int i = 0; i < block_length; i++)
    {
        L[i] = 2 * y[i] / noise_variance;
    }

    // Initialize variable-to-check messages using the L values
    for (int variable_node_index = 0; variable_node_index < block_length; variable_node_index++)
    {
        for (int q_edge_index = 0; q_edge_index < row_weight; q_edge_index++)
        {
            int check_node_index = node_connection_matrix[variable_node_index + block_length][q_edge_index];
            variable_to_check_msgs[variable_node_index][check_node_index] = L[variable_node_index];
        }
    }

    double total_message[1023] = {0.0};
    int x_estimated[1023] = {0};

    // Iterative decoding process
    for (int iteration_index = 0; iteration_index < iteration_max; iteration_index++)
    {

        // Bottom-up phase (Check node processing)
        for (int check_node_index = 0; check_node_index < block_length; check_node_index++)
        {

            // Compute initial messages to be sent to variable nodes
            double *chk_to_var_updates = CHK(variable_to_check_msgs, node_connection_matrix, row_weight, check_node_index);

            // Create two temporary arrays to hold messages
            double chk_to_var_updates_0to15[16] = {0.0}, chk_to_var_updates_16to31[16] = {0.0};

            // Populate the temporary arrays with the computed initial messages
            for (int i = 0; i < 16; i++)
            {
                chk_to_var_updates_0to15[i] = chk_to_var_updates[i];
                chk_to_var_updates_16to31[i] = chk_to_var_updates[i];
            }

            // Combine messages using Sum-Product Algorithm (SPA)
            for (int i = 0; i < 4; i++)
            {
                chk_to_var_updates_0to15[i + 8] = SPA(chk_to_var_updates_0to15[i + 8], chk_to_var_updates_0to15[i + 12]);
                chk_to_var_updates_16to31[i] = SPA(chk_to_var_updates_16to31[i], chk_to_var_updates_16to31[i + 4]);
            }

            for (int i = 0; i < 2; i++)
            {
                chk_to_var_updates_0to15[i + 8] = SPA(chk_to_var_updates_0to15[i + 8], chk_to_var_updates_0to15[i + 10]);
                chk_to_var_updates_16to31[i] = SPA(chk_to_var_updates_16to31[i], chk_to_var_updates_16to31[i + 2]);
            }

            // Further combination of messages
            chk_to_var_updates_0to15[8] = SPA(chk_to_var_updates_0to15[8], chk_to_var_updates_0to15[9]);
            chk_to_var_updates_16to31[0] = SPA(chk_to_var_updates_16to31[0], chk_to_var_updates_16to31[1]);

            // Create backup of the temporary arrays
            double update_r_edge_val_temp_0[16] = {0.0}, update_r_edge_val_temp_1[16] = {0.0};
            for (int i = 0; i < 16; i++)
            {

                update_r_edge_val_temp_0[i] = chk_to_var_updates_0to15[i];
                update_r_edge_val_temp_1[i] = chk_to_var_updates_16to31[i];
            }

            // Process edges with indices less than 16
            for (int r_edge_index = 0; r_edge_index < 16; r_edge_index++)
            {

                // Exclude message from the current edge and update with corresponding value from q_edge
                int variable_node_index = r_edge_index % 2 == 0 ? node_connection_matrix[check_node_index][r_edge_index + 1] : node_connection_matrix[check_node_index][r_edge_index - 1];
                update_r_edge_val_temp_0[r_edge_index >> 1] = variable_to_check_msgs[variable_node_index][check_node_index];

                // Combine messages using a tree-reduction algorithm
                for (int i = 0; i < 4; i++)
                {
                    update_r_edge_val_temp_0[i] = SPA(update_r_edge_val_temp_0[i], update_r_edge_val_temp_0[i + 4]);
                }
                for (int i = 0; i < 2; i++)
                {
                    update_r_edge_val_temp_0[i] = SPA(update_r_edge_val_temp_0[i], update_r_edge_val_temp_0[i + 2]);
                }
                update_r_edge_val_temp_0[0] = SPA(update_r_edge_val_temp_0[0], update_r_edge_val_temp_0[1]);
                update_r_edge_val_temp_0[0] = SPA(update_r_edge_val_temp_0[0], update_r_edge_val_temp_0[8]);

                // Update r_edge matrix
                variable_node_index = node_connection_matrix[check_node_index][r_edge_index];
                check_to_variable_msgs[check_node_index][variable_node_index] = update_r_edge_val_temp_0[0];

                // Restore the original values to the temporary arrays
                for (int i = 0; i < 16; i++)
                {
                    update_r_edge_val_temp_0[i] = chk_to_var_updates_0to15[i];
                }
            }

            // Process edges with indices 16 and above
            for (int r_edge_index = 16; r_edge_index < row_weight; r_edge_index++)
            {

                // Exclude message from the current edge and update with corresponding value from q_edge
                int variable_node_index = r_edge_index % 2 == 0 ? node_connection_matrix[check_node_index][r_edge_index + 1] : node_connection_matrix[check_node_index][r_edge_index - 1];
                update_r_edge_val_temp_1[r_edge_index >> 1] = variable_to_check_msgs[variable_node_index][check_node_index];

                // Combine messages using a tree-reduction algorithm
                for (int i = 0; i < 4; i++)
                {
                    update_r_edge_val_temp_1[i + 8] = SPA(update_r_edge_val_temp_1[i + 8], update_r_edge_val_temp_1[i + 12]);
                }
                for (int i = 0; i < 2; i++)
                {
                    update_r_edge_val_temp_1[i + 8] = SPA(update_r_edge_val_temp_1[i + 8], update_r_edge_val_temp_1[i + 10]);
                }
                update_r_edge_val_temp_1[8] = SPA(update_r_edge_val_temp_1[8], update_r_edge_val_temp_1[9]);
                update_r_edge_val_temp_1[0] = SPA(update_r_edge_val_temp_1[0], update_r_edge_val_temp_1[8]);

                // Update r_edge matrix
                variable_node_index = node_connection_matrix[check_node_index][r_edge_index];
                check_to_variable_msgs[check_node_index][variable_node_index] = update_r_edge_val_temp_1[0];

                // Restore the original values to the temporary arrays
                for (int i = 0; i < 16; i++)
                {
                    update_r_edge_val_temp_1[i] = chk_to_var_updates_16to31[i];
                }
            }

            // Release the dynamically allocated memory
            free(chk_to_var_updates);
        }

        // Top-down phase (Variable node processing)
        for (int variable_node_index = 0; variable_node_index < block_length; variable_node_index++)
        {
            // Compute the combined message from all connected check nodes to the current variable node
            double val = VAR(check_to_variable_msgs, node_connection_matrix, row_weight, variable_node_index, block_length);

            // Update the messages sent from the variable node to each connected check node
            for (int q_edge_index = 0; q_edge_index < row_weight; q_edge_index++)
            {
                int check_node_index = node_connection_matrix[variable_node_index + block_length][q_edge_index];

                // Update the message for this edge. The new message is a combination of the intrinsic information (L-values)
                // and the extrinsic information from all other connected check nodes (combined_message), except the current check node.
                variable_to_check_msgs[variable_node_index][check_node_index] = L[variable_node_index] + val - check_to_variable_msgs[check_node_index][variable_node_index];
            }
        }

        // Termination phase
        for (int variable_node_index = 0; variable_node_index < block_length; variable_node_index++)
        {
            // Initialize with the log-likelihood ratio (LLR) based on the received sequence
            total_message[variable_node_index] = L[variable_node_index];

            // Sum up the messages received from all connected check nodes for the current variable node
            for (int r_edge_index = 0; r_edge_index < row_weight; r_edge_index++)
            {
                int check_node_index = node_connection_matrix[variable_node_index + block_length][r_edge_index];
                total_message[variable_node_index] += check_to_variable_msgs[check_node_index][variable_node_index];
            }

            // Decide the estimated bit value (0 or 1) for the current variable node, based on the sign of its total message
            x_estimated[variable_node_index] = (total_message[variable_node_index] > 0.0) ? 0 : 1;
        }

        // Syndrome check phase
        int zero_counter = 0, cols = block_length, rows = block_length, syndrome[1023] = {0};
        for (int i = 0; i < rows; i++)
        {
            // Compute the syndrome
            for (int j = 0; j < row_weight; j++)
            {
                int index = node_connection_matrix[i][j];
                syndrome[i] += x_estimated[index];
            }
            syndrome[i] = syndrome[i] % 2;

            // Count the number of zeros in the syndrome
            if (syndrome[i] == 0)
            {
                zero_counter++;
            }
            else
            {
                break;
            }
        }

        // If all entries in the syndrome are zeros, x_estimated is a valid codeword, terminate the iterations
        if (zero_counter == cols)
        {
            break;
        }
    }

    // Copy the estimated codeword to the output array
    for (int i = 0; i < dimension; i++)
    {
        estimate_codeword[i] = x_estimated[i];
    }

    return estimate_codeword;
}

// CHK function: This function processes messages sent from variable nodes to a specific check node
// identified by 'check_node_index'. It combines pairs of messages using the Sum-Product Algorithm (SPA).
double *CHK(double **variable_to_check_msgs, int **node_connection_matrix, int row_weight, int check_node_index)
{
    // Allocate memory for combined messages.
    double *val = (double *)calloc(row_weight / 2, sizeof(double));

    // Iterate over pairs of variable nodes connected to the check node.
    for (int i = 0; i < row_weight / 2; i++)
    {

        // Using ldpc_H, fetch the indices of the pair of variable nodes
        // that are connected to the check node 'check_node_index'.
        int x_index0 = node_connection_matrix[check_node_index][2 * i];
        int x_index1 = node_connection_matrix[check_node_index][2 * i + 1];

        // Retrieve the messages sent from these variable nodes to the check node.
        double L0 = variable_to_check_msgs[x_index0][check_node_index];
        double L1 = variable_to_check_msgs[x_index1][check_node_index];

        // Combine the retrieved messages using the Sum-Product Algorithm (SPA).
        val[i] = SPA(L0, L1);
    }

    // Return the array containing the combined messages.
    return val;
}

// Function: SPA
// Sum-Product Algorithm (SPA) function, which is used for combining messages in the LDPC decoder
inline double SPA(double L1, double L2)
{
    // Compute the delta term (small correction)
    double delta = log((1 + exp(-1 * fabs(L1 + L2))) / (1 + exp(-1 * fabs(L1 - L2))));

    // Compute the minimum of the magnitudes of L1 and L2
    double min_val = fmin(fabs(L1), fabs(L2));

    // Return the combined log-likelihood ratio (LLR) based on the signs of L1 and L2
    return (L1 * L2 > 0.0) ? min_val + delta : -min_val + delta;
    // return (L1 * L2 > 0.0) ? min_val : -min_val;
}

// VAR function: Computes the outgoing message from a specified variable node
// identified by 'variable_node_index'. The function aggregates messages
// coming from all check nodes connected to this variable node.
double VAR(double **check_to_variable_msgs, int **node_connection_matrix, int row_weight, int variable_node_index, int block_length)
{

    double val = 0.0;

    // Iterate over each connection of the variable node and accumulate the messages
    // from the connected check nodes.
    for (int r_edge_index = 0; r_edge_index < row_weight; r_edge_index++)
    {
        int check_node_index = node_connection_matrix[variable_node_index + block_length][r_edge_index];
        val += check_to_variable_msgs[check_node_index][variable_node_index];
    }

    return val;
}

// Generate a sequence of length 'truncation_length' using a specified rule
int *generated_sequence(int truncation_length)
{
    int *u_tx = (int *)malloc((truncation_length) * sizeof(int));
    // Set the first element to 1
    u_tx[0] = 1;
    for (int i = 1; i < (truncation_length); i++)
        u_tx[i] = 0;

    // Apply the specified rule to generate the sequence starting from the 7th element
    for (int i = 6; i < truncation_length; i++)
        u_tx[i] = (u_tx[i - 6] + u_tx[i - 5]) % 2;

    // Return the generated sequence
    return u_tx;
}

// Function to calculate noise variance given the Signal-to-Noise Ratio (SNR)
double noise_variance(double SNR)
{
    double R = 781.0 / 1023.0;
    return 1 / (2 * R * pow(10, SNR / 10));
}

// Function to initialize the random number generator with a given seed
long *rand_init(long *idum, long seed)
{
    *idum = seed;
    return idum;
}

// Function to generate a uniform random number in the range (0, 1) using a linear congruential generator
double ran1(long *idum)
{
    int j;
    long k;
    static long iy = 0;
    static long iv[NTAB];
    double temp;

    if (*idum <= 0 || !iy)
    {
        if (-(*idum) < 1)
            *idum = 1;
        else
            *idum = -(*idum);
        for (j = NTAB + 7; j >= 0; j--)
        {
            k = (*idum) / IQ;
            *idum = IA * (*idum - k * IQ) - IR * k;
            if (*idum < 0)
                *idum += IM;
            if (j < NTAB)
                iv[j] = *idum;
        }
        iy = iv[0];
    }
    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;
    if (*idum < 0)
        *idum += IM;
    j = iy / NDIV;
    iy = iv[j];
    iv[j] = *idum;

    // Normalize and return the random number
    if ((temp = AM * iy) > RNMX)
        return RNMX;
    else
        return temp;
}

// Function to generate two normally distributed random numbers using the Box-Muller transform
void normal(long *idum, double *n, double sigma)
{
    double x1, x2, s;

    // Generate two independent uniform random variables in the range (-1, 1)
    do
    {
        x1 = ran1(idum);
        x2 = ran1(idum);
        x1 = 2.0 * x1 - 1.0;
        x2 = 2.0 * x2 - 1.0;
        s = pow(x1, 2) + pow(x2, 2); // s is the square of the distance from the origin
    } while (s >= 1.0);              // Repeat if the point is outside the unit circle

    // Use the Box-Muller transform to generate two normally distributed numbers
    n[0] = sigma * x1 * sqrt((-2) * log(s) / s);
    n[1] = sigma * x2 * sqrt((-2) * log(s) / s);
}
