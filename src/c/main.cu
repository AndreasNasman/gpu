#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define BIN_WIDTH 0.25
#define BLOCK_DIM 256
#define COVERAGE 180
#define LINE_LENGTH 30

#define BINS_TOTAL (COVERAGE * (int)(1 / BIN_WIDTH))

typedef struct Galaxy
{
    float declination;
    float declination_cos;
    float declination_sin;
    float right_ascension;
} Galaxy;

__global__ void adjust_galaxy_set(Galaxy *galaxy_set, int n);
__global__ void collect_histograms(Galaxy *D_galaxy_set, Galaxy *R_galaxy_set, int *DD_histogram_collected, int *DR_histogram_collected, int *RR_histogram_collected, int n);
__global__ void measure_galaxy_distribution(int *DD_histogram, int *DR_histogram, int *RR_histogram, float *distribution, int n);
void accumulate_histograms(int *DD_histogram_collected, int *DR_histogram_collected, int *RR_histogram_collected, int *DD_histogram, int *DR_histogram, int *RR_histogram, int n);
void read_file(FILE *file_pointer, const char *DELIMITER, Galaxy *galaxy_set, int n);
void write_file_int(FILE *file_pointer, int *content, int n);
void write_file_float(FILE *file_pointer, float *content, int n);

int main()
{
    /* READING REAL GALAXIES FILE */
    FILE *file_pointer = fopen("./input-data/real-galaxies.txt", "r");
    const char *DELIMITER = "\t";

    // Reads number of galaxies to process (defined on the first line of the file).
    char line[LINE_LENGTH];
    fgets(line, LINE_LENGTH, file_pointer);
    const int GALAXIES_TOTAL = atoi(line);

    Galaxy *real;
    cudaMallocManaged(&real, GALAXIES_TOTAL * sizeof(Galaxy));

    read_file(file_pointer, DELIMITER, real, GALAXIES_TOTAL);

    /* READING RANDOM GALAXIES FILE */
    file_pointer = fopen("./input-data/random-galaxies.txt", "r");
    DELIMITER = " ";

    // Checks that number of galaxies is equal in both files.
    fgets(line, LINE_LENGTH, file_pointer);
    if (GALAXIES_TOTAL != atoi(line))
    {
        printf("Both files should have equal number of galaxies!");
        return 1;
    }

    Galaxy *random;
    cudaMallocManaged(&random, GALAXIES_TOTAL * sizeof(Galaxy));

    read_file(file_pointer, DELIMITER, random, GALAXIES_TOTAL);

    /* ADJUSTING GALAXY SETS */
    int GRID_DIM = ceilf(GALAXIES_TOTAL / (float)BLOCK_DIM);
    adjust_galaxy_set<<<GRID_DIM, BLOCK_DIM>>>(real, GALAXIES_TOTAL);
    adjust_galaxy_set<<<GRID_DIM, BLOCK_DIM>>>(random, GALAXIES_TOTAL);

    /* COLLECTING HISTOGRAMS */
    const int COLLECTED_HISTOGRAM_SIZE = GRID_DIM * BINS_TOTAL;

    int *DD_histogram_collected, *DR_histogram_collected, *RR_histogram_collected;
    cudaMallocManaged(&DD_histogram_collected, COLLECTED_HISTOGRAM_SIZE * sizeof(int));
    cudaMallocManaged(&DR_histogram_collected, COLLECTED_HISTOGRAM_SIZE * sizeof(int));
    cudaMallocManaged(&RR_histogram_collected, COLLECTED_HISTOGRAM_SIZE * sizeof(int));

    collect_histograms<<<GRID_DIM, BLOCK_DIM>>>(real, random, DD_histogram_collected, DR_histogram_collected, RR_histogram_collected, GALAXIES_TOTAL);
    cudaDeviceSynchronize();

    /* ACCUMULATING HISTOGRAMS */
    int *DD_histogram, *DR_histogram, *RR_histogram;
    cudaMallocManaged(&DD_histogram, BINS_TOTAL * sizeof(int));
    cudaMallocManaged(&DR_histogram, BINS_TOTAL * sizeof(int));
    cudaMallocManaged(&RR_histogram, BINS_TOTAL * sizeof(int));

    accumulate_histograms(DD_histogram_collected, DR_histogram_collected, RR_histogram_collected, DD_histogram, DR_histogram, RR_histogram, COLLECTED_HISTOGRAM_SIZE);

    /* DETERMINING DISTRIBUTION */
    float *distribution;
    cudaMallocManaged(&distribution, BINS_TOTAL * sizeof(float));

    GRID_DIM = ceilf(BINS_TOTAL / (float)BLOCK_DIM);
    measure_galaxy_distribution<<<GRID_DIM, BLOCK_DIM>>>(DD_histogram, DR_histogram, RR_histogram, distribution, BINS_TOTAL);
    cudaDeviceSynchronize();

    /* WRITING RESULTS TO FILE */
    system("mkdir -p results");

    file_pointer = fopen("results/DD_histogram.txt", "w");
    write_file_int(file_pointer, DD_histogram, BINS_TOTAL);

    file_pointer = fopen("results/RR_histogram.txt", "w");
    write_file_int(file_pointer, RR_histogram, BINS_TOTAL);

    file_pointer = fopen("results/Distribution.txt", "w");
    write_file_float(file_pointer, distribution, BINS_TOTAL);

    /* CLEAN UP */
    fclose(file_pointer);

    cudaFree(real);
    cudaFree(random);
    cudaFree(DD_histogram);
    cudaFree(DR_histogram);
    cudaFree(RR_histogram);

    printf("Galaxy distribution successfully calculated!\n");

    return 0;
}

__device__ float arcminutes_to_radians(float arcminute_value)
{
    return (M_PI * arcminute_value) / (60 * 180);
}

__global__ void adjust_galaxy_set(Galaxy *galaxy_set, int n)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < n; i += stride)
    {
        float declination = arcminutes_to_radians(galaxy_set[i].declination);
        galaxy_set[i].declination = declination;
        galaxy_set[i].declination_cos = cosf(declination);
        galaxy_set[i].declination_sin = sinf(declination);

        galaxy_set[i].right_ascension = arcminutes_to_radians(galaxy_set[i].right_ascension);
    }
}

__device__ float angle_between_galaxies(Galaxy first_galaxy, Galaxy second_galaxy)
{
    float x = first_galaxy.declination_sin * second_galaxy.declination_sin +
              first_galaxy.declination_cos * second_galaxy.declination_cos *
                  cosf(first_galaxy.right_ascension - second_galaxy.right_ascension);

    // Checks that x is within the boundaries of [-1.0f, 1.0f].
    x = fminf(1.0f, fmaxf(-1.0f, x));

    return acosf(x);
}

__device__ float radians_to_degrees(float radian_value)
{
    return radian_value * (180 / M_PI);
}

__device__ void update_bin(int *bin, float angle, int incrementor)
{
    int index = floorf(radians_to_degrees(angle) / BIN_WIDTH);
    atomicAdd(&bin[index], incrementor);
}

__global__ void collect_histograms(Galaxy *D_galaxy_set, Galaxy *R_galaxy_set, int *DD_histogram_collected, int *DR_histogram_collected, int *RR_histogram_collected, int n)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    // Shared arrays are used to reduce duration of atomic adding when updating results.
    // Each block uses its own set of arrays of histograms, equal to the size of the actual histogram arrays.
    __shared__ int shared_DD_histogram[BINS_TOTAL];
    __shared__ int shared_DR_histogram[BINS_TOTAL];
    __shared__ int shared_RR_histogram[BINS_TOTAL];
    for (int i = threadIdx.x; i < BINS_TOTAL; i += blockDim.x)
    {
        shared_DD_histogram[i] = 0;
        shared_DR_histogram[i] = 0;
        shared_RR_histogram[i] = 0;
    }
    __syncthreads();

    for (int i = 0; i < n; i += 1)
        for (int j = index; j < n; j += stride)
        {
            // Every pair of D-R galaxy is compared.
            float angle = angle_between_galaxies(D_galaxy_set[i], R_galaxy_set[j]);
            update_bin(shared_DR_histogram, angle, 1);

            // D-D and R-R galaxy pairs are only compared from the same starting index forward.
            // If both indexes are the same, the relevant bin is incremented by one.
            if (j == i)
            {
                angle = 0;
                update_bin(shared_DD_histogram, angle, 1);
                update_bin(shared_RR_histogram, angle, 1);
            }
            // When one of the indexes is greater, the relevant bin is incremented by two.
            // This is the same as doing the comparison twice, thus saving execution time.
            else if (j > i)
            {
                angle = angle_between_galaxies(D_galaxy_set[i], D_galaxy_set[j]);
                update_bin(shared_DD_histogram, angle, 2);

                angle = angle_between_galaxies(R_galaxy_set[i], R_galaxy_set[j]);
                update_bin(shared_RR_histogram, angle, 2);
            }
        }
    __syncthreads();

    // The section corresponding to the current block is updated in the collected histogram arrays.
    for (int i = threadIdx.x; i < BINS_TOTAL; i += blockDim.x)
    {
        DD_histogram_collected[blockIdx.x * BINS_TOTAL + i] = shared_DD_histogram[i];
        DR_histogram_collected[blockIdx.x * BINS_TOTAL + i] = shared_DR_histogram[i];
        RR_histogram_collected[blockIdx.x * BINS_TOTAL + i] = shared_RR_histogram[i];
    }
}

__global__ void measure_galaxy_distribution(int *DD_histogram, int *DR_histogram, int *RR_histogram, float *distribution, int n)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < n; i += stride)
    {
        if (RR_histogram[i] == 0)
            continue;

        distribution[i] = (DD_histogram[i] - 2.0f * DR_histogram[i] + RR_histogram[i]) / RR_histogram[i];
    }
}

void accumulate_histograms(int *DD_histogram_collected, int *DR_histogram_collected, int *RR_histogram_collected, int *DD_histogram, int *DR_histogram, int *RR_histogram, int n)
{
    for (int i = 0; i < n; i += 1)
    {
        DD_histogram[i % BINS_TOTAL] += DD_histogram_collected[i];
        DR_histogram[i % BINS_TOTAL] += DR_histogram_collected[i];
        RR_histogram[i % BINS_TOTAL] += RR_histogram_collected[i];
    }
}

void read_file(FILE *file_pointer, const char *DELIMITER, Galaxy *galaxy_set, int n)
{
    char line[LINE_LENGTH];
    const int DECLINATION_INDEX = 1;
    const int RIGHT_ASCENSION_INDEX = 0;

    for (int i = 0; i < n; i += 1)
    {
        fgets(line, LINE_LENGTH, file_pointer);

        char *token = strtok(line, DELIMITER);

        int index = 0;
        while (token != NULL)
        {
            float arcminute_value = atof(token);

            if (index == DECLINATION_INDEX)
                galaxy_set[i].declination = arcminute_value;
            else if (index == RIGHT_ASCENSION_INDEX)
                galaxy_set[i].right_ascension = arcminute_value;

            index += 1;
            token = strtok(NULL, DELIMITER);
        }
    }
}

void write_file_int(FILE *file_pointer, int *content, int n)
{
    for (int i = 0; i < n; i += 1)
        fprintf(file_pointer, "%d\n", content[i]);
}

void write_file_float(FILE *file_pointer, float *content, int n)
{
    for (int i = 0; i < n; i += 1)
        fprintf(file_pointer, "%f\n", content[i]);
}
