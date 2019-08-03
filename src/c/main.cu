#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define BIN_WIDTH 0.25
#define BLOCK_DIM 256
#define LINE_LENGTH 30

typedef struct Galaxy
{
    float declination;
    float right_ascension;
} Galaxy;

typedef struct GalaxySet
{
    Galaxy *galaxies;
} GalaxySet;

__global__ void build_histograms(GalaxySet real, GalaxySet random, int *DD_histogram, int *DR_histogram, int *RR_histogram, int n);
__global__ void galaxy_distribution(int *DD_histogram, int *DR_histogram, int *RR_histogram, float *distribution, int n);
void read_file(FILE *file_pointer, const char *DELIMITER, Galaxy *galaxy_set, int n);
void write_file_int(FILE *file_pointer, int *content, int n);
void write_file_float(FILE *file_pointer, float *content, int n);

int main()
{
    /* READING REAL GALAXIES FILE */
    FILE *file_pointer = fopen("./input-data/real-galaxies.txt", "r");
    const char *DELIMITER = "\t";

    // Reads number of lines to process (defined on the first line of the file).
    char line[LINE_LENGTH];
    fgets(line, LINE_LENGTH, file_pointer);
    const int NUMBER_OF_LINES = atoi(line);

    GalaxySet real;
    cudaMallocManaged(&real.galaxies, NUMBER_OF_LINES * sizeof(Galaxy));

    read_file(file_pointer, DELIMITER, real.galaxies, NUMBER_OF_LINES);

    /* READING RANDOM GALAXIES FILE */
    file_pointer = fopen("./input-data/random-galaxies.txt", "r");
    DELIMITER = " ";

    // Checks that number of lines is equal in both files.
    fgets(line, LINE_LENGTH, file_pointer);
    if (NUMBER_OF_LINES != atoi(line))
    {
        printf("Both files should have equal number of lines!");
        return 1;
    }

    GalaxySet random;
    cudaMallocManaged(&random.galaxies, NUMBER_OF_LINES * sizeof(Galaxy));

    read_file(file_pointer, DELIMITER, random.galaxies, NUMBER_OF_LINES);

    /* BUILDING HISTOGRAMS */
    const int COVERAGE = 180; // degrees
    const int NUMBER_OF_BINS = COVERAGE * (1 / BIN_WIDTH);

    // Defines number of blocks to use.
    const int GRID_DIM = ceil(NUMBER_OF_LINES / (float)BLOCK_DIM);

    int *DD_histogram_temp, *DR_histogram_temp, *RR_histogram_temp;
    cudaMallocManaged(&DD_histogram_temp, NUMBER_OF_BINS * GRID_DIM * sizeof(int));
    cudaMallocManaged(&DR_histogram_temp, NUMBER_OF_BINS * GRID_DIM * sizeof(int));
    cudaMallocManaged(&RR_histogram_temp, NUMBER_OF_BINS * GRID_DIM * sizeof(int));

    build_histograms<<<GRID_DIM, BLOCK_DIM>>>(real, random, DD_histogram_temp, DR_histogram_temp, RR_histogram_temp, NUMBER_OF_LINES);
    cudaDeviceSynchronize();

    int *DD_histogram, *DR_histogram, *RR_histogram;
    cudaMallocManaged(&DD_histogram, NUMBER_OF_BINS * sizeof(int));
    cudaMallocManaged(&DR_histogram, NUMBER_OF_BINS * sizeof(int));
    cudaMallocManaged(&RR_histogram, NUMBER_OF_BINS * sizeof(int));

    for (int i = 0; i < NUMBER_OF_BINS * GRID_DIM; i += 1)
    {
        DD_histogram[i % NUMBER_OF_BINS] += DD_histogram_temp[i];
        DR_histogram[i % NUMBER_OF_BINS] += DR_histogram_temp[i];
        RR_histogram[i % NUMBER_OF_BINS] += RR_histogram_temp[i];
    }

    /* DETERMINING DISTRIBUTION */
    float *distribution;
    cudaMallocManaged(&distribution, NUMBER_OF_BINS * sizeof(float));

    galaxy_distribution<<<GRID_DIM, BLOCK_DIM>>>(DD_histogram, DR_histogram, RR_histogram, distribution, NUMBER_OF_BINS);
    cudaDeviceSynchronize();

    /* WRITING RESULTS TO FILE */
    system("mkdir -p results");

    file_pointer = fopen("results/DD_histogram.txt", "w");
    write_file_int(file_pointer, DD_histogram, NUMBER_OF_BINS);

    file_pointer = fopen("results/RR_histogram.txt", "w");
    write_file_int(file_pointer, RR_histogram, NUMBER_OF_BINS);

    file_pointer = fopen("results/Distribution.txt", "w");
    write_file_float(file_pointer, distribution, NUMBER_OF_BINS);

    /* CLEAN UP */
    fclose(file_pointer);

    cudaFree(real.galaxies);
    cudaFree(random.galaxies);
    cudaFree(DD_histogram);
    cudaFree(DR_histogram);
    cudaFree(RR_histogram);

    printf("Done!\n\n");

    return 0;
}

__device__ float angle_between_two_galaxies(Galaxy first_galaxy, Galaxy second_galaxy)
{
    float x = sinf(first_galaxy.declination) * sinf(second_galaxy.declination) +
              cosf(first_galaxy.declination) * cosf(second_galaxy.declination) *
                  cosf(first_galaxy.right_ascension - second_galaxy.right_ascension);

    // Checks that x is within the boundaries of [-1.0f, 1.0f].
    return acosf(fmin(1.0f, fmax(-1.0f, x)));
}

__device__ float radians_to_degrees(float radian_value)
{
    return radian_value * (180 / M_PI);
}

__device__ void update_bin(int *bin, float angle, int incrementor)
{
    int index = floor(radians_to_degrees(angle) / BIN_WIDTH);
    atomicAdd(&bin[index], incrementor);
}

__global__ void build_histograms(GalaxySet real, GalaxySet random, int *DD_histogram, int *DR_histogram, int *RR_histogram, int n)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    __shared__ int shared_DD_histogram[720];
    __shared__ int shared_DR_histogram[720];
    __shared__ int shared_RR_histogram[720];
    for (int i = 0; i < 720; i += 1)
    {
        shared_DD_histogram[i] = 0;
        shared_DR_histogram[i] = 0;
        shared_RR_histogram[i] = 0;
    }
    __syncthreads();

    float angle;
    for (int i = 0; i < n; i += 1)
        for (int j = index; j < n; j += stride)
        {
            // Every pair of real-random galaxy is compared.
            angle = angle_between_two_galaxies(real.galaxies[i], random.galaxies[j]);
            update_bin(shared_DR_histogram, angle, 1);

            // Real-real and random-random galaxy pairs are only compared from the same starting index forward.
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
                angle = angle_between_two_galaxies(real.galaxies[i], real.galaxies[j]);
                update_bin(shared_DD_histogram, angle, 2);

                angle = angle_between_two_galaxies(random.galaxies[i], random.galaxies[j]);
                update_bin(shared_RR_histogram, angle, 2);
            }
        }
    __syncthreads();

    for (int i = 0; i < 720; i += 1)
    {
        DD_histogram[720 * blockIdx.x + i] = shared_DD_histogram[i];
        DR_histogram[720 * blockIdx.x + i] = shared_DR_histogram[i];
        RR_histogram[720 * blockIdx.x + i] = shared_RR_histogram[i];
    }
}

__global__ void galaxy_distribution(int *DD_histogram, int *DR_histogram, int *RR_histogram, float *distribution, int n)
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

float arcminutes_to_radians(float arcminute_value)
{
    return (M_PI * arcminute_value) / (60 * 180);
}

void read_file(FILE *file_pointer, const char *DELIMITER, Galaxy *galaxies, int n)
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
                galaxies[i].declination = arcminutes_to_radians(arcminute_value);
            else if (index == RIGHT_ASCENSION_INDEX)
                galaxies[i].right_ascension = arcminutes_to_radians(arcminute_value);

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
