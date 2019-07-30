#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define BIN_WIDTH 0.25
#define BLOCK_SIZE 256
#define LINE_LENGTH 30

typedef struct Galaxy
{
    double declination;
    double right_ascension;
} Galaxy;

typedef struct GalaxySet
{
    Galaxy *galaxies;
} GalaxySet;

__global__ void build_histograms(GalaxySet real, GalaxySet random, int *DD_histogram, int *DR_histogram, int *RR_histogram, int n);
__global__ void galaxy_distribution(int *DD_histogram, int *DR_histogram, int *RR_histogram, int n, double *distribution);
void read_file(FILE *filePointer, const char *DELIMITER, int n, Galaxy *galaxy_set);
void write_file_int(FILE *filePointer, int *content, int n);
void write_file_double(FILE *filePointer, double *content, int n);

int main()
{
    /* READING REAL GALAXIES FILE */
    FILE *filePointer = fopen("./input-data/real-galaxies.txt", "r");
    const char *DELIMITER = "\t";

    // Reads number of lines to process (defined on the first line of the file).
    char line[LINE_LENGTH];
    fgets(line, LINE_LENGTH, filePointer);
    const int NUMBER_OF_LINES = atoi(line);

    GalaxySet real;
    cudaMallocManaged(&real.galaxies, NUMBER_OF_LINES * sizeof(Galaxy));

    read_file(filePointer, DELIMITER, NUMBER_OF_LINES, real.galaxies);

    /* READING RANDOM GALAXIES FILE */
    filePointer = fopen("./input-data/random-galaxies.txt", "r");
    DELIMITER = " ";

    // Checks that number of lines is equal in both files.
    fgets(line, LINE_LENGTH, filePointer);
    if (NUMBER_OF_LINES != atoi(line))
    {
        printf("Both files should have equal number of lines!");
        return 1;
    }

    GalaxySet random;
    cudaMallocManaged(&random.galaxies, NUMBER_OF_LINES * sizeof(Galaxy));

    read_file(filePointer, DELIMITER, NUMBER_OF_LINES, random.galaxies);

    /* BUILDING HISTOGRAMS */
    const int COVERAGE = 180; // degrees
    const int NUMBER_OF_BINS = COVERAGE * (1 / BIN_WIDTH);

    // Defines number of blocks to use.
    const int NUMBER_OF_BLOCKS = (NUMBER_OF_LINES + BLOCK_SIZE - 1) / BLOCK_SIZE;

    int *DD_histogram, *DR_histogram, *RR_histogram;
    cudaMallocManaged(&DD_histogram, NUMBER_OF_BINS * sizeof(int));
    cudaMallocManaged(&DR_histogram, NUMBER_OF_BINS * sizeof(int));
    cudaMallocManaged(&RR_histogram, NUMBER_OF_BINS * sizeof(int));

    build_histograms<<<NUMBER_OF_BLOCKS, BLOCK_SIZE>>>(real, random, DD_histogram, DR_histogram, RR_histogram, NUMBER_OF_LINES);

    cudaDeviceSynchronize();

    /* DETERMINING DISTRIBUTION */
    double *distribution;
    cudaMallocManaged(&distribution, NUMBER_OF_BINS * sizeof(double));
    galaxy_distribution<<<NUMBER_OF_BLOCKS, BLOCK_SIZE>>>(DD_histogram, DR_histogram, RR_histogram, NUMBER_OF_BINS, distribution);

    cudaDeviceSynchronize();

    /* WRITING RESULTS TO FILE */
    system("mkdir -p results");

    filePointer = fopen("results/DD_histogram.txt", "w");
    write_file_int(filePointer, DD_histogram, NUMBER_OF_BINS);

    filePointer = fopen("results/RR_histogram.txt", "w");
    write_file_int(filePointer, RR_histogram, NUMBER_OF_BINS);

    filePointer = fopen("results/Distribution.txt", "w");
    write_file_double(filePointer, distribution, NUMBER_OF_BINS);

    /* CLEAN UP */
    fclose(filePointer);

    cudaFree(real.galaxies);
    cudaFree(random.galaxies);
    cudaFree(DD_histogram);
    cudaFree(DR_histogram);
    cudaFree(RR_histogram);

    printf("Done!\n\n");

    return 0;
}

__device__ double angle_between_two_galaxies(Galaxy first_galaxy, Galaxy second_galaxy)
{
    // Checks for duplications (exists in the real galaxy file), which would otherwise make the algorithm return 'null' or 'nan'.
    if (first_galaxy.declination == second_galaxy.declination && first_galaxy.right_ascension == second_galaxy.right_ascension)
        return 0;

    return acos(
        sin(first_galaxy.declination) * sin(second_galaxy.declination) +
        cos(first_galaxy.declination) * cos(second_galaxy.declination) *
            cos(first_galaxy.right_ascension - second_galaxy.right_ascension));
}

__device__ double radians_to_degrees(double radian_value)
{
    return radian_value * (180 / M_PI);
}

__device__ void update_bin(int *bin, double angle, int incrementor)
{
    int index = floor(radians_to_degrees(angle) / BIN_WIDTH);
    atomicAdd(&bin[index], incrementor);
}

__global__ void build_histograms(GalaxySet real, GalaxySet random, int *DD_histogram, int *DR_histogram, int *RR_histogram, int n)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    double angle;
    for (int i = 0; i < n; i += 1)
        for (int j = index; j < n; j += stride)
        {
            // Every pair of real-random galaxy is compared.
            angle = angle_between_two_galaxies(real.galaxies[i], random.galaxies[j]);
            update_bin(DR_histogram, angle, 1);

            // Real-real and random-random galaxy pairs are only compared from the same starting index forward.
            // If both indexes are the same, the relevant bin is incremented by one.
            if (j == i)
            {
                angle = 0;
                update_bin(DD_histogram, angle, 1);
                update_bin(RR_histogram, angle, 1);
            }
            // When one of the indexes is greater, the relevant bin is incremented by two.
            // This is the same as doing the comparison twice, thus saving execution time.
            else if (j > i)
            {
                angle = angle_between_two_galaxies(real.galaxies[i], real.galaxies[j]);
                update_bin(DD_histogram, angle, 2);

                angle = angle_between_two_galaxies(random.galaxies[i], random.galaxies[j]);
                update_bin(RR_histogram, angle, 2);
            }
        }
}

__global__ void galaxy_distribution(int *DD_histogram, int *DR_histogram, int *RR_histogram, int n, double *distribution)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < n; i += stride)
    {
        if (RR_histogram[i] == 0)
            continue;

        distribution[i] = (double)(DD_histogram[i] - 2 * DR_histogram[i] + RR_histogram[i]) / RR_histogram[i];
    }
}

double arcminutes_to_radians(double arcminute_value)
{
    return (M_PI * arcminute_value) / (60 * 180);
}

void read_file(FILE *filePointer, const char *DELIMITER, int n, Galaxy *galaxies)
{
    char line[LINE_LENGTH];
    const int DECLINATION_INDEX = 1;
    const int RIGHT_ASCENSION_INDEX = 0;

    for (int i = 0; i < n; i += 1)
    {
        fgets(line, LINE_LENGTH, filePointer);

        char *token = strtok(line, DELIMITER);

        int index = 0;
        while (token != NULL)
        {
            double arcminuteValue = atof(token);

            if (index == DECLINATION_INDEX)
            {
                galaxies[i].declination = arcminutes_to_radians(arcminuteValue);
            }
            else if (index == RIGHT_ASCENSION_INDEX)
            {
                galaxies[i].right_ascension = arcminutes_to_radians(arcminuteValue);
            }

            index += 1;
            token = strtok(NULL, DELIMITER);
        }
    }
}

void write_file_int(FILE *filePointer, int *content, int n)
{
    for (int i = 0; i < n; i += 1)
        fprintf(filePointer, "%d\n", content[i]);
}

void write_file_double(FILE *filePointer, double *content, int n)
{
    for (int i = 0; i < n; i += 1)
        fprintf(filePointer, "%f\n", content[i]);
}
