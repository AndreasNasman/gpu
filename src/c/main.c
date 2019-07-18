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

void read_file(FILE *fpointer, const char *DELIMITER, int n, Galaxy *galaxy_set);
__global__ void build_histograms(GalaxySet real, GalaxySet random, int *DD_histogram, int *DR_histogram, int *RR_histogram, int n);

int main()
{
    // REAL GALAXIES FILE
    FILE *fpointer = fopen("./input-data/real-galaxies.txt", "r");
    const char *DELIMITER = "\t";

    // Reads number of lines to process (defined on the first line of the file).
    char line[LINE_LENGTH];
    fgets(line, LINE_LENGTH, fpointer);
    const int NUMBER_OF_LINES = atoi(line);

    GalaxySet real;
    cudaMallocManaged(&real.galaxies, NUMBER_OF_LINES * sizeof(Galaxy));

    read_file(fpointer, DELIMITER, NUMBER_OF_LINES, real.galaxies);

    // RANDOM GALAXIES FILE
    fpointer = fopen("./input-data/random-galaxies.txt", "r");
    DELIMITER = " ";

    // Checks that number of lines is equal in both files.
    fgets(line, LINE_LENGTH, fpointer);
    if (NUMBER_OF_LINES != atoi(line))
    {
        printf("Both files should have equal number of lines!");
        return 1;
    }

    GalaxySet random;
    cudaMallocManaged(&random.galaxies, NUMBER_OF_LINES * sizeof(Galaxy));

    read_file(fpointer, DELIMITER, NUMBER_OF_LINES, random.galaxies);

    fclose(fpointer);

    // BUILDING OF HISTOGRAMS
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

    /* DEBUG */

    long int totalSizeDD = 0;
    long int totalSizeDR = 0;
    long int totalSizeRR = 0;
    for (int i = 0; i < NUMBER_OF_BINS; i += 1)
    {
        totalSizeDD += DD_histogram[i];
        totalSizeDR += DR_histogram[i];
        totalSizeRR += RR_histogram[i];
    }
    printf("\nTotal size DD: %li\n", totalSizeDD);
    printf("Total size DR: %li\n", totalSizeDR);
    printf("Total size RR: %li\n", totalSizeRR);

    /* DEBUG */

    // CLEAN UP
    cudaFree(real.galaxies);
    cudaFree(random.galaxies);
    cudaFree(DD_histogram);
    cudaFree(DR_histogram);
    cudaFree(RR_histogram);

    printf("Done!\n\n");

    return 0;
}

double arcminutes_to_radians(double arcminute_value)
{
    return (M_PI * arcminute_value) / (60 * 180);
}

void read_file(FILE *fpointer, const char *DELIMITER, int n, Galaxy *galaxies)
{
    char line[LINE_LENGTH];
    const int DECLINATION_INDEX = 0;
    const int RIGHT_ASCENSION_INDEX = 1;

    for (int i = 0; i < n; i += 1)
    {
        fgets(line, LINE_LENGTH, fpointer);

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

__device__ double angle_between_two_galaxies(Galaxy first_galaxy, Galaxy second_galaxy)
{
    return acos(
        sin(first_galaxy.declination) * sin(second_galaxy.declination) +
        cos(first_galaxy.declination) * cos(second_galaxy.declination) *
            cos(first_galaxy.right_ascension - second_galaxy.right_ascension));
}

__device__ double radians_to_degrees(double radian_value)
{
    return radian_value * (180 / M_PI);
}

__device__ void update_bin(int *bin, double angle)
{
    int index = floor(radians_to_degrees(angle) / BIN_WIDTH);
    atomicAdd(&bin[index], 1);
}

__global__ void build_histograms(GalaxySet real, GalaxySet random, int *DD_histogram, int *DR_histogram, int *RR_histogram, int n)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    double angle;
    for (int i = 0; i < n; i += 1)
        for (int j = index; j < n; j += stride)
        {
            angle = angle_between_two_galaxies(real.galaxies[i], random.galaxies[j]);
            update_bin(DR_histogram, angle);

            if (j > i)
            {
                angle = angle_between_two_galaxies(real.galaxies[i], real.galaxies[j]);
                update_bin(DD_histogram, angle);

                angle = angle_between_two_galaxies(random.galaxies[i], random.galaxies[j]);
                update_bin(RR_histogram, angle);
            }
        }
}
