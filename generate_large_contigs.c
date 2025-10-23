#include <mkl.h>

#include <stdint.h>
#include <stdio.h>

int main(void)
{
    const char* FILE_NAME = "large_contigs.fasta";
    // Initialize MKL Random Number Generator
    VSLStreamStatePtr stream;
    vslNewStream(&stream, VSL_BRNG_SFMT19937, 0x12345678);

    FILE* fp = fopen(FILE_NAME, "w");
    if (fp == NULL) {
        printf("Error opening file '%s'\n", FILE_NAME);
        return 1;
    }
    const size_t NUM_RAND_NUMBER_CACHE = 1024 * 1024;
    const size_t CONTIG_SIZE = 32ULL * 1024 * 1024 * 1024;
    const size_t NUM_CONTIGS = 4;
    const size_t LINE_LEN = 64;
    // Generate 1M random numbers for later use
    int* rand_nums = (int*)calloc((NUM_RAND_NUMBER_CACHE), sizeof(int));
    size_t rand_offset = 0;
    size_t total_bases = 0;
    char line[LINE_LEN + 1];
    line[LINE_LEN] = '\0';
    for (size_t i = 0; i < NUM_CONTIGS; i++) {
        fprintf(fp, ">contig%d\n", i);
        // Generate 32 GiB ATCGN sequence
        // That is, 64 chars per line for 512 M lines.
        for (size_t j = 0; j < CONTIG_SIZE / LINE_LEN; j++) {
            for (int k = 0; k < LINE_LEN; k++) {
                if (rand_offset == 0) {
                    // Refill random number cache
                    viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, NUM_RAND_NUMBER_CACHE, rand_nums, 0, 5);
                    rand_offset = NUM_RAND_NUMBER_CACHE;
                    printf("Generated %lu/%lu random numbers\n", total_bases, CONTIG_SIZE * NUM_CONTIGS);
                }

                int r = rand_nums[--rand_offset];
                line[k] = "AGCTN"[r];
                total_bases += LINE_LEN;
            }
            fprintf(fp, "%s\n", line);
        }
    }
    vslDeleteStream(&stream);
    fclose(fp);
    return 0;
}
