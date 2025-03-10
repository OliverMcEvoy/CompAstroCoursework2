#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>

/**
 * Get a random number within a range
 * @param seed A pointer to the seed.
 * @param a The multiplier in the linear congruential generator formula.
 * @param c The increment in the linear congruential generator formula (in this instance it will always be 0, but included for completeness).
 * @param m The modulus in the linear congruential generator formula.
 * @param lower_bound Lower bound of the desired random number range.
 * @param upper_bound Upper bound of the desired random number range.
 * @return A random double within the specified range [lower_bound, upper_bound].
 */
double random_number_with_mapping(int64_t *seed, int16_t a, int16_t c, int32_t m, double lower_bound, double upper_bound)
{
    // Generate and update the seed using linear congruential generator formula.
    (*seed) = (a * (*seed) + c) % m;

    // Map the generated number to the desired range.
    double normalized = (double)(*seed) / m;
    return (normalized * (upper_bound - lower_bound)) + lower_bound;
}

/**
 * Function to calculate the probability of mu according to the given formula.
 * @param mu The input for which we want to calculate the probability.
 * @return The calculated probability to be used wether to accept a given mu or not, see report for derivation
 */
double probability_of_mu(double mu)
{
    // 3 and 8 have to be specified as doubles to avoid the interpretation of the math on them being done first and therefore having a non-double * a double
    return (double)3 / (double)8 * (1 + pow(mu, 2));
}

/**
 * Function to perform rejection sampling according to the given parameters.
 * @param num_samples The number of samples to generate.
 * @param seed Initial seed for random number generator.
 * @param a, c, m Parameters for the random number generator.
 * @param to_low, to_high Mapping range for random numbers.
 */
void rejection_sampling(int32_t num_samples, int64_t *seed, int16_t *a, int16_t *c, int32_t *m, double to_low, double to_high, FILE *file_ptr)
{
    for (int32_t i = 0; i < num_samples; i++)
    {
        double mu = random_number_with_mapping(seed, *a, *c, *m, -1, 1);
        double prob_mu = probability_of_mu(mu);

        double acception_threshold = random_number_with_mapping(seed, *a, *c, *m, 0, 1);

        fprintf(file_ptr, "%f,%s\n", mu, acception_threshold < prob_mu ? "accepted" : "rejected");
    }
}

int main()
{
    int64_t seed = 1;
    int32_t m = 2147483647;
    int16_t a = 16807, c = 0;

    FILE *file_ptr;
    if ((file_ptr = fopen("results.csv", "w")) == NULL)
    {
        printf("Error opening file!\n");
        return 1;
    }

    rejection_sampling((int32_t)20000, &seed, &a, &c, &m, 0, 3.14, file_ptr);

    fclose(file_ptr);
    return 0;
}