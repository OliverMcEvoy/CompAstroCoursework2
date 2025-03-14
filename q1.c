#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

/**
 * Get a random number within a range
 * @param random_parameters The parameters that are used for the Linear congreguntion function.
 * @param lower_bound Lower bound of the desired random number range.
 * @param upper_bound Upper bound of the desired random number range.
 * @return A random double within the specified range [lower_bound, upper_bound].
 */
double random_number_with_mapping(int64_t *random_parameters, double lower_bound, double upper_bound)
{
    // Extract values, not neccessary but I think it improves readility.
    int64_t *seed = &random_parameters[0];
    int64_t *a = &random_parameters[1];
    int64_t *c = &random_parameters[2];
    int64_t *m = &random_parameters[3];

    // Generate and update the seed using linear congruential generator formula.
    (*seed) = ((*a) * (*seed) + (*c)) % (*m);

    // Map the generated number to the desired range.
    double normalised = (double)(*seed) / (double)(*m);
    return (normalised * (upper_bound - lower_bound)) + lower_bound;
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
void rejection_sampling(int32_t num_samples, int64_t *random_parameters, double to_low, double to_high, FILE *rejection_method_results)
{
    for (int32_t i = 0; i < num_samples; i++)
    {
        double mu = random_number_with_mapping(random_parameters, -1, 1);
        double prob_mu = probability_of_mu(mu);

        double acception_threshold = random_number_with_mapping(random_parameters, 0, 1);

        // Ensure that enough accepted samples are gotten
        if (acception_threshold > prob_mu)
            --i;

        fprintf(rejection_method_results, "%f,%s\n", mu, acception_threshold < prob_mu ? "accepted" : "rejected");
    }
}

/**
 * Function to generate a lookup array to match a random y value to a given mu value
 */
double *generate_y_to_mu_lookup_array(int32_t length_of_array)
{

    // Open the file
    FILE *y_to_mu_lookup_csv;
    y_to_mu_lookup_csv = fopen("y_to_mu_lookup.csv", "r");

    double increment = (double)2 / ((double)length_of_array - 1);
    // y in the first column and then mu in the second column.

    // double inverse_distribution[number_of_samples][2];
    double *y_to_mu_lookup_array;
    y_to_mu_lookup_array = (double *)malloc(length_of_array * 2 * sizeof(double));

    // To save some processing power if the look up table is already created we can just read it in.
    if (y_to_mu_lookup_csv != NULL)
    {
        for (int32_t i = 0; i < length_of_array; i++)
        {
            // Elemnent in an array is accessed via row index * total columns + column index you want, it will always be expressed as this to limit ambiguity
            // int count = fscanf(inverse_distribution_file, "%lf,%lf", inverse_distribution[i * 2 + 0], inverse_distribution[i * 2 + 1]);

            int count = fscanf(y_to_mu_lookup_csv, "%lf,%lf", &y_to_mu_lookup_array[i * 2 + 0], &y_to_mu_lookup_array[i * 2 + 1]);
            // int count = 2;
            // I'll delete the file and run it again if theres an issue with the file, TODO: implement something that if this happens more than once an error is raised.
            if (count != 2)
            {
                printf("ERROR: fscanf matching error line %d\n", i);
                // Note I could try open the file again in write more but the specific behavior depends on the operating system and that I want to avoid.
                fclose(y_to_mu_lookup_csv);
                remove("y_to_mu_lookup.csv");

                free(y_to_mu_lookup_array);
                return generate_y_to_mu_lookup_array(length_of_array);
            }
        }
        fclose(y_to_mu_lookup_csv);
        // If no errors horray! and we can return the lookup array.
        return y_to_mu_lookup_array;
    }

    else
    {
        y_to_mu_lookup_csv = fopen("y_to_mu_lookup.csv", "w");
    }

    // First iteration
    // Elemnent in an array is accessed via row index * total columns + column index you want, it will always be expressed as this to limit ambiguity
    y_to_mu_lookup_array[0] = -1;
    y_to_mu_lookup_array[0 * 2 + 1] = (pow(y_to_mu_lookup_array[0], 3) + 3 * y_to_mu_lookup_array[0] + 4) / (double)8;
    fprintf(y_to_mu_lookup_csv, "%lf,%lf\n", y_to_mu_lookup_array[0 * 2 + 0], y_to_mu_lookup_array[0 * 2 + 1]);

    printf("generating the inverse distribution...\n");
    for (int32_t i = 1; i < length_of_array; i++)
    {
        // Generate the next mu
        y_to_mu_lookup_array[i * 2 + 0] = y_to_mu_lookup_array[(i - 1) * 2 + 0] + increment;
        // Get the corrosponding y value
        y_to_mu_lookup_array[i * 2 + 1] = (pow(y_to_mu_lookup_array[i * 2 + 0], 3) + 3 * y_to_mu_lookup_array[i * 2 + 0] + 4) / (double)8;

        fprintf(y_to_mu_lookup_csv, "%lf,%lf\n", y_to_mu_lookup_array[i * 2 + 0], y_to_mu_lookup_array[i * 2 + 1]);
    }
    printf("done generating the inverse distribution \n");
    fclose(y_to_mu_lookup_csv);

    return y_to_mu_lookup_array;
}

double get_mu_from_y(double y, const double *lookup_table, int32_t size_of_table)
{
    int32_t low = 0, high = size_of_table - 1;

    // Binary search for the closest y values.
    while (low <= high)
    {
        int32_t mid = low + (high - low) / 2;
        double mid_y = lookup_table[2 * mid + 1];

        if (mid_y < y)
        {
            low = mid + 1;
        }
        else
        {
            high = mid - 1;
        }
    }

    // Take the midpoint of the two closest mu values
    return (lookup_table[2 * high] + lookup_table[2 * low]) * 0.5;
}

void direct_mapping(int32_t number_of_samples, int64_t *random_parameters)
{

    int32_t size_of_table = 1000;
    double *y_to_mu_lookup_array = generate_y_to_mu_lookup_array(size_of_table);

    double y, mu;
    FILE *efficient_results;
    if ((efficient_results = fopen("efficient_results.csv", "w")) == NULL)
    {
        printf("Error opening file!\n");
    }

    for (int32_t i = 0; i < number_of_samples; i++)
    {

        y = random_number_with_mapping(random_parameters, 0, 1);

        // As this is a 1D array, the accessing of the array is simply the row index i
        mu = get_mu_from_y(y, y_to_mu_lookup_array, size_of_table);

        fprintf(efficient_results, "%f,%f\n", mu, y);
    }

    fclose(efficient_results);
    free(y_to_mu_lookup_array);
}

int main()
{
    int64_t seed = 1, a = 16807, c = 0, m = 2147483647;
    int64_t rand_parameters[] = {seed, a, c, m};
    int32_t number_of_samples = 50000;

    FILE *file_ptr;
    if ((file_ptr = fopen("results.csv", "w")) == NULL)
    {
        printf("Error opening file!\n");
        return 1;
    }

    // Measure time for rejection_sampling
    clock_t start_rejection = clock(); // Start timer
    rejection_sampling(number_of_samples, rand_parameters, 0, 3.14, file_ptr);
    clock_t end_rejection = clock();                                                    // End timer
    double time_rejection = (double)(end_rejection - start_rejection) / CLOCKS_PER_SEC; // Calculate time in seconds
    printf("Time taken for rejection_sampling: %.6f seconds\n", time_rejection);

    // Measure time for direct_mapping
    clock_t start_direct = clock(); // Start timer
    direct_mapping(number_of_samples, rand_parameters);
    clock_t end_direct = clock();                                              // End timer
    double time_direct = (double)(end_direct - start_direct) / CLOCKS_PER_SEC; // Calculate time in seconds
    printf("Time taken for direct_mapping: %.6f seconds\n", time_direct);

    // Compare the two methods
    if (time_rejection < time_direct)
    {
        printf("rejection_sampling was faster by %.6f seconds\n", time_direct - time_rejection);
    }
    else
    {
        printf("direct_mapping was faster by %.6f seconds\n", time_rejection - time_direct);
    }

    fclose(file_ptr);
    return 0;
}
