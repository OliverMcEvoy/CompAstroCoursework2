#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// While I dislike using global varibales it makes sense for pi, I try not to use MP_PI as a refrence for pi as it is not defined in standard C
double pi = 3.141592653589793;

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
    // Note I decided to use a pointer to the seed as if I implement processing over multiple threads, pointers are thread safe while static varibles are typically not
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

double *generate_lookup_array(int32_t params[3], const char *filename,
                              double (*compute_value1)(double),
                              double (*compute_value2)(double))
{
    // Extract the parameters from the array
    int32_t length_of_array = params[0];
    double initial_value = params[1];
    double final_value = params[2];

    // Open the file
    FILE *lookup_file = fopen(filename, "r");

    double increment = (final_value - initial_value) / ((double)length_of_array - 1.0);

    // Determine the number of columns based on whether the second function is provided
    int columns = (compute_value2 == NULL) ? 2 : 3;

    // Allocate memory for the lookup array (2 or 3 columns of data)
    double *lookup_array = (double *)malloc(length_of_array * columns * sizeof(double));
    if (!lookup_array)
    {
        fprintf(stderr, "ERROR: Memory allocation failed for %s\n", filename);
        return NULL;
    }

    // If file exists, read data from it
    if (lookup_file != NULL)
    {
        for (int32_t i = 0; i < length_of_array; i++)
        {
            int count = (columns == 2) ? fscanf(lookup_file, "%lf,%lf", &lookup_array[i * 2], &lookup_array[i * 2 + 1]) : fscanf(lookup_file, "%lf,%lf,%lf", &lookup_array[i * 3], &lookup_array[i * 3 + 1], &lookup_array[i * 3 + 2]);
            if (count != columns)
            {
                printf("ERROR: fscanf matching error in %s on line %d\n", filename, i);
                fclose(lookup_file);
                remove(filename);
                free(lookup_array);
                return generate_lookup_array(params, filename, compute_value1, compute_value2);
            }
        }
        fclose(lookup_file);
        printf("Successfully loaded lookup table from %s\n", filename);
        return lookup_array;
    }

    // Otherwise, create a new file and generate the lookup table
    lookup_file = fopen(filename, "w");
    if (!lookup_file)
    {
        fprintf(stderr, "ERROR: Could not create file %s\n", filename);
        free(lookup_array);
        return NULL;
    }

    // First iteration
    lookup_array[0] = initial_value;
    lookup_array[1] = compute_value1(lookup_array[0]);
    if (columns == 3)
    {
        lookup_array[2] = compute_value2(lookup_array[0]);
        fprintf(lookup_file, "%lf,%lf,%lf\n", lookup_array[0], lookup_array[1], lookup_array[2]);
    }
    else
    {
        fprintf(lookup_file, "%lf,%lf\n", lookup_array[0], lookup_array[1]);
    }

    printf("Generating lookup table and saving to %s...\n", filename);
    for (int32_t i = 1; i < length_of_array; i++)
    {
        // Generate the next x value
        lookup_array[i * columns] = lookup_array[(i - 1) * columns] + increment;
        // Compute the corresponding y value using the first function
        lookup_array[i * columns + 1] = compute_value1(lookup_array[i * columns]);

        // If a second function is provided, compute the corresponding z value
        if (columns == 3)
        {
            lookup_array[i * columns + 2] = compute_value2(lookup_array[i * columns]);
            fprintf(lookup_file, "%lf,%lf,%lf\n", lookup_array[i * columns], lookup_array[i * columns + 1], lookup_array[i * columns + 2]);
        }
        else
        {
            fprintf(lookup_file, "%lf,%lf\n", lookup_array[i * columns], lookup_array[i * columns + 1]);
        }
    }
    printf("Done generating lookup table, saved to %s\n", filename);
    fclose(lookup_file);

    return lookup_array;
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
void rejection_sampling(int32_t num_samples, int64_t *random_parameters, double to_low, double to_high)
{
    FILE *rejection_method_results;
    if ((rejection_method_results = fopen("results.csv", "w")) == NULL)
    {
        printf("Error opening file!\n");
    }

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

    fclose(rejection_method_results);
}

double calculate_y_from_mu(double y)
{
    return (pow(y, 3) + 3 * y + 4) / 8.0;
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

    int32_t lookup_table_params[3] = {2000, -1, 1};
    double *y_to_mu_lookup_array = generate_lookup_array(lookup_table_params, "y_to_mu_lookup.csv", *calculate_y_from_mu, NULL);

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
        mu = get_mu_from_y(y, y_to_mu_lookup_array, lookup_table_params[0]);

        fprintf(efficient_results, "%f,%f\n", mu, y);
    }

    fclose(efficient_results);
    free(y_to_mu_lookup_array);
}

double get_u_from_x(double x)
{
    double r = 10.0;
    return (2.0 / pi) * asin(x / r);
}

double get_y_from_X(double x)
{
    double r = 10.0;
    return sqrt(r * r - x * x);
}

void photon_scattering(int64_t total_photon_count)
{
    // pi hard coded to the max percision a double provides

    // Size of table , start value, final value
    int32_t lookup_table_params[3] = {2000, -10, 10};
    double *displacement_lookup_array = generate_lookup_array(lookup_table_params, "displacement_direction_lookup.csv", *get_u_from_x, *get_y_from_X);

    // Implement the for loop from the lecture notes.
}

int main()
{

    //
    // As seed, a,c,m will all be used in operations with each other I am declaring them as the same type to ensure type safe code.
    int64_t seed = 1, a = 16807, c = 0, m = 2147483647;
    int64_t rand_parameters[] = {seed, a, c, m};

    // Questinon 1.

    int32_t number_of_samples = 500000;
    // Measure time for rejection_sampling
    clock_t start_rejection = clock(); // Start timer
    rejection_sampling(number_of_samples, rand_parameters, 0, 3.14);
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

    // Question 2.

    int64_t total_photon_count = 1000;

    photon_scattering(total_photon_count);

    return 0;
}
