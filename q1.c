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
                              double (*compute_value2)(double, double))
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
        lookup_array[2] = compute_value2(lookup_array[0], lookup_array[1]);
        fprintf(lookup_file, "%lf,%lf,%lf\n", lookup_array[0], lookup_array[1], lookup_array[2]);
    }
    else
    {
        fprintf(lookup_file, "%lf,%lf\n", lookup_array[0], lookup_array[1]);
    }

    printf("Generating lookup table and saving to %s...\n", filename);
    for (int32_t i = 1; i < length_of_array; i++)
    {
        // Generate the next z value
        lookup_array[i * columns] = lookup_array[(i - 1) * columns] + increment;
        // Compute the corresponding y value using the first function
        lookup_array[i * columns + 1] = compute_value1(lookup_array[i * columns]);

        // If a second function is provided, compute the corresponding z value
        if (columns == 3)
        {
            lookup_array[i * columns + 2] = compute_value2(lookup_array[i * columns], lookup_array[i * columns + 1]);
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

double get_z_from_t(double t)
{
    double r = 10.0;
    if (t < 0)
        t += 10;

    // Edge case where at the end of the circle
    else if (t = 20)
        return 10;

    else if (t > 0)
        t -= 10;

    return r * (2.0 / pi) * asin(t / r);
}

double get_y_from_t_and_z(double t, double z)
{
    double r = 10.0;

    if (t < 0)
        return sqrt(r * r - z * z);

    else if (t > 0)
        return -sqrt(r * r - z * z);
}

void photon_scattering_lookup(double t, double *lookup_table, int32_t size_of_table, double *dz, double *dy)
{
    int32_t low = 0, high = size_of_table - 1;

    // Binary search
    while (low <= high)
    {
        int32_t mid = low + (high - low) / 2;
        double mid_t = lookup_table[3 * mid];

        if (mid_t < t)
        {
            low = mid + 1;
        }
        else
        {
            high = mid - 1;
        }
    }

    // update dz and dy
    *dz = (lookup_table[3 * low + 1] + lookup_table[3 * high + 1]) * 0.5;
    *dy = (lookup_table[3 * low + 2] + lookup_table[3 * high + 2]) * 0.5;
}

void photon_scattering(int64_t total_photon_count, int64_t *random_paramteres)
{

    // Size of table , start value, final value
    int32_t size_of_table = 25000;
    int32_t lookup_table_params[3] = {size_of_table, -20, 20};
    double *displacement_lookup_array = generate_lookup_array(lookup_table_params, "displacement_direction_lookup.csv", *get_z_from_t, *get_y_from_t_and_z);

    FILE *photon_scattering;
    if ((photon_scattering = fopen("photon_scattering_results.csv", "w")) == NULL)
    {
        printf("Error opening file!\n");
    }

    for (int64_t i = 0; i < total_photon_count; i++)
    {
        // Initial values
        double y = 0, z = 0, dz = 0, dy = 0;

        int q = 0;

        // loop round until photon escapes
        while (z >= 0 && z <= 200)
        {
            double t = random_number_with_mapping(random_paramteres, -20, 20);

            // pass in the memory address of dz and dy to update them, If I was to generalise to N dimensions a array would be better but this works alright for now
            photon_scattering_lookup(t, displacement_lookup_array, size_of_table, &dz, &dy);

            z += dz;
            y += dy;
            q++;
        }

        // Check which direction the photon escaped form
        if (z < 0)
        {
            --i;
        }

        else
        {
            // As the material is infinte in the y plane, assume if it didnt escape in z < 0, it escaped in z >0.
            // Then from this calculate the angle from the origin

            // First, adjustment for overshooting. we want where it crosses z =200

            double overshoot_z = z - 200;
            // As the ratios will remain same we can find the overshoot y as a ratio releative to the z overshoot
            double overshoot_y = overshoot_z * dy / dz;

            // Apply the overshoot correction
            y -= overshoot_y;

            // Find the angle.
            double theta = atan(fabs(y) / 200);

            fprintf(photon_scattering, "%f,%f,%d\n", y, theta, q);
        }
    }
    fclose(photon_scattering);
}

int main()
{

    //
    // As seed, a,c,m will all be used in operations with each other I am declaring them as the same type to ensure type safe code.
    int64_t seed = 1, a = 16807, c = 0, m = 2147483647;
    int64_t rand_parameters[] = {seed, a, c, m};

    // Questinon 1.

    int32_t number_of_samples = 500000;

    clock_t start_rejection = clock();
    rejection_sampling(number_of_samples, rand_parameters, 0, 3.14);
    clock_t end_rejection = clock();
    double time_rejection = (double)(end_rejection - start_rejection) / CLOCKS_PER_SEC;

    printf("Time taken for rejection_sampling: %.6f seconds\n", time_rejection);

    clock_t start_direct = clock();
    direct_mapping(number_of_samples, rand_parameters);
    clock_t end_direct = clock();
    double time_direct = (double)(end_direct - start_direct) / CLOCKS_PER_SEC;
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

    int64_t total_photon_count = 1000000;

    clock_t start_photon_scattering = clock();
    photon_scattering(total_photon_count, rand_parameters);
    clock_t end_photon_scattering = clock();
    double time_photon_scattering = (double)(end_photon_scattering - start_photon_scattering) / CLOCKS_PER_SEC;

    printf("Photon scattering completed in %.6f\n", time_photon_scattering);

    return 0;
}
