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

double binary_search_2_columns(double y, const double *lookup_table, int32_t size_of_table)
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
        mu = binary_search_2_columns(y, y_to_mu_lookup_array, lookup_table_params[0]);

        fprintf(efficient_results, "%f,%f\n", mu, y);
    }

    fclose(efficient_results);
    free(y_to_mu_lookup_array);
}

double get_z_from_u(double u)
{
    // In theory if the distance travelled
    double r = 1.0;
    if (u < 0)
        u += 1;

    // Edge case where at the end of the circle ( 1.9999 as floating points)
    else if (u > 1.999999)
        return 1;

    else if (u > 0)
        u -= 1;

    return r * (2.0 / pi) * asin(u / r);
}

double get_y_from_u_and_z(double u, double z)
{
    double r = 1.0;

    if (u < 0)
        return sqrt(r * r - z * z);

    else if (u > 0)
        return -sqrt(r * r - z * z);
}

// While it is possible to merge the lookups into one, I ran into an issue where in C only a single value is returned.
// I think this function is distinct enough to warrent its existance. but I think the functionality of binary_search_2_columns could be rewrote to use pointers Or maybe a more general lookup taking X pointers in or something
// TODO: See comment above? maybe, not sure if its the most important thing in the world, and there is other piorities first.
void binary_search_3_columns(double u, double *lookup_table, int32_t size_of_table, double *dz, double *dy)
{
    int32_t low = 0, high = size_of_table - 1;

    // Binary search
    while (low <= high)
    {
        int32_t mid = low + (high - low) / 2;
        double mid_u = lookup_table[3 * mid];

        if (mid_u < u)
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

// Function to find bin index based on theta
int32_t find_bin_index(double theta, double *bin_edges, int32_t number_of_bins)
{
    for (int32_t i = 0; i < number_of_bins; i++)
    {
        if (theta >= bin_edges[i] && theta < bin_edges[i + 1])
        {
            return i;
        }
    }
}

double calculate_b_from_theta(double theta)
{
    // 2 theta * 1/2pi sin(4 pi theta). I dislike doing division in C, but this generation of the lookup table is not a bottle neck one bit so ah well.
    // Also its an odd function so it functions as intended whenever theta < 0
    return 2 * theta + (double)(1 / (2 * pi)) * sin(4 * pi * theta);
}
void photon_scattering(int64_t total_photon_count, int64_t *random_parameters, int do_rayleigh_scattering, double mean_free_path, const char *results_file_name, const char *binned_file_name)
{
    // Size of table, start value, final value
    int32_t size_of_table = 5000, number_of_bins = 10;
    int32_t lookup_table_params[3] = {size_of_table, -2, 2};

    double escape_z = 200;
    double *displacement_lookup_array = generate_lookup_array(lookup_table_params, "displacement_direction_lookup.csv", *get_z_from_u, *get_y_from_u_and_z);

    FILE *photon_scattering_results, *binned_intensity;
    if ((photon_scattering_results = fopen(results_file_name, "w")) == NULL ||
        (binned_intensity = fopen(binned_file_name, "w")) == NULL)
    {
        printf("Error opening file!\n");
        return;
    }

    // I should only do this if rayliegh scattering, ah well.
    int32_t size_of_scattering_table = 2000;
    int32_t scattering_lookup_table_params[3] = {size_of_scattering_table, -2, 2};
    double *scattering_lookup_array = generate_lookup_array(scattering_lookup_table_params, "scattering_angle_lookup.csv", *calculate_b_from_theta, NULL);

    int32_t bin_counts[number_of_bins];
    double bin_midpoints[number_of_bins], bin_intensity[number_of_bins], bin_edges[number_of_bins + 1];

    // Set up the bins and their mid points
    for (int32_t i = 0; i <= number_of_bins; i++)
    {
        bin_edges[i] = (-pi / 2.0) + pi * i / number_of_bins;
    }

    for (int32_t i = 0; i < number_of_bins; i++)
    {
        double theta_mid = (bin_edges[i] + bin_edges[i + 1]) / 2.0;
        bin_midpoints[i] = cos(theta_mid);
    }

    // There has to be a better way to do this, I do not like this but it fixes a bug.
    for (int32_t i = 0; i < number_of_bins; i++)
    {
        bin_counts[i] = 0;
    }

    for (int32_t i = 0; i < total_photon_count; i++)
    {
        // Initial values
        double y = 0, z = 0, dz = 0, dy = 0, t = 0, sum_t = 0, theta = 0, u = 0;

        int q = 0;

        // Start it towards the z direction if Rayleigh scattering
        if (do_rayleigh_scattering)
        {
            z += -log(random_number_with_mapping(random_parameters, 0, 1)) * mean_free_path;
            u += 2;
        }

        // Loop until photon escapes
        while (z >= 0 && z <= escape_z)
        {

            // Getting a bit too indented for my liking but it will have to do.
            if (do_rayleigh_scattering)
            {
                // Adjust theta by a probabilistic amount
                double b = random_number_with_mapping(random_parameters, -4, 4);
                theta = binary_search_2_columns(b, scattering_lookup_array, size_of_scattering_table);
                u += theta;

                // Check for a complete rotation
                if (u > 2)
                    u -= 4;
                if (u < -2)
                    u += 4;
            }
            else
            {
                u = random_number_with_mapping(random_parameters, -2, 2);
            }

            // Pass in the memory address of dz and dy to update them, If I was to generalise to N dimensions a array would be better but this works alright for now
            binary_search_3_columns(u, displacement_lookup_array, size_of_table, &dz, &dy);

            // Calculate the distance travelled using a direct map. TODO: another lookup table for this one?
            t = -log(random_number_with_mapping(random_parameters, 0, 1)) * mean_free_path;

            sum_t += t;
            z += dz * t;
            y += dy * t;
            q++;
        }

        // Compute and print average t
        if (q > 0)
        {
            double avg_t = sum_t / (double)q;
            // printf("Average t: %f\n", avg_t);
        }

        // Check which direction the photon escaped from
        if (z < 0)
        {
            --i;
        }
        else
        {
            // As the material is infinite in the y plane, assume if it didn't escape in z < 0, it escaped in z > 0.
            // Then from this calculate the angle from the origin

            // First, adjustment for overshooting. we want where it crosses z = 200
            double overshoot_z = z - escape_z;
            // As the ratios will remain the same we can find the overshoot y as a ratio relative to the z overshoot
            double overshoot_y = overshoot_z * dy / dz;

            // Apply the overshoot correction
            y -= overshoot_y;

            // Find the angle.
            double theta_lab = atan(y / escape_z);
            double intensity = cos(theta_lab);

            fprintf(photon_scattering_results, "%f,%f,%d\n", y, intensity, q);

            int32_t bin_index = find_bin_index(theta_lab, bin_edges, number_of_bins);
            bin_counts[bin_index]++;
        }
    }
    fclose(photon_scattering_results);

    // Write binned intensity results
    for (int32_t i = 0; i < number_of_bins; i++)
    {
        bin_intensity[i] = (double)bin_counts[i] * fabs(bin_midpoints[i]) / (double)total_photon_count;
        fprintf(binned_intensity, "%lf,%lf,%lf,%lf\n", bin_midpoints[i], bin_intensity[i], bin_edges[i], bin_edges[i + 1]);
    }
    fclose(binned_intensity);
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

    int64_t total_photon_count = 100000;

    clock_t start_photon_scattering = clock();
    photon_scattering(total_photon_count, rand_parameters, 0, 10, "photon.csv", "photon_bins.csv");
    clock_t end_photon_scattering = clock();
    double time_photon_scattering = (double)(end_photon_scattering - start_photon_scattering) / CLOCKS_PER_SEC;

    printf("Photon scattering completed in %.6f\n", time_photon_scattering);

    // Question 3

    photon_scattering(total_photon_count, rand_parameters, 1, 10, "rayleigh_blue.csv", "rayleigh_blue_bins.csv");
    photon_scattering(total_photon_count, rand_parameters, 1, 0.1, "rayleigh_other.csv", "rayleigh_other_bins.csv");

    return 0;
}
