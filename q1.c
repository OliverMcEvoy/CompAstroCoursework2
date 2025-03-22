#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <omp.h>
#include <string.h>

// While I dislike using global varibales it makes sense for pi, I try not to use MP_PI as a refrence for pi as it is not defined in standard C.
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
    // Extract values, not neccessary but I think it improves readility. when the code gets compiled down this is moot code.
    // Note I decided to use a pointer to the seed as if I implement processing over multiple threads, pointers are thread safe while static varibles are typically not
    int64_t *seed = &random_parameters[0];
    int64_t *a = &random_parameters[1];
    int64_t *c = &random_parameters[2];
    int64_t *m = &random_parameters[3];

    // This code calcuates and updates the seed.
    (*seed) = ((*a) * (*seed) + (*c)) % (*m);

    // Map to range.
    // If i was to hardcode 1 / m I would get a efficncy boost, but at the loss of adaptablity so I will not.
    double normalised = (double)(*seed) / (double)(*m);
    return (normalised * (upper_bound - lower_bound)) + lower_bound;
}
/**
 * Generate a lookup array based on a given function (or 2)
 * @param params an array consisting of Table size, start value, finish value
 * @param filename Name of the file, Note: remeber to include .csv in the file name
 * @param compute_value1 Compute the corrosponding value for the 2nd column
 * @param compute_value2 Compute the corrosponding value for the 3rd column, Set to NULL if not needed
 */
double *generate_lookup_array(int32_t params[3], const char *filename,
                              double (*compute_value1)(double),
                              double (*compute_value2)(double, double),
                              double (*compute_value3)(double, double))
{
    // Extract the parameters from the array
    int32_t length_of_array = params[0];
    double initial_value = params[1];
    double final_value = params[2];

    // Open the file in read mode to load a previously saved table
    FILE *lookup_file = fopen(filename, "r");

    double increment = (final_value - initial_value) / ((double)length_of_array - 1.0);

    // Determine the number of columns based on whether the second function is provided. could add a use case if a 2nd is provided and a 3rd is not, but not needed atm
    int32_t columns = (compute_value2 == NULL) ? 2 : 4;

    double *lookup_array = (double *)malloc(length_of_array * columns * sizeof(double));

    // If file already
    if (lookup_file != NULL)
    {
        for (int32_t i = 0; i < length_of_array; i++)
        {
            // Read the csv diffiently depending on the amount of columns
            int count = (columns == 2) ? fscanf(lookup_file, "%lf,%lf",
                                                &lookup_array[i * columns],
                                                &lookup_array[i * columns + 1])
                                       : fscanf(lookup_file, "%lf,%lf,%lf,%lf",
                                                &lookup_array[i * columns],
                                                &lookup_array[i * columns + 1],
                                                &lookup_array[i * columns + 2],
                                                &lookup_array[i * columns + 3]);
            if (count != columns)
            {
                printf("ERROR: fscanf matching error in %s on line %d\n", filename, i);
                fclose(lookup_file);
                remove(filename);
                free(lookup_array);
                return generate_lookup_array(params, filename, compute_value1, compute_value2, compute_value3);
            }
        }
        fclose(lookup_file);
        printf("Loading table %s\n", filename);
        return lookup_array;
    }

    // Otherwise, create a new file and generate the lookup table
    lookup_file = fopen(filename, "w");

    // First iteration handle seperately, due to initial value tomfoolery.
    lookup_array[0] = initial_value;
    lookup_array[1] = compute_value1(lookup_array[0]);
    if (columns == 4)
    {
        lookup_array[2] = compute_value2(lookup_array[0], lookup_array[1]);
        lookup_array[3] = compute_value3(lookup_array[0], lookup_array[1]);
        fprintf(lookup_file, "%lf,%lf,%lf,%lf\n", lookup_array[0], lookup_array[1], lookup_array[2], lookup_array[3]);
    }
    else
    {
        fprintf(lookup_file, "%lf,%lf\n", lookup_array[0], lookup_array[1]);
    }

    // printf("Lookup table saved to %s...\n", filename);
    for (int32_t i = 1; i < length_of_array; i++)
    {
        // Generate the value for the first column
        lookup_array[i * columns] = lookup_array[(i - 1) * columns] + increment;
        // Compute the corresponding y value using the first function
        lookup_array[i * columns + 1] = compute_value1(lookup_array[i * columns]);

        // If a second function is provided, compute the corresponding value in the third column
        if (columns == 4)
        {
            lookup_array[i * columns + 2] = compute_value2(lookup_array[i * columns], lookup_array[i * columns + 1]);
            lookup_array[i * columns + 3] = compute_value3(lookup_array[i * columns], lookup_array[i * columns + 1]);
            fprintf(lookup_file, "%lf,%lf,%lf,%lf\n", lookup_array[i * columns], lookup_array[i * columns + 1], lookup_array[i * columns + 2], lookup_array[i * columns + 3]);
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

// The inverse culumitive distribution function for y
double calculate_y_from_mu(double mu)
{
    return (pow(mu, 3) + 3 * mu + 4) * 0.125;
}

double calculate_z_from_u(double u)
{
    // Deal with square root of 1.
    if (u > 0.99999)
        return 1.0;
    return copysign(1.0 - sqrt(1.0 - u * u), u);
}

/**
 * Compute y using a deterministic azimuthal angle
 * @param u A value in range [-2, 2]
 * @param z The calculated z-coordinate
 * @return y-coordinate on the sphere
 */
double calculate_y_from_u_and_z(double u, double z)
{
    double phi = 50.0 * pi * u;   // Deterministic azimuthal angle
    double r = sqrt(1.0 - z * z); // Radius at that z slice
    return r * sin(phi);
}

/**
 * Compute x using a deterministic azimuthal angle
 * @param u A value in range [-2, 2]
 * @param z The calculated z-coordinate
 * @return x-coordinate on the sphere
 */
double calculate_x_from_u_and_z(double u, double z)
{
    double phi = 50.0 * pi * u;   // Deterministic azimuthal angle
    double r = sqrt(1.0 - z * z); // Radius at that z slice
    return r * cos(phi);
}
void binary_search_2_columns(double y, const double *lookup_table, int32_t size_of_table, double *mu)
{
    // set boundaries at the extremes of the table
    int32_t low = 0, high = size_of_table;

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
    *mu = (lookup_table[2 * high] + lookup_table[2 * low]) * 0.5;
}
// I think this function is distinct enough to warrent its existance. but I think the functionality of binary_search_2_columns could be rewrote to maybe be a more general lookup taking X pointers in or something
// TODO: See comment above? maybe, not sure if its the most important thing in the world, and there is other piorities first.
void binary_search_4_columns(double u, double *lookup_table, int32_t size_of_table, double *dz, double *dy, double *dx)
{
    int32_t low = 0, high = size_of_table;
    while (low <= high)
    {
        int32_t mid = low + (high - low) / 2;
        double mid_u = lookup_table[4 * mid];

        if (mid_u < u)
            low = mid + 1;
        else
            high = mid - 1;
    }

    // update the memory addresses directly
    *dz = (lookup_table[4 * low + 1] + lookup_table[4 * high + 1]) * 0.5;
    *dy = (lookup_table[4 * low + 2] + lookup_table[4 * high + 2]) * 0.5;
    *dx = (lookup_table[4 * low + 3] + lookup_table[4 * high + 3]) * 0.5;
}
// // See report for derivation. NOT NEEDED ANYMORE THANK GOD
// double calculate_b_from_theta(double theta)
// {
//     // 2 theta * 1/2pi sin(4 pi theta). I dislike doing division in C, but the generation of the lookup table is not a bottle neck one bit so ah well.
//     // Also its an odd function so it functions as intended whenever theta < 0
//     return 2 * theta + (double)(1 / (2 * pi)) * sin(4 * pi * theta);
// }
/**
 * Function to calculate the probability of mu according to the given formula.
 * @param mu The input for which we want to calculate the probability.
 * @return The calculated probability to be used wether to accept a given mu or not, see report for derivation
 */
double probability_of_mu(double mu)
{
    // Expressing 3/8 as 0.375 to avoid floating point division
    return 0.375 * (1 + pow(mu, 2));
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

// The logic for generating random numbers using direct mapping and a lookup table
void direct_mapping(int32_t number_of_samples, int64_t *random_parameters)
{

    int32_t lookup_table_params[3] = {2000, -1, 1};
    double *y_to_mu_lookup_array = generate_lookup_array(lookup_table_params, "y_to_mu_lookup.csv", *calculate_y_from_mu, NULL, NULL);

    double y, mu;
    FILE *efficient_results;
    efficient_results = fopen("efficient_results.csv", "w");

    for (int32_t i = 0; i < number_of_samples; i++)
    {

        y = random_number_with_mapping(random_parameters, 0, 1);

        // As this is a 1D array, the accessing of the array is simply the row index i
        binary_search_2_columns(y, y_to_mu_lookup_array, lookup_table_params[0], &mu);

        fprintf(efficient_results, "%f,%f\n", mu, y);
    }

    fclose(efficient_results);
    free(y_to_mu_lookup_array);
}

// Function to find bin index based on theta
int32_t find_bin_index(double theta_origin, double *bin_edges, int32_t number_of_bins)
{
    for (int32_t i = 0; i < number_of_bins; i++)
    {
        if (theta_origin >= bin_edges[i] && theta_origin < bin_edges[i + 1])
        {
            return i;
        }
    }
}

// Function to normalise a 3D vector
void normalise(double v[3])
{
    double inverse_magnitude = 1 / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    v[0] *= inverse_magnitude;
    v[1] *= inverse_magnitude;
    v[2] *= inverse_magnitude;
}

// Function to multiply a vector by a scalar
void multiply_vector(double v[3], double x)
{
    v[0] *= x;
    v[1] *= x;
    v[2] *= x;
}

// Compute the dot product of two 3D vectors
double dot_product(double v1[3], double v2[3])
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// Compute the cross product of two 3D vectors (result stored in result vector)
void cross_product(double v1[3], double v2[3], double result[3])
{
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

// Rotate a vector around an axis by a given angle
void rotate_vector(double v[3], double axis[3], double angle)
{
    normalise(axis);

    // Rodrigues' rotation formula
    double dot = dot_product(v, axis);
    double cross[3];
    cross_product(v, axis, cross);

    // Apply rotation
    for (int i = 0; i < 3; i++)
    {
        v[i] = v[i] * cos(angle) + cross[i] * sin(angle) + axis[i] * dot * (1 - cos(angle));
    }
}

// Does the photon scatter or not
int32_t photon_scatters(double absorption_probability, int64_t *random_parameters)
{
    double scatter = random_number_with_mapping(random_parameters, 0, 1);

    if (scatter > absorption_probability)
        return 1;
    // Else isn't actually needed but readability and stuff.
    else
        return 0;
}

void photon_scattering(int64_t total_photon_count,
                       int64_t *random_parameters,
                       int do_rayleigh_scattering,
                       double mean_free_path,
                       double absorption_probability,
                       const char *binned_angle_file_name,
                       const char *binned_position_file_name)
{
    // Size of table, start value, final value
    int32_t size_of_displacement_table = 25000, number_of_angle_bins = 10, number_of_position_bins = 4000;
    int32_t lookup_table_params[3] = {size_of_displacement_table, -1, 1};

    double escape_z = 200;
    double *displacement_lookup_array = generate_lookup_array(lookup_table_params, "displacement_direction_lookup.csv", *calculate_z_from_u, *calculate_y_from_u_and_z, *calculate_x_from_u_and_z);

    FILE *binned_angle_results, *binned_position_results;
    binned_angle_results = fopen(binned_angle_file_name, "w");
    binned_position_results = fopen(binned_position_file_name, "w");

    // I should only do this if Rayleigh scattering, ah well.
    int32_t size_of_scattering_table = 25000;
    int32_t scattering_lookup_table_params[3] = {size_of_scattering_table, -1, 1};
    double *scattering_lookup_array = generate_lookup_array(scattering_lookup_table_params, "rayleigh_scattering_lookup.csv", *calculate_y_from_mu, NULL, NULL);

    // Angle-based binning setup
    int32_t angle_bin_counts[number_of_angle_bins];
    double angle_bin_edges[number_of_angle_bins + 1], angle_bin_midpoints[number_of_angle_bins], angle_bin_intensity[number_of_angle_bins];

    double position_bin_edges[number_of_position_bins + 1];

    // Counts for x, y, z, and q
    int32_t x_bin_counts[number_of_position_bins];
    int32_t y_bin_counts[number_of_position_bins];
    int32_t z_bin_counts[number_of_position_bins];
    int32_t q_bin_counts[number_of_position_bins];
    // Set up the angle bins for equal solid angle
    for (int i = 0; i <= number_of_angle_bins; i++)
    {
        double fraction = (double)i / number_of_angle_bins;
        angle_bin_edges[i] = acos(1.0 - fraction);
    }

    // Initialize bin counts
    for (int i = 0; i < number_of_angle_bins; i++)
    {
        angle_bin_counts[i] = 0;
    }

    // Compute midpoints
    for (int i = 0; i < number_of_angle_bins; i++)
    {
        angle_bin_midpoints[i] = 0.5 * (angle_bin_edges[i] + angle_bin_edges[i + 1]);
    }

    // Set up the position bins and their edges
    double position_bin_width = 1.0; // Each bin represents an integer range (e.g., -2000 to -1999, -1999 to -1998, etc.)
    for (int32_t i = 0; i <= number_of_position_bins; i++)
    {
        position_bin_edges[i] = -2000 + i * position_bin_width;
    }

    for (int32_t i = 0; i < number_of_position_bins; i++)
    {
        x_bin_counts[i] = 0;
        y_bin_counts[i] = 0;
        z_bin_counts[i] = 0;
        q_bin_counts[i] = 0;
    }

#pragma omp parallel for
    for (int32_t i = 0; i < total_photon_count; i++)
    {
        // Initial values
        double position[3] = {0, 0, 0};  // [x, y, z]
        double direction[3] = {0, 0, 0}; // [dx, dy, dz]
        double t = 0, u = 0, theta = 0, theta_displacement = 0, theta_lab = 0;
        int32_t q = 0;
        double b_vector[3] = {0, 0, 0};

        // Start it towards the z direction if Rayleigh scattering, and also initialize 2 more vectors
        if (do_rayleigh_scattering)
        {
            position[2] += -log(random_number_with_mapping(random_parameters, 0, 1)) * mean_free_path;
            theta = 1;
            if (!photon_scatters(absorption_probability, random_parameters))
                continue;
        }

        // Loop until photon escapes
        while (position[2] >= 0 && position[2] <= escape_z)
        {
            if (do_rayleigh_scattering)
            {
                if (!photon_scatters(absorption_probability, random_parameters))
                    continue;
                // Find what u is meant to be
                double b = random_number_with_mapping(random_parameters, 0, 1);
                binary_search_2_columns(b, scattering_lookup_array, size_of_scattering_table, &theta_displacement);

                b_vector[0] = position[0];
                b_vector[1] = position[1];
                b_vector[2] = position[2];
                normalise(b_vector);
                multiply_vector(b_vector, theta_displacement);

                // Calculate the angle between b_vector and the z-axis.
                double z_axis[3] = {0, 0, 1};
                double dot = dot_product(b_vector, z_axis);
                double angle = acos(dot);

                // Calculate the rotation axis (cross product of b_vector and z-axis).
                double rotation_axis[3];
                cross_product(b_vector, z_axis, rotation_axis);
                normalise(rotation_axis);

                // Rotate the direction vector
                rotate_vector(direction, rotation_axis, angle);

                // Technically I could always work in terms of u but I think this is more readable.
                u = theta_displacement;

                // Pass in the memory address of direction array to update dx, dy, dz
                binary_search_4_columns(u, displacement_lookup_array, size_of_displacement_table, &direction[0], &direction[1], &direction[2]);

                // Calculate the distance travelled using a direct map
                t = -log(random_number_with_mapping(random_parameters, 0, 1)) * mean_free_path;

                position[0] += direction[0] * t;
                position[1] += direction[1] * t;
                position[2] += direction[2] * t;
                q++;
            }
            else
            {
                u = random_number_with_mapping(random_parameters, -1, 1);

                // Pass in the memory address of direction array to update dx, dy, dz
                binary_search_4_columns(u, displacement_lookup_array, size_of_displacement_table, &direction[2], &direction[1], &direction[0]);

                // Calculate the distance travelled using a direct map
                t = -log(random_number_with_mapping(random_parameters, 0, 1)) * mean_free_path;

                if (!photon_scatters(absorption_probability, random_parameters))
                    continue;

                position[2] += direction[2] * t;
                position[1] += direction[1] * t;
                position[0] += direction[0] * t;
                q++;
            }
        }

        if (position[2] < 0)
        {
            --i;
            continue;
        }

        int32_t x_bin_index = (int32_t)(position[0] + 2000);
        int32_t y_bin_index = (int32_t)(position[1] + 2000);
        int32_t z_bin_index = (int32_t)(position[2] + 2000);
        theta_lab = atan(position[1] / escape_z);

        // #pragma omp critical
        {
            // Angle-based binning
            int32_t angle_bin_index = find_bin_index(theta_lab, angle_bin_edges, number_of_angle_bins);
            angle_bin_counts[angle_bin_index]++;

            // Position-based binning.
            if (x_bin_index >= 0 && x_bin_index < number_of_position_bins)
            {
                x_bin_counts[x_bin_index]++;
            }
            if (y_bin_index >= 0 && y_bin_index < number_of_position_bins)
            {
                y_bin_counts[y_bin_index]++;
            }
            if (z_bin_index >= 0 && z_bin_index < number_of_position_bins)
            {
                z_bin_counts[z_bin_index]++;
            }

            // Use z_bin_index for q counts (or adjust as needed)
            if (z_bin_index >= 0 && z_bin_index < number_of_position_bins)
            {
                q_bin_counts[z_bin_index] += q;
            }
        }
    }

    for (int32_t i = 0; i < number_of_angle_bins; i++)
    {
        angle_bin_intensity[i] = (double)angle_bin_counts[i] * cos(angle_bin_midpoints[i]) / (double)total_photon_count;
        fprintf(binned_angle_results, "%lf,%lf,%lf,%lf\n", angle_bin_midpoints[i], angle_bin_intensity[i], angle_bin_edges[i], angle_bin_edges[i + 1]);
    }
    fclose(binned_angle_results);

    // Write position-based binning results
    for (int32_t i = 0; i < number_of_position_bins; i++)
    {
        fprintf(binned_position_results, "%lf,%lf,%d,%d,%d,%d\n", position_bin_edges[i], position_bin_edges[i + 1], x_bin_counts[i], y_bin_counts[i], z_bin_counts[i], q_bin_counts[i]);
    }
    fclose(binned_position_results);

    free(displacement_lookup_array);
    free(scattering_lookup_array);
}
// Time Macro as the code was getting a bit messy doing this at the start and end of every function.
#define MEASURE_TIME(func, ...) ({                                                                  \
    clock_t start = clock();                                                                        \
    func(__VA_ARGS__);                                                                              \
    clock_t end = clock();                                                                          \
    double time_taken = (double)(end - start) / CLOCKS_PER_SEC;                                     \
    printf("Total compute time across all threads used for %s: %.6f seconds\n", #func, time_taken); \
    time_taken;                                                                                     \
})

int main()
{
    int64_t seed = 1, a = 16807, c = 0, m = 2147483647;
    int64_t rand_parameters[] = {seed, a, c, m};

    // Question 1
    int32_t number_of_random_samples = 500000;
    double time_rejection = MEASURE_TIME(rejection_sampling, number_of_random_samples, rand_parameters, 0, 3.14);
    double time_direct = MEASURE_TIME(direct_mapping, number_of_random_samples, rand_parameters);

    if (time_rejection < time_direct)
    {
        printf("rejection_sampling was faster by %.6f seconds\n", time_direct - time_rejection);
    }
    else
    {
        printf("direct_mapping was faster by %.6f seconds\n", time_rejection - time_direct);
    }

    int32_t total_photon_count = 1000000;

    // Question 2
    MEASURE_TIME(photon_scattering, total_photon_count, rand_parameters, 0, 200 / 10, 0, "photon_bins.csv", "photon_position_bins.csv");

    // Question 3
    MEASURE_TIME(photon_scattering, total_photon_count, rand_parameters, 1, 200 / 10, 0, "rayleigh_blue_bins.csv", "rayleigh_blue_position_bins.csv");

    MEASURE_TIME(photon_scattering, total_photon_count, rand_parameters, 1, 200 / 0.1, 0, "rayleigh_other_bins.csv", "rayleigh_other_position_bins.csv");

    return 0;
}