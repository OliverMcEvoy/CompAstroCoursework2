#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <omp.h>

// While I dislike using global varibales it makes sense for pi, I try not to use MP_PI as a refrence for pi as it is not defined in standard C.
double pi = 3.141592653589793;

//// GENERATE MAPPED RANDOMN NUMBER AND GENERATE LOOKUP TABLE HELPER FUNCTIONS ////

/**
 * Get a random number within a range.
 *
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
 * Generate a lookup array based on a given function (or 3!).
 *
 * @param params an array consisting of Table size, start value, finish value.
 * @param filename Name of the file, Note: remeber to include .csv in the file name.
 * @param compute_value1 Compute the corrosponding value for the 2nd column.
 * @param compute_value2 Compute the corrosponding value for the 3rd column, Set to NULL if not needed.
 */
double *generate_lookup_array(int32_t params[3], const char *filename,
                              double (*compute_value1)(double),
                              double (*compute_value2)(double, double),
                              double (*compute_value3)(double, double))
{
    // Extract the parameters from the array.
    int32_t length_of_array = params[0];
    double initial_value = (double)params[1];
    double final_value = (double)params[2];

    // Open in read mode to try load a previously saved table, if it is present.
    FILE *lookup_file = fopen(filename, "r");

    double increment = (final_value - initial_value) / ((double)length_of_array - 1.0);

    // Determine the number of columns based on whether the second function is provided. could add a use case if a 2nd is provided and a 3rd is not, but not needed atm
    int32_t columns = (compute_value2 == NULL) ? 2 : 4;

    double *lookup_array = (double *)malloc(length_of_array * columns * sizeof(double));

    // If file already present read from it.
    if (lookup_file != NULL)
    {
        for (int32_t i = 0; i < length_of_array; i++)
        {
            // Read the csv diffiently depending on the amount of columns.
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
                // Error sanity check. probably could include more error checking.
                printf("ERROR: fscanf matching error in %s on line %d\n", filename, i);
                fclose(lookup_file);
                remove(filename);
                free(lookup_array);
                return generate_lookup_array(params, filename, compute_value1, compute_value2, compute_value3);
            }
        }
        fclose(lookup_file);

        // typically if each question was isolated I would have this printed, but with the amount of lookup tables it didn't seem practical.
        // printf("Loading table %s\n", filename);
        return lookup_array;
    }

    // No file present ( or issues with current), make a new one.
    lookup_file = fopen(filename, "w");

    // First iteration handle seperately,
    lookup_array[0] = initial_value;
    lookup_array[1] = compute_value1(lookup_array[0]);

    // TODO Probably should make this work for N columns and pass in an array containing pointers to functions.
    // Print to the file the first iteration.
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

    // Generate values and print out the rest of the array.
    for (int32_t i = 1; i < length_of_array; i++)
    {
        // Generate the value for the first column.
        lookup_array[i * columns] = lookup_array[(i - 1) * columns] + increment;
        // Compute the corresponding y value using the first function.
        lookup_array[i * columns + 1] = compute_value1(lookup_array[i * columns]);

        if (columns == 4)
        {
            // If a second function is provided, compute the corresponding value in the third and fourth column usisng the provided function then print to file
            lookup_array[i * columns + 2] = compute_value2(lookup_array[i * columns], lookup_array[i * columns + 1]);
            lookup_array[i * columns + 3] = compute_value3(lookup_array[i * columns], lookup_array[i * columns + 1]);
            fprintf(lookup_file, "%lf,%lf,%lf,%lf\n", lookup_array[i * columns], lookup_array[i * columns + 1], lookup_array[i * columns + 2], lookup_array[i * columns + 3]);
        }
        else
        {
            // Otherwise just print to file.
            fprintf(lookup_file, "%lf,%lf\n", lookup_array[i * columns], lookup_array[i * columns + 1]);
        }
    }

    // typically if each question was isolated I would have this printed, but with the amount of lookup tables throughout this programme it didn't seem practical.
    // printf("Done generating lookup table, saved to %s\n", filename);
    fclose(lookup_file);

    return lookup_array;
}

//// FUNCTIONS THAT WILL BE PASSED INTO LOOKUP TABLE ////

/**
 * The inverse culumitive distribution function for y.
 */
double calculate_y_from_mu(double mu)
{
    return (pow(mu, 3) + 3 * mu + 4) * 0.125;
}

// Get z from a ranom variable u, see report for derivation.
double calculate_z_from_u(double u)
{
    return 2 * u - 1;
}

/**
 * Compute y from u and z.
 */
double calculate_y_from_u_and_z(double u, double z)
{
    double phi = 1500.0 * pi * u, r;
    // Deal with z = 1 causing nans.
    double r = (z > 0.99999) ? 0.0 : sqrt(1.0 - z * z);
    return r * sin(phi);
}

/**
 * Compute x from z and u.
 */
double calculate_x_from_u_and_z(double u, double z)
{
    double phi = 1500.0 * pi * u, r;
    // Deal with z = 1 causing nans.
    double r = (z > 0.99999) ? 0.0 : sqrt(1.0 - z * z);
    return r * cos(phi);
}

//// THE DIFFERENT BINARY SEARCHES USED ////

void binary_search_2_columns(double y, const double *lookup_table, int32_t size_of_table, double *mu)
{
    // set boundaries at the extremes of the table.
    int32_t low = 0, high = size_of_table - 1;

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

    // Take the midpoint of the two closest mu values.
    *mu = (lookup_table[2 * high] + lookup_table[2 * low]) * 0.5;
}

// I think this function is distinct enough to warrent its existance. but I think the functionality of binary_search_2_columns could be rewrote to maybe be a more general lookup taking X pointers in or something
// TODO: See comment above? maybe, not sure if its the most important thing in the world, and there is other piorities first.
void binary_search_4_columns(double u, double *lookup_table, int32_t size_of_table, double *dz, double *dy, double *dx)
{
    int32_t low = 0, high = size_of_table - 1;
    while (low <= high)
    {
        int32_t mid = low + (high - low) / 2;
        double mid_u = lookup_table[4 * mid];

        if (mid_u < u)
            low = mid + 1;
        else
            high = mid - 1;
    }

    // update the memory addresses directly.
    *dz = (lookup_table[4 * low + 1] + lookup_table[4 * high + 1]) * 0.5;
    *dy = (lookup_table[4 * low + 2] + lookup_table[4 * high + 2]) * 0.5;
    *dx = (lookup_table[4 * low + 3] + lookup_table[4 * high + 3]) * 0.5;
}

/**
 * As z varies from -1 to 1 we need to make modifications to the binary search.
 */
void binary_search_b(double b, double *lookup_table, int32_t size_of_table, double *dy, double *dx)
{
    int32_t low, high, mid;

    // Determine the search range and order based on the sign of b.
    if (b >= 0)
    {
        // Search the upper half of the table (positive b) in ascending order.
        low = size_of_table / 2 + 1;
        high = size_of_table - 1;
    }
    else
    {
        // Search the lower half of the table (negative b) in descending order.
        low = 0;
        high = size_of_table / 2 - 1;
    }

    // Binary search depending on where b is.
    while (low <= high)
    {
        mid = low + (high - low) / 2;
        double mid_b = lookup_table[4 * mid + 1];

        if (mid_b < b)
            low = mid + 1;
        else
            high = mid - 1;
    }

    // Update dy and dx based on the closest match
    *dy = (lookup_table[4 * low + 2] + lookup_table[4 * high + 2]) * 0.5;
    *dx = (lookup_table[4 * low + 3] + lookup_table[4 * high + 3]) * 0.5;
}

/**
 * Function to calculate the probability of mu according to the given formula.
 * @param mu The input for which we want to calculate the probability.
 * @return The calculated probability to be used wether to accept a given mu or not, see report for derivation
 */
double probability_of_mu(double mu)
{
    // Expressing 3/8 as 0.375 to avoid floating point division.
    return 0.375 * (1 + pow(mu, 2));
}
/**
 * Function to perform rejection sampling according to the given parameters.
 *
 * @param num_samples The number of samples to generate.
 * @param seed Initial seed for random number generator.
 * @param a, c, m Parameters for the random number generator.
 */
void rejection_sampling(int32_t num_samples, int64_t *random_parameters)
{
    // Create 100 bins for mu values between -1 and 1.
    int32_t num_bins = 100;

    // Calloc used as technically a non zero chance a value will not fall into a bin.
    int32_t *bins = calloc(num_bins, sizeof(int32_t));
    double min_mu = -1.0;
    double max_mu = 1.0;
    double bin_width = (max_mu - min_mu) / num_bins;

    FILE *rejection_method_results = fopen("rejection_sampling_binned_results.csv", "w");

    // Increment until the total samples are reached.
    int32_t accepted_count = 0;
    while (accepted_count < num_samples)
    {
        // Get mu and other values.
        double mu = random_number_with_mapping(random_parameters, -1, 1);
        double prob_mu = probability_of_mu(mu);
        double acceptance_threshold = random_number_with_mapping(random_parameters, 0, 1);

        // logic if the mu is accepted.
        if (acceptance_threshold <= prob_mu)
        {
            // Bin the accepted sample.
            int32_t bin_index = (int32_t)((mu - min_mu) / bin_width);

            bins[bin_index]++;
            accepted_count++;
        }
    }

    // Write bin results to file.
    for (int32_t i = 0; i < num_bins; i++)
    {
        double bin_start = min_mu + i * bin_width;
        double bin_end = bin_start + bin_width;
        double normalised_value = (double)bins[i] / num_samples * (double)num_bins * 0.5;
        fprintf(rejection_method_results, "%f,%f,%f\n", bin_start, bin_end, normalised_value);
    }

    fclose(rejection_method_results);
    free(bins);
}

/**
 * Generate a random distribtion via direct mapping.
 *
 * @param number_of_samples This is the amount of samples taken.
 * @param random_parameters A pointer to the random parameters used for the random number generator.
 * @param output_filename The name of the file to save in.
 */
void direct_mapping(int32_t number_of_samples, int64_t *random_parameters, const char *output_filename)
{
    int32_t lookup_table_params[3] = {201, -1, 1};
    double *y_to_mu_lookup_array = generate_lookup_array(lookup_table_params, "y_to_mu_lookup.csv", *calculate_y_from_mu, NULL, NULL);

    // Create 100 bins for mu values between -1 and 1.
    int32_t num_bins = 100;
    // Theres 'technically' no gaurentee each bin will get filled. esp if a low number of samples are used, so Ill use calloc.
    int32_t *bins = calloc(num_bins, sizeof(int32_t));

    // Constants for binning.
    const double min_mu = -1.0;
    const double max_mu = 1.0;
    const double mu_range = max_mu - min_mu;
    const double bin_width = mu_range / num_bins;

    double y, mu;
    for (int32_t i = 0; i < number_of_samples; i++)
    {
        y = random_number_with_mapping(random_parameters, 0, 1);
        binary_search_2_columns(y, y_to_mu_lookup_array, lookup_table_params[0], &mu);

        // Determine which bin this mu value falls into.
        int32_t bin_index = (int32_t)((mu - min_mu) / bin_width);

        // And update the index of the bin.
        bins[bin_index]++;
    }

    FILE *efficient_results;
    efficient_results = fopen(output_filename, "w");

    // Write bin results to file.
    for (int32_t i = 0; i < num_bins; i++)
    {
        double bin_start = min_mu + i * bin_width;
        double bin_end = min_mu + (i + 1) * bin_width;
        double normalised_value = (double)bins[i] / number_of_samples * (double)num_bins * 0.5;

        fprintf(efficient_results, "%f,%f,%f\n", bin_start, bin_end, normalised_value);
    }

    fclose(efficient_results);
    free(y_to_mu_lookup_array);
    free(bins);
}

/**
 * Generate the expected probablity distribution of mu and save it as bins.
 */
void generate_expected_density_bins_for_plot_for_report()
{
    const int num_bins = 100;
    const double min_mu = -1.0;
    const double max_mu = 1.0;
    const double bin_width = (max_mu - min_mu) / num_bins;

    FILE *file = fopen("expected_density_bins.csv", "w");
    if (!file)
    {
        perror("Failed to create output file");
        return;
    }

    // Not sure if its counts as the Python plotting programme having to deal with the column names as a pitfall of the C programme so decided not to include them.
    // All the plotting programmes do is plot.
    // fprintf(file, "bin_start,bin_end,expected_density\n");

    for (int i = 0; i < num_bins; i++)
    {
        double bin_start = min_mu + i * bin_width;
        double bin_end = bin_start + bin_width;
        double bin_center = (bin_start + bin_end) / 2.0;
        double density = probability_of_mu(bin_center);

        fprintf(file, "%f,%f,%f\n", bin_start, bin_end, density);
    }

    fclose(file);
}

/**
 *   Function to find bin index based on theta between the scattered photon and the origin.
 *
 *  @param theta_origin Angle between the y value and the z axis.
 *  @param bin_edges where each bin ends and starts.
 *  @param number_of_bins How many total bins is there.
 */
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

/**
 * Normalise a vector.
 *
 * @param v vector
 */
void normalise_vector(double *v)
{
    // Done this way to reduce floating point division.
    double inverse_magnitude = 1 / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    v[0] *= inverse_magnitude;
    v[1] *= inverse_magnitude;
    v[2] *= inverse_magnitude;
}

/**
 * Multiple each of the values of a vector by a scalar.
 *
 * @param v Vector.
 * @param x Scalar.
 */
void multiply_vector(double *v, double x)
{
    v[0] *= x;
    v[1] *= x;
    v[2] *= x;
}

/**
 * Compute the dot product of two vectors.
 *
 * @param v1 Vector 1.
 * @param v2 Vector 2.
 */
double dot_product(double v1[3], double v2[3])
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

/**
 * Compute the cross product of two 3D vectors (result stored in result vector)
 *
 * @param v1 Vector 1.
 * @param v2 Vector 2.
 * @param result The vector resulting as a cross product between vector 1 and 2.
 */
void cross_product(double v1[3], double v2[3], double *result)
{

    // printf("cross product\n");
    // printf("v vectors %lf,%lf,%lf\n", v1[0], v1[1], v1[2]);
    // printf("v vectors %lf,%lf,%lf\n", v2[0], v2[1], v2[2]);
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

/**
 * Rotate a vector around an axis by a given angle using Rodrigues' rotation formula.
 *
 * @param v A vector to be rotated.
 * @param axis The axis on which to rotate the vector around.
 * @param angle The angle between the current position of the vector and the new desired positon.
 */
void rotate_vector(double *v, double axis[3], double angle)
{
    normalise_vector(axis);

    // Compute values needed for the rotation formula.
    double cos_theta = cos(angle);
    double sin_theta = sin(angle);
    double one_minus_cos_theta = 1 - cos_theta;
    double dot = dot_product(v, axis);

    // Compute the cross product.
    double cross[3];
    cross_product(v, axis, cross);

    // Apply Rodrigues' rotation formula onto the vector v
    v[0] = v[0] * cos_theta + cross[0] * sin_theta + axis[0] * dot * one_minus_cos_theta;
    v[1] = v[1] * cos_theta + cross[1] * sin_theta + axis[1] * dot * one_minus_cos_theta;
    v[2] = v[2] * cos_theta + cross[2] * sin_theta + axis[2] * dot * one_minus_cos_theta;
}

/**
 * Does the photon scatter or not
 *
 * @param absorption_probability The chance the photon will absorb.
 * @param random_parameters The random parameters for the random number generator.
 */
int32_t does_photon_scatter(double absorption_probability, int64_t *random_parameters)
{
    double scatter = random_number_with_mapping(random_parameters, 0, 1);

    if (scatter > absorption_probability)
        return 1;
    // Else isn't actually needed but readability and stuff.
    else
        return 0;
}

/**
 * Function to handle question 2 and 3 of the report.
 *
 * @param total_photon_count The total amount of photons to simmulate.
 * @param random_parameters A pointer to the inital parameters of the random number generator.
 * @param do_rayleigh_scattering A toggle to simulate rayleigh scattering or not, 0 to not 1 to do.
 * @param mean_free_path The mean free path of an photon.
 * @param absorption_probability The probability that a photon will be absorbed.
 * @param binned_angle_file_name The name of the file to save the distribution of the resulting angles.
 * @param binned_position_file_name The name of the file to save the distribution of the resulting positions.
 */
void photon_scattering(int32_t total_photon_count,
                       int64_t *random_parameters,
                       int do_rayleigh_scattering,
                       double mean_free_path,
                       double absorption_probability,
                       const char *binned_angle_file_name,
                       const char *binned_position_file_name)
{
    // Declare some inital values.
    double escape_z = 200;
    int32_t size_of_displacement_table = 50000, size_of_scattering_table = 15000, number_of_angle_bins = 10, number_of_position_bins = 4000;
    // Size of table, start value, final value.
    int32_t lookup_table_params[3] = {size_of_displacement_table, 0, 1};

    // Generate the displacement lookup array.
    double *displacement_lookup_array = generate_lookup_array(lookup_table_params, "displacement_direction_lookup.csv", *calculate_z_from_u, *calculate_y_from_u_and_z, *calculate_x_from_u_and_z);

    FILE *binned_angle_results, *binned_position_results;
    binned_angle_results = fopen(binned_angle_file_name, "w");
    binned_position_results = fopen(binned_position_file_name, "w");

    // Size of table, start value, final value.
    int32_t scattering_lookup_table_params[3] = {size_of_scattering_table, -1, 1};
    // I should only do this if Rayleigh scattering, but some compilers does not like that a varibale could be refrenced that has not been delcared (even though it wont actually be accesed).
    double *scattering_lookup_array = generate_lookup_array(scattering_lookup_table_params, "rayleigh_scattering_lookup.csv", *calculate_y_from_mu, NULL, NULL);

    // Note, whenever there is a possibility of a value being read later in the code, but it is not explicility set ( a non zero chance it could be read while not being defined) I will use calloc to ensure the array is initialised with values of 0.
    // calloc is not needed when everything is being set immediately.

    // Angle-based binning.
    int32_t *angle_bin_counts = calloc(number_of_angle_bins, sizeof(int32_t));
    double *angle_bin_edges = malloc((number_of_angle_bins + 1) * sizeof(double));
    double *angle_bin_midpoints = malloc(number_of_angle_bins * sizeof(double));
    double *angle_bin_intensity = calloc(number_of_angle_bins, sizeof(double));

    // Position-based binning.
    double *position_bin_edges = malloc((number_of_position_bins + 1) * sizeof(double));
    int32_t *x_bin_counts = calloc(number_of_position_bins, sizeof(int32_t));
    int32_t *y_bin_counts = calloc(number_of_position_bins, sizeof(int32_t));
    int32_t *z_bin_counts = calloc(number_of_position_bins, sizeof(int32_t));
    int32_t *q_bin_counts = calloc(number_of_position_bins, sizeof(int32_t));

    // Set up the angle bins for equal solid angle.
    for (int i = 0; i <= number_of_angle_bins; i++)
    {
        angle_bin_edges[i] = asin((double)i / number_of_angle_bins);
    }

    // Compute midpoints
    for (int32_t i = 0; i < number_of_angle_bins; i++)
    {
        angle_bin_midpoints[i] = 0.5 * (angle_bin_edges[i] + angle_bin_edges[i + 1]);
    }

    // Set up the position bins and their edges.
    double position_bin_width = 1.0;
    for (int32_t i = 0; i <= number_of_position_bins; i++)
    {
        position_bin_edges[i] = -2000 + i * position_bin_width;
    }

    // Variable to count total scattering events
    int32_t total_scattering_events = 0;

#pragma omp parallel for
    for (int32_t i = 0; i < total_photon_count; i++)
    {
    // A tag for the goto statement to refrence.
    generate_photon:
        // Initial values.
        double position[3] = {0, 0, 0};  // [x, y, z]
        double direction[3] = {0, 0, 0}; // [dx, dy, dz]
        double t = 0, u = 0, b = 0, theta = 0;
        int32_t scattering_count = 0;
        double b_vector[3] = {0, 0, 0};

        // Start it towards the z direction if Rayleigh scattering.
        if (do_rayleigh_scattering == 1)
        {
            // Calculate the distance travelled using a direct map.
            position[2] += -log(random_number_with_mapping(random_parameters, 0, 1)) * mean_free_path;

            // If photon gets absorpbed continue to the next iteration without decreasing the increment.
            if (!does_photon_scatter(absorption_probability, random_parameters))
                continue;
        }
        else
            // If its regular photon scattering as its isotropic we dont want to count the first random direction as a scattering event.
            --scattering_count;

        // Loop until photon escapes.
        while (position[2] >= 0 && position[2] <= escape_z)
        {
            // If photon gets absorpbed continue to the next iteration without decreasing the increment.
            if (!does_photon_scatter(absorption_probability, random_parameters))
                continue;

            if (do_rayleigh_scattering == 1)
            {
                // Find c, see report for more info on the mapping.
                double c = random_number_with_mapping(random_parameters, 0, 1);
                binary_search_2_columns(c, scattering_lookup_array, size_of_scattering_table, &b);

                // Set up the rotation axis, it is found by taking a cross product of the current photon position and the z_axis.
                double z_axis[3] = {0.0, 0.0, 1.0};
                double rotation_axis[3] = {0, 0, 0};
                cross_product(position, z_axis, rotation_axis);

                binary_search_b(b, displacement_lookup_array, size_of_displacement_table, &direction[1], &direction[0]);
                direction[2] = b;

                // Only rotate if rotation axis is non-zero vector, if it is zero, then it is already parallel and it is not needed to rotate.
                if (rotation_axis[0] != 0 || rotation_axis[1] != 0 || rotation_axis[2] != 0)
                {

                    // Calculate theta as the angle between z-axis and normalised position vector.
                    double position_normalised[3] = {position[0], position[1], position[2]};
                    normalise_vector(position_normalised);
                    double theta = acos(dot_product(position_normalised, z_axis));
                    normalise_vector(rotation_axis);

                    // Update the direction vector using Rodrigues' rotation formula.
                    rotate_vector(direction, rotation_axis, theta);
                }
                else if (scattering_count > 2)
                    printf("this shouldnt happen");

                // Calculate the distance travelled using a direct map.
                t = -log(random_number_with_mapping(random_parameters, 0, 1)) * mean_free_path;
                // Update the position based on the direction.
                position[0] += direction[0] * t;
                position[1] += direction[1] * t;
                position[2] += direction[2] * t;
                scattering_count++;
            }
            else
            {
                // Not rayleigh scattering.

                // Calculate the distance travelled using a direct map.
                t = -log(random_number_with_mapping(random_parameters, 0, 1)) * mean_free_path;
                u = random_number_with_mapping(random_parameters, 0, 1);

                // Pass in the memory address of direction array to update dx, dy, dz.
                binary_search_4_columns(u, displacement_lookup_array, size_of_displacement_table, &direction[2], &direction[1], &direction[0]);

                position[2] += direction[2] * t;
                position[1] += direction[1] * t;
                position[0] += direction[0] * t;
                scattering_count++;
            }
        }

        if (position[2] < 0)
            // While I could just decrement the index and use continue, this is more thread safe.
            goto generate_photon;

// Add to the total scattering events
#pragma omp atomic update
        total_scattering_events += scattering_count;

        // int manipulation so division isn't too bad.
        int32_t x_bin_index = (int32_t)(position[0] + 2000);
        int32_t y_bin_index = (int32_t)(position[1] + 2000);
        int32_t z_bin_index = (int32_t)(position[2] + 2000);

        // Angle-based binning, fabs to take the modules as it is only the angle we are intrested in.
        theta = atan(fabs(position[1] / escape_z));
        int32_t angle_bin_index = find_bin_index(theta, angle_bin_edges, number_of_angle_bins);

// Use atomic to update the bin counts to ensure some thread safe stuff.
#pragma omp atomic update
        angle_bin_counts[angle_bin_index]++;

        // Position-based binning using atomic increment.
        if (x_bin_index >= 0 && x_bin_index < number_of_position_bins)
        {
#pragma omp atomic update
            x_bin_counts[x_bin_index]++;
        }
        if (y_bin_index >= 0 && y_bin_index < number_of_position_bins)
        {
#pragma omp atomic update
            y_bin_counts[y_bin_index]++;
        }
        if (z_bin_index >= 0 && z_bin_index < number_of_position_bins)
        {
#pragma omp atomic update
            z_bin_counts[z_bin_index]++;
        }
    }

    // Calculate and print the average q.
    double average_q = (double)total_scattering_events / total_photon_count;
    printf("Average count of scattering events for an escaping photon: %lf\n", average_q);

    // Print the angle based binning to a csv.
    for (int32_t i = 0; i < number_of_angle_bins; i++)
    {
        // Calculate the intensity using the midpoint and mu = cos theta as requested by the coursework problem/
        angle_bin_intensity[i] = (double)angle_bin_counts[i] * cos(angle_bin_midpoints[i]) / (double)total_photon_count;
        fprintf(binned_angle_results, "%lf,%lf,%lf,%lf\n", angle_bin_midpoints[i], angle_bin_intensity[i], angle_bin_edges[i], angle_bin_edges[i + 1]);
    }
    fclose(binned_angle_results);

    // Print the position based binning to a csv.
    for (int32_t i = 0; i < number_of_position_bins; i++)
    {
        // Normalise by diving by total photon count..
        double x_bin_normalised = (double)x_bin_counts[i] / (double)total_photon_count;
        double y_bin_normalised = (double)y_bin_counts[i] / (double)total_photon_count;
        double z_bin_normalised = (double)z_bin_counts[i] / (double)total_photon_count;
        double q_bin_normalised = (double)q_bin_counts[i] / (double)total_photon_count;

        fprintf(binned_position_results, "%lf,%lf,%lf,%lf,%lf,%lf\n",
                position_bin_edges[i], position_bin_edges[i + 1],
                x_bin_normalised, y_bin_normalised, z_bin_normalised, q_bin_normalised);
    }
    fclose(binned_position_results);

    // Free everything.
    free(displacement_lookup_array);
    free(scattering_lookup_array);
    free(angle_bin_counts);
    free(angle_bin_edges);
    free(angle_bin_midpoints);
    free(position_bin_edges);
    free(x_bin_counts);
    free(y_bin_counts);
    free(z_bin_counts);
    free(q_bin_counts);
}

// Time Macro as the code was getting a bit messy doing this at the start and end of every function.
// When only 1 thread is passed in thread will be printed instead of threads.
#define MEASURE_TIME(func, num_threads, ...) ({                                                         \
    clock_t start = clock();                                                                            \
    func(__VA_ARGS__);                                                                                  \
    clock_t end = clock();                                                                              \
    double thread_time = (double)(end - start) / CLOCKS_PER_SEC;                                        \
    double total_time = thread_time / num_threads;                                                      \
    printf("Total compute time for %s: %.6f seconds (%.6f seconds total compute time, %.0f %s used)\n", \
           #func, total_time, thread_time, num_threads, (num_threads == 1) ? "thread" : "threads");     \
    thread_time;                                                                                        \
})

int main()
{
    int64_t seed = 1, a = 16807, c = 0, m = 2147483647;
    int64_t rand_parameters[] = {seed, a, c, m};

    // Question 1
    printf("\nQ1.\n");

    int32_t number_of_random_samples = 500000;

    double time_rejection = MEASURE_TIME(rejection_sampling, 1.0, number_of_random_samples, rand_parameters);
    double time_direct = MEASURE_TIME(direct_mapping, 1.0, number_of_random_samples, rand_parameters, "direct_mapping_binned_results_50000.csv");

    if (time_rejection < time_direct)
    {
        printf("rejection_sampling was faster by %.6f seconds\n", time_direct - time_rejection);
    }
    else
    {
        printf("direct_mapping was faster by %.6f seconds\n", time_rejection - time_direct);
    }

    // For the nice graphs.
    direct_mapping(10000, rand_parameters, "direct_mapping_binned_results_10000.csv");
    direct_mapping(2500, rand_parameters, "direct_mapping_binned_results_2500.csv");
    generate_expected_density_bins_for_plot_for_report();

    printf("\nQ2.\n");

    // Question 2.

    double thread_count = (double)omp_get_max_threads();
    int32_t total_photon_count = 1000000;
    MEASURE_TIME(photon_scattering, thread_count, total_photon_count, rand_parameters, 0, 200 / 10, 0, "photon_bins.csv", "photon_position_bins.csv");

    // Question 3.
    printf("\nQ3.\n");
    printf("Blue light\n");
    MEASURE_TIME(photon_scattering, thread_count, total_photon_count, rand_parameters, 1, 200 / 10, 0, "rayleigh_blue_bins.csv", "rayleigh_blue_position_bins.csv");
    printf("\nOther colours\n");
    MEASURE_TIME(photon_scattering, thread_count, total_photon_count, rand_parameters, 1, 200 / 0.1, 0, "rayleigh_other_bins.csv", "rayleigh_other_position_bins.csv");

    printf("\n");
    return 0;
}