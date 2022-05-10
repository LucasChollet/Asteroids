#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "Utils/Constants.h"
#include "Utils/Gravity.h"
#include "Utils/Integrators.h"
#include "Utils/SpaceMechanics.h"

/**
 * Write to file
 * @param filename
 * @param container1
 * @param container2
 */
static void dump_file(std::string const& filename, std::vector<vec3> const& container1, std::vector<vec3> const& container2) {
    std::ofstream s{filename};
    for (unsigned i{}; i < std::size(container1); ++i)
        s << container1[i].x << ' ' << container1[i].y << ' ' << container1[i].z << ' ' << container2[i].x << ' ' << container2[i].y << ' ' << container2[i].z << '\n';
}

/**
 * Write to file
 * @param filename
 * @param container
 * @param sample_rate_in_year
 */
static void dump_file(std::string const& filename, std::vector<double> const& container, auto sample_rate_in_year) {
    std::ofstream s{filename};
    for (unsigned i{}; i < std::size(container); ++i)
        s << double(i) * sample_rate_in_year << ' ' << container[i] << '\n';
}

int main() {
    // Time of simulation
    double constexpr final_time = 11.9 * year;
    // Step
    double constexpr step = .5;
    // Number of iteration
    unsigned constexpr iteration_nb = 1 + unsigned(final_time / step);
    // Sample rate (at which we save the data)
    unsigned constexpr sample_rate = 26;
    // Number of samples
    unsigned constexpr sample_nb = 1 + iteration_nb / sample_rate;

    double time = 0;

    // Convert initial keplerian elements to cartesian
    Body VW_266 = kep2cart(aA, eA, iA, wA, WA, time, M0A, T0A);

    // Create vectors in which data will be stored
    std::vector<double> semi_major_axis(sample_nb);
    std::vector<vec3> positions(sample_nb);
    std::vector<vec3> positions_jup(sample_nb);
    std::vector<double> eccentricity(sample_nb);
    std::vector<double> inclination(sample_nb);

    std::cout << std::setprecision(2) << std::fixed;
    for (unsigned i = 0; i < iteration_nb; ++i) {
        // Apply Runge-Kutta on the system
        VW_266 = runge_kutta4(VW_266, time, step, derivate);

        // Store data
        if (i % sample_rate == 0) {
            auto const sample_index = i / sample_rate;
            positions[sample_index] = VW_266.position;
            positions_jup[sample_index] = kep2cart(aJ, eJ, iJ, wJ, WJ, time, M0J, T0J).position;
            semi_major_axis[sample_index] = compute_semi_major_axis(VW_266);
            eccentricity[sample_index] = compute_eccentricity(VW_266);
            inclination[sample_index] = to_degrees(compute_inclination(VW_266));
        }

        if (i % (iteration_nb / 1000) == 0)
            std::cout << "\r" << 100 * double(i) / iteration_nb << " %" << std::flush;

        time += step;
    }
    std::cout << "\r100.00 %" << std::endl;

    // Write data to file
    dump_file("position.txt", positions_jup, positions);
    dump_file("semi_major_axis.txt", semi_major_axis, double(sample_rate) * step / year);
    dump_file("eccentricity.txt", eccentricity, double(sample_rate) * step / year);
    dump_file("inclination.txt", inclination, double(sample_rate) * step / year);
}