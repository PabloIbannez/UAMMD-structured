#include <iostream>
#include <cmath>
#include <random>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <fstream>

// Function to generate a random rotation matrix with a small angle
// given a maximum rotation angle maxAngle.
Eigen::Matrix3d randomSmallRotation(std::mt19937 &gen, double maxAngle) {
    std::uniform_real_distribution<double> distUniform(0.0, 1.0);
    // Generate a random unit vector for the rotation axis
    double u = distUniform(gen);
    double costheta = 1.0 - 2.0*u;
    if (costheta > 1.0)  {costheta =  1.0;}
    if (costheta < -1.0) {costheta = -1.0;}
    double theta = std::acos(costheta);
    double phi = 2.0 * M_PI * distUniform(gen);
    double x = std::sin(theta)*std::cos(phi);
    double y = std::sin(theta)*std::sin(phi);
    double z = std::cos(theta);

    // Generate a small angle from a uniform distribution in [-maxAngle, maxAngle]
    std::uniform_real_distribution<double> distAngle(-maxAngle, maxAngle);
    double angle = distAngle(gen);

    // Construct rotation matrix using Rodrigues' rotation formula
    Eigen::Vector3d axis(x, y, z);
    axis.normalize();
    Eigen::Matrix3d K;
    K << 0, -axis.z(), axis.y(),
         axis.z(), 0, -axis.x(),
        -axis.y(), axis.x(), 0;

    Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + std::sin(angle)*K + (1.0 - std::cos(angle))*(K*K);
    return R;
}

// Compute the potential energy E = (K/4)*(3 - tr(B^T A R))
double energy(const Eigen::Matrix3d &A, const Eigen::Matrix3d &B, const Eigen::Matrix3d &R, double K) {
    double trace_val = (B.transpose() * A * R).trace();
    double E = (K/4.0) * (3.0 - trace_val);
    return E;
}

double polymerEnergy(const std::vector<Eigen::Matrix3d> &particles, const Eigen::Matrix3d &R, double K) {
    double E = 0.0;
    for (int i = 0; i < particles.size()-1; i++) {
        E += energy(particles[i], particles[i+1], R, K);
    }
    return E;
}

int main(int argc, char **argv) {

    int N;

    double T;

    double K;
    double q0, q1, q2, q3;

    int steps;
    int equilibration;
    int output;

    std::string filename;

    if (argc != 12) {
        std::cerr << "Usage: " << argv[0] << " N T K q0 q1 q2 q3 steps equilibration output filename" << std::endl;
        return 1;
    } else {
        N = std::stoi(argv[1]);

        T = std::stod(argv[2]);

        K = std::stod(argv[3]);
        q0 = std::stod(argv[4]);
        q1 = std::stod(argv[5]);
        q2 = std::stod(argv[6]);
        q3 = std::stod(argv[7]);

        steps         = std::stoi(argv[8]);
        equilibration = std::stoi(argv[9]);
        output        = std::stoi(argv[10]);

        filename = argv[11];
    }

    steps         = steps * N; // Convert steps to MC steps
    equilibration = equilibration * N; // Convert equilibration to MC steps

    std::cerr << "N = " << N << std::endl;
    std::cerr << "T = " << T << std::endl;
    std::cerr << "K = " << K << std::endl;
    std::cerr << "qR = " << q0 << " " << q1 << " " << q2 << " " << q3 << std::endl;
    std::cerr << "steps = " << steps << std::endl;
    std::cerr << "equilibration = " << equilibration << std::endl;
    std::cerr << "output = " << output << std::endl;
    std::cerr << "filename = " << filename << std::endl;

    Eigen::Quaterniond qR(q0, q1, q2, q3);
    Eigen::Matrix3d R = qR.toRotationMatrix();

    double beta = 1.0 / T;

    // Initialize configurations (start with identity matrices)
    std::vector<Eigen::Matrix3d> particles(N);
    for (int i = 0; i < N; i++) {
        particles[i] = Eigen::Matrix3d::Identity();
    }

    // RNG Setup
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distUniform(0.0, 1.0);

    // Initial Energy
    double E_current = polymerEnergy(particles, R, K);

    // MC parameters
    double maxAngle = 0.1; // maximum rotation angle for trial moves
    // Equilibration phase
    for (int i = 0; i < equilibration; i++) {
        // Choose which particle to rotate
        int particleIndex = std::floor(distUniform(gen) * N);

        Eigen::Matrix3d trial = particles[particleIndex];
        trial = randomSmallRotation(gen, maxAngle) * trial;

        Eigen::Matrix3d old = particles[particleIndex];
        particles[particleIndex] = trial;

        double E_new = polymerEnergy(particles, R, K);
        double dE = E_new - E_current;

        // Metropolis acceptance
        if (dE < 0.0 || distUniform(gen) < std::exp(-beta * dE)) {
            // Accept move
            E_current = E_new;
        } else {
            // Reject move
            particles[particleIndex] = old;
        }
    }

    std::ofstream file;
    file.open(filename);

    // Production phase: after equilibration, sample configurations
    // Here we just run and possibly print out some configurations.
    // In practice, you might store these configurations or compute averages.
    for (int i = 0; i < steps; i++) {
        // Choose which particle to rotate
        int particleIndex = static_cast<int>(distUniform(gen) * N);

        Eigen::Matrix3d trial = particles[particleIndex];
        trial = randomSmallRotation(gen, maxAngle) * trial;

        Eigen::Matrix3d old = particles[particleIndex];
        particles[particleIndex] = trial;

        double E_new = polymerEnergy(particles, R, K);
        double dE = E_new - E_current;

        // Metropolis acceptance
        if (dE < 0.0 || distUniform(gen) < std::exp(-beta * dE)) {
            // Accept move
            E_current = E_new;
        } else {
            // Reject move
            particles[particleIndex] = old;
        }

        // Output
        if (i % output == 0) {
            std::cerr << "Completed " << float(i) / steps * 100 << "% of the simulation." << std::endl;
            file << N << std::endl;
            file << std::endl;
            // Format: id type x y z q0 q1 q2 q3
            // Convert rotation matrix to quaternion
            for (int j = 0; j < N; j++) {
                Eigen::Quaterniond q(particles[j]);
                file << j << " 0 0 0 0 " << q.w() << " " << q.x() << " " << q.y() << " " << q.z() << std::endl;
            }
        }
    }

    return 0;
}

