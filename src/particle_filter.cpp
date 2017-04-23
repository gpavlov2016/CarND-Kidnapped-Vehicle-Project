/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>

#include "particle_filter.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace std;

void ParticleFilter::print_best_particle()
{
    double highest_weight = 0.0;
    Particle best_particle;
    for (int i = 0; i < num_particles; ++i) {
        if (particles[i].weight > highest_weight) {
            highest_weight = particles[i].weight;
            best_particle = particles[i];
        }
    }

    cout << "best x: " << best_particle.x << ", y: " << best_particle.y << ", theta: " << best_particle.theta << ", w: " << best_particle.weight << endl;
}
void ParticleFilter::print_all_particles()
{
    Particle best_particle;
    for (int i = 0; i < num_particles; ++i) {
        best_particle = particles[i];
        cout << "p" << i << " x: " << best_particle.x << ", y: " << best_particle.y << ", theta: " << best_particle.theta << ", w: " << best_particle.weight << endl;
    }

}


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    num_particles = 100;
    weights.resize(num_particles,0);

	default_random_engine gen;
	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

    for (int i=0; i < num_particles; i++) 
    {
		double sample_x, sample_y, sample_theta;
        sample_x = dist_x(gen);
        sample_y = dist_y(gen);
        sample_theta = dist_theta(gen);	 

        Particle particle;
        particle.id = i;
        particle.x = sample_x;
        particle.y = sample_y;
        particle.theta = sample_theta;
        particle.weight = 1;
        particles.push_back(particle);
    }
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    default_random_engine gen;

    for (int i=0; i < num_particles; i++) 
    {
        Particle& particle = particles[i];
        double theta_0 = particle.theta;
        
        double theta_1 = theta_0 + delta_t * yaw_rate;
        normal_distribution<double> dist_theta(theta_1, std_pos[2]);
        particle.theta = dist_theta(gen);

        double x = particle.x + (velocity / yaw_rate)*(sin(theta_1) - sin(theta_0));
        normal_distribution<double> dist_x(x, std_pos[0]);
        particle.x = dist_x(gen);

        double y = particle.y + (velocity / yaw_rate)*(cos(theta_0) - cos(theta_1));
        normal_distribution<double> dist_y(y, std_pos[1]);
        particle.y = dist_y(gen);
    }
    //cout << "particles[0].x: " << particles[0].x << ", 1: " << particles[1].x << endl;
    //print_best_particle();
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

    double sigmax = std_landmark[0];
    double sigmay = std_landmark[1];
    double w_sum = 0;
    for (int i=0; i < num_particles; i++) 
    {
        Particle& particle = particles[i];
        double x = particle.x;
        double y = particle.y;
        double theta = particle.theta;
        //cout << "particle x: " << x << ", y: " << y << ", theta: " << theta << endl;;

        double weight = 1.0;
        for (int o=0; o<observations.size(); o++) 
        {
            LandmarkObs obs = observations[o];
            //Convert to map coordinates:
            double o_x = x + obs.x*cos(theta)-obs.y*sin(theta);
            double o_y = y + obs.x*sin(theta)+obs.y*cos(theta);
            //cout << "observation x: " << obs.x << ", y: " << obs.y << endl;
            //cout << "rotated observation o_x: " << o_x << ", o_y: " << o_y << endl;
            //define pseudo observation vector:
            double min_dist = sensor_range;

            double d_x_min = 0;
            double d_y_min = 0;
            //loop over number of landmarks and calculate ranges:
            for (unsigned int l = 0; l < map_landmarks.landmark_list.size(); ++l) 
            {
                double l_x = map_landmarks.landmark_list[l].x_f;
                double l_y = map_landmarks.landmark_list[l].y_f;
                double d_x = l_x - o_x;
                double d_y = l_y - o_y;
                double dist = sqrt(d_x*d_x + d_y*d_y);
                if (dist < min_dist) 
                {
                    min_dist = dist;
                    d_x_min = d_x;
                    d_y_min = d_y;
                }
                //cout << "landmark x: " << l_x << ", y: " << l_y << endl;
            }

            double dist = sqrt(d_x_min*d_x_min + d_y_min*d_y_min);
            //cout << "dist: " << dist << ", sensor_range:" << sensor_range << endl;
            if (dist < sensor_range) 
            {
                //cout << "d_x_min: " << d_x_min << ", d_y_min: " << d_y_min << endl;
                double cur_weight = (1/(2*M_PI*sigmax*sigmay)) * exp(- ( (d_x_min*d_x_min)/(2*sigmax*sigmax) + (d_y_min*d_y_min)/(2*sigmay*sigmay) ));
                //cout << "cur_weight: " << cur_weight << ", weight: " << weight << endl;
                weight *= cur_weight;
            }
        }
        particle.weight = weight;
        weights[i] = weight;
        //print_all_particles();
    }
    //print_all_particles();
}
#include <random>
void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<int> d(weights.begin(), weights.end());

    std::vector<Particle> resampled;
    for (int i=0; i < num_particles; i++) 
    {
        int rnd = d(gen);
        //cout << "rnd: " << rnd << ", x:" << particles[rnd].x << ", y: " << particles[rnd].y << endl;

        resampled.push_back(particles[rnd]);
    }
    particles = resampled;
    //cout << "particles[0].x: " << particles[0].x << ", 1: " << particles[1].x << endl;
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
