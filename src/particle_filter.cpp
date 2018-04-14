/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    // TODO: Set standard deviations for x, y, and theta.
    
    
        double std_x, std_y, std_theta;
    
        /// we choose 50 particles
        num_particles=  50;
    
        // save  the standard deviations
        std_x = std[0];
        std_y =  std[1];
        std_theta = std[2];
    
        default_random_engine gen;
    
        // This line creates a normal (Gaussian) distribution for x.
        normal_distribution<double> dist_x(x, std_x);
        //Create normal distributions for y and theta.
        normal_distribution<double> dist_y(y, std_y);
        normal_distribution<double> dist_theta(theta, std_theta);
    
        for (int i = 0; i < num_particles; i++) {
        
            double sample_x, sample_y, sample_theta;
            // TODO: Sample  and from these normal distrubtions like this:
            sample_x = dist_x(gen);
            // where "gen" is the random engine initialized earlier.
            sample_y = dist_y(gen);
            sample_theta = dist_theta(gen);
            
            Particle p ;
            
            // for each particle initialize the x,y and theta
            p.id = i;
            p.x = sample_x;
            p.y = sample_y;
            p.theta = sample_theta;
            
            // set initial weight to 1.0
            p.weight = 1.0;
            particles.push_back(p);
            weights.push_back(p.weight);
        }
        is_initialized =  true;
   

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    double std_x = std_pos[0];
    double std_y =  std_pos[1];
    double std_theta = std_pos[2];
    default_random_engine gen;
   
    for (int i = 0; i < num_particles; i++) {
        Particle p = particles[i];
            
        double pred_x;
        double pred_y;
        double pred_theta;
        if (abs(yaw_rate) > 0.0001)  {
            pred_x= p.x + (velocity / yaw_rate) * (sin(p.theta + (yaw_rate * delta_t))  - sin(p.theta));
            pred_y= p.y + (velocity / yaw_rate) * (cos(p.theta) -cos(p.theta + (yaw_rate * delta_t)));
            pred_theta = p.theta + yaw_rate * delta_t;
        } else {
            pred_x = p.x + velocity * cos(p.theta) * delta_t;
            pred_y = p.y + velocity * sin(p.theta) * delta_t;
            pred_theta = p.theta;
        }
            
        // This line creates a normal (Gaussian) distribution for x.
        normal_distribution<double> dist_x(pred_x, std_x);
            
        // TODO: Create normal distributions for y and theta.
        normal_distribution<double> dist_y(pred_y, std_y);
            
        normal_distribution<double> dist_theta(pred_theta, std_theta);
        particles[i].x = dist_x(gen);
        
        // where "gen" is the random engine initialized earlier.
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations, double sensor_range) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    
    // find the closest neighbor landmark using distance calculation , then assig the predicted landmark 's id to observed land mark
    for (int i =0; i < observations.size(); i++)  {
        double lowest_dist = sensor_range * sqrt(2);
        int closest_landmark_id = -1;
        for (int k =0 ; k < predicted.size(); k++)  {
            double current_dist  = dist(observations[i].x, observations[i].y,predicted[k].x,predicted[k].y);
            if (current_dist < lowest_dist) {
                lowest_dist =  current_dist;
                closest_landmark_id = predicted[k].id;
            }
        }
        observations[i].id = closest_landmark_id;
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    
    /*This variable is used for normalizing weights of all particles and bring them in the range
     of [0, 1]*/
    double weight_normalizer = 0.0;
    
    
    for (int i = 0; i < num_particles ; i++) {
        Particle p =  particles[i];
        
        
        std::vector<LandmarkObs> transformedObservations ;
        
         // transform the observed landmarks from car coordinates to map coordinates
         for (int j =0; j < observations.size(); j++)  {
             LandmarkObs landmark  = observations[j];
             LandmarkObs transformed_obs;
             transformed_obs.id = j;
             transformed_obs.x = p.x + (cos(p.theta) * landmark.x) - (sin(p.theta) * landmark.y);
             // transform to map y coordinate
             transformed_obs.y = p.y + (sin(p.theta) * landmark.x) + (cos(p.theta) * landmark.y) ;
             transformedObservations.push_back(transformed_obs);
         }
        
        // only pick the landm marks from the map which is within the sensor range, as only these can be
        // used for nearest neighbor calculation using sensor observations
        vector<LandmarkObs> predicted_landmarks;
        for (int k = 0; k < map_landmarks.landmark_list.size(); k++) {
            Map::single_landmark_s current_landmark = map_landmarks.landmark_list[k];
            if ((fabs((p.x - current_landmark.x_f)) <= sensor_range) && (fabs((p.y - current_landmark.y_f)) <= sensor_range)) {
                predicted_landmarks.push_back(LandmarkObs {current_landmark.id_i, current_landmark.x_f, current_landmark.y_f});
            }
        }
        
        // assign nearest neighbor
        dataAssociation(predicted_landmarks, transformedObservations, sensor_range);
        
        particles[i].weight = 1.0;
       
        
        double sigma_x = std_landmark[0];
        double sigma_y = std_landmark[1];
        double sigma_x_2 = pow(sigma_x, 2);
        double sigma_y_2 = pow(sigma_y, 2);
        // calculate normalization term
        double gauss_norm = (1.0/(2.0 * M_PI * sigma_x * sigma_y));
        
        // adjust the weights of the landmarks according to the probability of their occurrance
        for (int j =0; j < transformedObservations.size(); j++)  {
             double trans_obs_x = transformedObservations[j].x;
             double trans_obs_y = transformedObservations[j].y;
             int trans_obs_id = transformedObservations[j].id;
             for (int k =0; k < predicted_landmarks.size(); k++)  {
                 double pred_landmark_x = predicted_landmarks[k].x;
                 double pred_landmark_y = predicted_landmarks[k].y;
                 int pred_landmark_id = predicted_landmarks[k].id;
                 if (trans_obs_id == pred_landmark_id)  {
                     // if nearest neighbor, then  adjust the weight by the probablity
                     double multi_prob = gauss_norm * exp(-1.0 * ((pow((trans_obs_x - pred_landmark_x), 2)/(2.0 * sigma_x_2)) + (pow((trans_obs_y - pred_landmark_y), 2)/(2.0 * sigma_y_2))));
                     particles[i].weight *= multi_prob;
                 }
             }
         }
        weight_normalizer += particles[i].weight;
    }
    for (int i = 0; i < particles.size(); i++) {
        particles[i].weight /= weight_normalizer;
        weights[i] = particles[i].weight;
    }
    
    
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    vector<Particle> resampled_particles;
    
    default_random_engine gen;
    
    //Generate random particle index
    uniform_int_distribution<int> particle_index(0, num_particles - 1);
    
    int current_index = particle_index(gen);
    
    
    
    double beta = 0.0;
    
    // re sample the particles with replacement
    double max_weight = 2.0 * *max_element(weights.begin(),weights.end());
    for (int i = 0 ; i < particles.size() ; i++) {
        uniform_real_distribution<double> random_weight(0.0, max_weight);
        beta += random_weight(gen);
        while (beta > weights[current_index]) {
            beta  -= weights[current_index];
            current_index = (current_index + 1) % num_particles;
        }
        resampled_particles.push_back(particles[current_index]);
    }
    particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
