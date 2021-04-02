/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using namespace std;
using std::normal_distribution;


// declare a random engine to be used across multiple and various method calls
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  // random number generator
  std::default_random_engine gen;
  
  // declaring the standard deviations 
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];
  
  // creatign a normal distrubution for the input values 
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);
  
  
  num_particles = 10;  // TODO: Set the number of particles
  
  //creating a vector of particles
  vector<Particle> particles_init(num_particles);
  
  for (int i= 0;i<num_particles; i++){
    
    
    
    // Declaring the Particle 
    Particle particle;
    
    // get the samples 
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    
    // weights to 1 
    particle.weight = 1.0;
    
    // add the particle to particle set 
    particles_init[i]=particle;
  }
  particles = particles_init;
  std::cout<<"Particles: "<<particles.size()<<" ,Initialization Done"<<"\n";

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  
  // random number generator
  std::default_random_engine gen;
  
    // declaring the standard deviations 
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2];
  
  // creatign a normal distrubution for the input values 
  normal_distribution<double> dist_x(0, std_x);
  normal_distribution<double> dist_y(0, std_y);
  normal_distribution<double> dist_theta(0, std_theta);
  
  for (int i= 0;i<num_particles; i++){
    
 	double yaw_angle = particles[i].theta + yaw_rate * delta_t;
    // prediction the equations
    particles[i].x = particles[i].x +velocity/yaw_rate*(sin(yaw_angle+yaw_rate*delta_t)-sin(yaw_angle));
    particles[i].y = particles[i].y +velocity/yaw_rate*(cos(yaw_angle)+cos(yaw_angle+yaw_rate*delta_t));
    particles[i].theta += yaw_rate * delta_t;
    
    // Adding noise 
    particles[i].x +=dist_x(gen) ;
    particles[i].y +=dist_y(gen) ;
    particles[i].theta +=dist_theta(gen) ;
   
  }
  
  // Calculate the prediction values 
  std::cout<<"Particles: "<<particles.size()<<" ,Prediction values Done"<<"\n";

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a 
   
   
   *   during the updateWeights phase.
   */
    // calculate the distance to each landmark and select the landmark having the least distance 
    double min_distance=numeric_limits<double>::max();
    double landmark_id_temp;

    for (unsigned int j = 0;j<observations.size();j++){
      for (unsigned int k=0; k<predicted.size();k++){
      	double distance_calc = sqrt(pow(predicted[k].x-observations[j].x,2)+pow(predicted[k].y-observations[j].x,2));
        //std::cout<<distance_calc<<"\n";
        
        if(distance_calc<min_distance){
        	min_distance = distance_calc;
          	//std::cout<<distance_calc<<"\n";
          	
          	landmark_id_temp = predicted[k].id;
        }
        
        
      }
      // set the particle association to the landmark id
      observations[j].id = landmark_id_temp;
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  // get the landmark list 
  //std::vector<single_landmark_s>landmarks = map_landmarks.landmark_list;
  
  // looping over the particles 
  std::cout<<particles.size()<<"\n";
  for (unsigned int i=0;i<particles.size();i++){
    
    double particle_x = particles[i].x;
    double particle_y = particles[i].y;
    //double particle_theta = particles[i].theta;
    
    // Select the landmarks inside the sensor_range 
    vector<LandmarkObs> predictions;
    
    for (unsigned int l=0;l<map_landmarks.landmark_list.size();l++){
    	
      float select_x = map_landmarks.landmark_list[l].x_f;
      float select_y = map_landmarks.landmark_list[l].y_f;
      int select_id = map_landmarks.landmark_list[l].id_i;
      double landmark_range = sqrt(pow(select_x-particle_x,2)+pow(select_y-particle_y,2));
      
      if(landmark_range<=sensor_range){
      	predictions.push_back(LandmarkObs{select_id,select_x,select_y});
      }
      
    
    }
    
    //set weight to 1
  	particles[i].weight = 1;
    
    // Select the landmarks inside the sensor_range 
    vector<LandmarkObs> MappedObs;
    
    
    // looping over observation 
    for (unsigned int j=0;j<observations.size();j++){
    	
      // transform the observation from car to map cordinates [Homogenous transformation]
      
      double map_x = particles[i].x + observations[j].x*cos(particles[i].theta)-observations[j].y*sin(particles[i].theta);
      double map_y = particles[i].y + observations[j].x*sin(particles[i].theta)+observations[j].y*cos(particles[i].theta);
      MappedObs.push_back(LandmarkObs{observations[j].id,map_x,map_y});
     
    }
    // Associate the Mpa corrdiante observation with the landmark observation in the sensor range
    dataAssociation(predictions,MappedObs);
    
    // Loop through the associated observations
    for (unsigned int j=0;j<MappedObs.size();j++){
      
      // Weight contribution calculation from bivariate gaussian distribution 
      // Fucntion for the bivariate gaussian distribution is adapted in the helper fucntions
      double sig_x  =std_landmark[0];
      double sig_y  =std_landmark[1];
      particles[i].weight*=calculateWeightsMultiVarGauss(sig_x,sig_y,particles[i].x,particles[i].y,MappedObs[j].x,MappedObs[j].y); 

    }
    cout<<"Weight of particle: "<<i<<": "<<particles[i].weight<<"\n";
    
  }
std::cout<<"Update Weigths Done"<<"\n";
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

	// Get weights and max weight.
    //cout<<num_particles<<"\n";
  // Get weights and max weight.
  vector<double> weights;
  double maxWeight = numeric_limits<double>::min();
  for(int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
    if ( particles[i].weight > maxWeight ) {
      maxWeight = particles[i].weight;
    }
  }

  // Creating distributions.
  uniform_real_distribution<double> distDouble(0.0, maxWeight);
  uniform_int_distribution<int> distInt(0, num_particles - 1);

  // Generating index.
  int index = distInt(gen);

  double beta = 0.0;

  // the wheel
  vector<Particle> resampledParticles;
  for(int i = 0; i < num_particles; i++) {
    beta += distDouble(gen) * 2.0;
    while( beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    resampledParticles.push_back(particles[index]);
  }

  particles = resampledParticles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}