#ifndef BASEPARTICLESIMULATOR_H
#define BASEPARTICLESIMULATOR_H

#include <vector>
#include <random>
#include <string>

#include "Configuration.h"
#include "Particle.h"
#include "Model.h"

namespace marqu{
  class BaseParticleSimulator{
    public:
      BaseParticleSimulator() : gen(std::random_device{}()){};
      ~BaseParticleSimulator();

      void set_seed(std::mt19937::result_type seed) {gen = std::mt19937(seed);};

      int getParticleNumber() const {return particleNumber;};

      virtual std::size_t getObservableCount() const = 0;
      virtual void observables(const Configuration &, std::vector<double> & out) = 0;
      const std::vector<double> & observableEstimate() const{return observableTracker;};

      virtual void discreteTimeStep(double dt) = 0;
      virtual double gillespieTimeStep() = 0;

      virtual void displayParticles() const = 0;
      virtual double compressionRate() {return 1;};

      // Only product states are supported for now
      void setInitialState(const std::string & orientations); 
      int initialize(int particleNumber, bool removeStatic = true); 

      void setModel(const Model & model);
      void setModel(Model && model);

    protected:
      int particleNumber = 0;
      int initParticleNumber = 1;
      double totalRate = 0;
      std::vector<double> observableTracker;
      std::vector<double> observableBuffer;

      virtual double eventRate(const Configuration & configuration) const; 
      //The sign refers to the resulting M+ or M-
      //Not const due to the random number generator
      virtual std::pair<Configuration,Sign> randomEvent(const Particle & particle); 

      // If (add == false) then subtract (when particles are moved or deleted)
      void updateObservable(const Particle & particle, bool add);

      void trackClear();
      void trackAdd(const Particle & particle);
      void trackRemove(const Particle & particle);
      void trackMove(const Particle & pIn, const Particle & pOut);
      void trackAnnihilate(int nAnnihilations, double configRate);

      virtual void clearParticles() = 0; 
      virtual void addParticle(Particle && particle) = 0;
      virtual void addParticle(const Particle & particle) = 0;

      std::mt19937 gen;
      std::uniform_real_distribution<double> uni_dist = std::uniform_real_distribution<double>(0,1); 

      Configuration * initConfig = nullptr;
      Model * model = nullptr;
  };
}

#endif
