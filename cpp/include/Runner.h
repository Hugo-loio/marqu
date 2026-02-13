#ifndef RUNNER_H
#define RUNNER_H

#include <chrono>
#include <optional>
#include <vector>

#include "BaseParticleSimulator.h"
#include "Histogram.h"

namespace marqu{
  class Runner{
    public:
      Runner(BaseParticleSimulator & sim) : simulator(sim) {}; 
      void run(double T, std::size_t nSamples);

      struct Options {
	int initialParticleNumber = 1;
	bool removeStaticConfigs = true;
	bool saveParticleNumber = true;
	bool saveCompressionRate = true;
	bool saveRuntime = true;
	bool saveMaxParticles = true;
	bool progress = true;
	int progressDivisor = 10;
	int stepMethod = 0; // 0 - Gillespie, 1 - Discrete
	std::size_t nBins = 100;
	double dt = 0.1; // In case of discrete time step method
      };

      Options options;

      // Sample averaged results
      std::vector<Histogram> avgHistObservables;
      std::optional<Histogram> avgHistParticleNumber;
      std::optional<Histogram> avgHistCompression;
      std::optional<double> avgMaxParticles;
      std::optional<double> avgRuntime;

    protected:
      void resultClear();
      void resultInitialize();
      void resetSample(double T);
      void addSample();

      void updateObservables();
      void updateParticleNumber();
      void updateCompressionRate();
      void updateRuntime();
      void updateMaxParticles();

      void gillespieTimeStep();
      void discreteTimeStep();

      BaseParticleSimulator & simulator;
      int nInitialParticles;

      //Sample local attributes
      std::chrono::steady_clock start;
      double tPrevious;
      double t;
      int maxParticles;
      double floatParticleNumber;
      double compRate; 
      vector<double> observables;
      std::size_t nObservables;
      std::vector<Histogram> histObservables;
      std::optional<Histogram> histParticleNumber;
      std::optional<Histogram> histCompression;
  };
}

#endif
