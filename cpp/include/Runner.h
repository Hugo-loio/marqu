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
	bool saveMaxParticles = true;
	bool saveRuntime = true;
	bool progress = true; 
	int progressDivisor = 10;
	int stepMethod = 0; // 0 - Gillespie, 1 - Discrete
	std::size_t nBins = 100;
	// This might not make much sense, have instead T/nBins
	double dt = 0.1; // In case of discrete time step method
      };

      Options options;

      // Sample averaged results
      std::vector<std::vector<double>> avgHistObservables;
      std::optional<std::vector<double>> avgHistParticleNumber;
      std::optional<std::vector<double>> avgHistCompression;
      std::optional<double> avgMaxParticles;
      std::optional<double> avgRuntime;
      std::vector<double> times;

    protected:
      void resultClear();
      void resultInitialize(double T);
      void resultRenormalize(std::size_t nSamples, 
	  std::size_t nInitialParticles);

      void resetSample(double T);
      void addSample();
      void printProgress(std::size_t count, std::size_t progressInterval, 
	  std::size_t nSamples);

      void addAvgHist(std::vector<double> & target, 
	  const Histogram<double> & source);
      void renormalize(std::vector<double> & hist, double norm);

      void update();
      void updateObservables();
      void updateParticleNumber();
      void updateCompressionRate();
      void updateMaxParticles();

      void timeStep();
      void gillespieTimeStep();
      void discreteTimeStep();

      BaseParticleSimulator & simulator;
      std::chrono::steady_clock::time_point runStart;

      //Sample local attributes
      std::chrono::steady_clock::time_point sampleStart;
      double tPrevious;
      double t;
      int maxParticles;
      double floatParticleNumber;
      double compRate; 
      std::vector<double> observables;
      std::size_t nObservables;
      std::vector<Histogram<double>> histObservables;
      std::optional<Histogram<double>> histParticleNumber;
      std::optional<Histogram<double>> histCompression;
  };
}

#endif
