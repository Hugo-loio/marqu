#include "Runner.h"

#include <algorithm>
#include <stdexcept>

//marqu::Runner::Runner(BaseParticleSimulator & sim) : simulator(sim){}

void marqu::Runner::run(double T, std::size_t nSamples){
  resultClear();
  resultInitialize();
  nInitialParticles = nsamples;

  std::vector<void (Runner::*)()> updaters;
  updaters.push_back(& Runner::updateObservables);
  if(options.saveParticleNumber) 
    updaters.push_back(& Runner::updateParticleNumber);
  if(options.saveCompressionRate) 
    updaters.push_back(& Runner::updateCompressionRate);
  if(options.saveRuntime) updaters.push_back(& Runner::updateRuntime);
  if(options.saveMaxParticles) updaters.push_back(& Runner::updateMaxParticles);

  void (Runner::*)() timeStep = nullptr;
  if(options.stepMethod == 0) timeStep = & Runner::gillespieTimeStep;
  else if(options.stepMethod == 1) timeStep = & Runner::discreteTimeStep;
  else throw std::invalid_argument("Invalid time step method option");

  for(int i = 0; i < nsamples; i++){
    nInitialParticles += simulator.initialize(options.initialParticleNumber, 
	options.removeStaticConfigs);

    resetSample(T);

    while(t < T){
      for (auto updater : updaters) (this->*updater)(); 
      tPrevious = t;
      (this->*timeStep)();
    }
    t = T;
    for (auto updater : updaters) (this->*updater)(); 
    auto end = std::chrono::high_resolution_clock::now();

    //TODO: Compute averaged results
    //In utils can have a function that converts marqu::histograms to NDArrays

    statistics[0] += (double)maxParticles/nsamples; // Max particle number
    statistics[1] += std::chrono::duration<double, std::milli>(end - start).count()/nsamples; // Elapsed time

    //avgHistogram += histogram.averagedWeights()/nsamples;
    Lasap::NDArray<double> avgHistogramMag1 = histogram_mag1.averagedWeights();
    Lasap::NDArray<double> avgHistogramCorr = histogram_corr.averagedWeights();
    for(int i = 0; i < nsavedts; i++){
      avgHistogram[std::vector<int>({0,i})] += avgHistogramMag1[i];
      avgHistogram[std::vector<int>({1,i})] += avgHistogramCorr[i];
    }
    avgHistogram_number += histogram_number.averagedWeights()/nsamples;
    avgHistogram_comp_rate += histogram_comp_rate.averagedWeights()/nsamples;

    //progress.print_progress(count);
  }
  avgHistogram /= nInitialParticles;

}

void marqu::Runner::resultClear(){
  avgHistObservables.clear();
  avgHistCompression.reset();
  avgHistParticleNumber.reset();
  avgMaxParticles.reset();
  avgRuntime.reset();
}

voir marqu::Runner::resultInitialize(){
  if(saveMaxParticles) avgMaxParticles.emplace(0);
  if(saveRuntime) avgRuntime.emplace(0);
}

void marqu::Runner::resetSample(double T){
  start = std::chrono::steady_clock::now();
  tPrevious = 0;
  t = 0;
  maxParticles = 0;
  compRate = 1;
  floatParticleNumber = 1;
  observables = simulator.observableEstimate();
  nObservables = observables.size();
  histObservables.clear();
  for(std::size_t i = 0; i < nObservables; i++){
    histObservables.emplace_back(0, T, options.nBins);
  }
  if(options.saveParticleNumber) histParticleNumber.emplace(0, T, options.nBins); 
  if(options.saveCompressionRate) histCompression.emplace(0, T, options.nBins); 
}

void marqu::Runner::updateObservables(){
  for(std::size_t i = 0; i < nObservables; i++){
    histObservables[i].add(tPrevious, t, observables[i]);
  }
  std::copy(observables.begin(), observables.end(),
      simulator.observableEstimate().begin());
}

void marqu::Runner::updateParticleNumber(){
  histParticleNumber.add(tPrevious, t, floatParticleNumber);
  floatParticleNumber = (double)simulator.getParticleNumber();
}

void marqu::Runner::updateCompressionRate(){
  histCompression.add(tPrevious, t, compRate);
  compRate = simulator.compressionRate();
}

void marqu::Runner::updateMaxParticles(){
  maxParticles = std::max(maxParticles, simulator.getParticleNumber());
}

void marqu::Runner::gillespieTimeStep(){
  t += simulator.gillespieTimeStep();
}

void marqu::Runner::discreteTimeStep(){
  t += options.dt;
  simulator.discreteTimeStep(dt);
}
