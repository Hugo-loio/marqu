#include "Runner.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <iostream>

void marqu::Runner::run(double T, std::size_t nSamples){
  resultClear();
  resultInitialize(T);
  std::size_t nInitialParticles = nSamples * options.initialParticleNumber;
  std::size_t progressInterval = std::ceil( static_cast<double>(nSamples)/options.progressDivisor);

  for(std::size_t i = 1; i <= nSamples; i++){
    nInitialParticles += simulator.initialize(options.initialParticleNumber, 
	options.removeStaticConfigs);

    resetSample(T);

    try {
      while(t < T){
	update();
	tPrevious = t;
	timeStep();
	//std::cout << std::setprecision(8) << " t = " << t << std::endl;
	//simulator.displayParticles();
      }
    } catch (const std::runtime_error & e) {
      if (std::string(e.what()) != "All particles are stationary") {
	throw; // Re-throw unrelated runtime errors
      }
    }
    t = T;
    update();

    addSample();
    printProgress(i, progressInterval, nSamples);
  }
  resultRenormalize(nSamples, nInitialParticles);
}

void marqu::Runner::resultClear(){
  avgHistObservables.clear();
  avgHistParticleNumber.reset();
  avgHistCompression.reset();
  avgMaxParticles.reset();
  avgRuntime.reset();
  histObservables.clear();
  times.clear();
}

void marqu::Runner::resultInitialize(double T){
  if(options.saveParticleNumber) avgHistParticleNumber.emplace(options.nBins, 0);
  if(options.saveParticleNumber) histParticleNumber.emplace(0, T, options.nBins);
  if(options.saveCompressionRate) avgHistCompression.emplace(options.nBins, 0);
  if(options.saveCompressionRate) histCompression.emplace(0, T, options.nBins);
  if(options.saveMaxParticles) avgMaxParticles.emplace(0);
  if(options.saveRuntime) avgRuntime.emplace(0);

  const double dt = T / options.nBins;
  times.resize(options.nBins);
  for(std::size_t i = 0; i < options.nBins; i++){
    times[i] = (0.5 + i)*dt;
  }

  runStart = std::chrono::steady_clock::now();
}

void marqu::Runner::resultRenormalize(std::size_t nSamples, 
    std::size_t nInitialParticles){
  double norm = 1.0 / nInitialParticles;
  for(std::size_t i = 0; i < nObservables; i++){
    renormalize(avgHistObservables[i], norm);
  }
  norm = 1.0/nSamples;
  if(options.saveParticleNumber) renormalize(*avgHistParticleNumber, norm);
  if(options.saveCompressionRate) renormalize(*avgHistCompression, norm);
  if(options.saveMaxParticles) *avgMaxParticles *= norm;
  if(options.saveRuntime) *avgRuntime *= norm;
}


void marqu::Runner::update(){
  updateObservables();
  if(options.saveParticleNumber) updateParticleNumber();
  if(options.saveCompressionRate) updateCompressionRate();
  if(options.saveMaxParticles) updateMaxParticles();
}


void marqu::Runner::resetSample(double T){
  sampleStart = std::chrono::steady_clock::now();
  tPrevious = 0;
  t = 0;
  maxParticles = 0;
  compRate = 1;
  floatParticleNumber = 1;
  observables = simulator.observableEstimate();
  nObservables = observables.size();
  //histObservables.clear();
  if(avgHistObservables.empty()){
    for(std::size_t i = 0; i < nObservables; i++){
      avgHistObservables.emplace_back(options.nBins, 0);
      histObservables.emplace_back(0, T, options.nBins);
    }
  }
  else{
    for(auto & h : histObservables) h.clear();
  }
  //if(options.saveParticleNumber) histParticleNumber.emplace(0, T, options.nBins); 
  if(options.saveParticleNumber) histParticleNumber->clear(); 
  //if(options.saveCompressionRate) histCompression.emplace(0, T, options.nBins); 
  if(options.saveCompressionRate) histCompression->clear(); 
}

void marqu::Runner::addSample(){
  for(std::size_t i = 0; i < nObservables; i++){
    addAvgHist(avgHistObservables[i], histObservables[i]);
  }
  if(options.saveParticleNumber) 
    addAvgHist(*avgHistParticleNumber, *histParticleNumber);
  if(options.saveCompressionRate) 
    addAvgHist(*avgHistCompression, *histCompression);
  if(options.saveMaxParticles) *avgMaxParticles += (double)maxParticles;
  if(options.saveRuntime) {
    auto end = std::chrono::steady_clock::now();
    *avgRuntime += 
      std::chrono::duration<double, std::milli>(end - sampleStart).count();
  }
}

void marqu::Runner::printProgress(std::size_t count, 
    std::size_t progressInterval, std::size_t nSamples){
  if(!options.progress) return;
  if(count % progressInterval == 0 || count == nSamples){
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
	std::chrono::steady_clock::now() - runStart).count();

    int hours = static_cast<int>(elapsed / (1000 * 3600));
    int minutes = static_cast<int>((elapsed / (1000 * 60)) % 60);
    int seconds = static_cast<int>((elapsed / 1000) % 60);

    std::cout << std::fixed << std::setprecision(0)
      << round(100.0 * count / nSamples) << " % at " << hours << ":"
      << std::setw(2) << std::setfill('0') << minutes << ":"
      << std::setw(2) << std::setfill('0') << seconds << std::endl;
  }
}

void marqu::Runner::addAvgHist(std::vector<double> & target,
    const Histogram<double> & source){
  //const std::vector<double> & avgHist = source.averagedWeights();
  std::vector<double> avgHist = source.averagedWeights();
  for (std::size_t i = 0; i < options.nBins; i++) {
    target[i] += avgHist[i];
  }
}

void marqu::Runner::renormalize(std::vector<double> & hist, double norm){
  for (std::size_t i = 0; i < options.nBins; i++) {
    hist[i] *= norm;
  }
}

void marqu::Runner::updateObservables(){
  for(std::size_t i = 0; i < nObservables; i++){
    histObservables[i].add(tPrevious, t, observables[i]);
  }
  //std::copy(simulator.observableEstimate.begin(), 
  //    simulator.observableEstimate.end(), observables.begin());
  observables = simulator.observableEstimate(); 
}

void marqu::Runner::updateParticleNumber(){
  histParticleNumber->add(tPrevious, t, floatParticleNumber);
  floatParticleNumber = (double)simulator.getParticleNumber();
}

void marqu::Runner::updateCompressionRate(){
  histCompression->add(tPrevious, t, compRate);
  compRate = simulator.compressionRate();
}

void marqu::Runner::updateMaxParticles(){
  maxParticles = std::max(maxParticles, simulator.getParticleNumber());
}

void marqu::Runner::timeStep(){
  if(options.stepMethod == 0) t += simulator.gillespieTimeStep();
  else if(options.stepMethod == 1){
    t += options.dt;
    simulator.discreteTimeStep(options.dt);
  }
  else throw std::invalid_argument("Invalid time step method option");
}

void marqu::Runner::gillespieTimeStep(){
  t += simulator.gillespieTimeStep();
}

void marqu::Runner::discreteTimeStep(){
  t += options.dt;
  simulator.discreteTimeStep(options.dt);
}
