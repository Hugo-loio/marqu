#include <stdexcept>
#include <iostream>

#include "Configuration.h"

marqu::Axis marqu::toAxis(char axis){
  switch (axis) {
    case 'x': return Axis::x;
    case 'y': return Axis::y;
    case 'z': return Axis::z;
    default:
	      throw std::invalid_argument("Invalid axis value: must be x, y, or z.");
  }
}

char marqu::toChar(Axis axis){
  switch (axis) {
    case Axis::x: return 'x';
    case Axis::y: return 'y';
    case Axis::z: return 'z';
    default: throw std::invalid_argument("Invalid Axis value");
  }
}

char marqu::toChar(Sign sign){
  switch (sign) {
    case  Sign::plus: return '+';
    case  Sign::minus: return '-';
    default: throw std::invalid_argument("Invalid Sign value");
  }
}

marqu::Sign marqu::toSign(char sign){
  switch (sign) {
    case '+': return Sign::plus;
    case '-': return Sign::minus;
    default:
	      throw std::invalid_argument("Invalid sign value: must be + or -.");
  }
}

int marqu::toInt(Sign sign, Axis axis){
  return static_cast<int>(sign) + 2*static_cast<int>(axis);
}

int marqu::leviCivita(Axis a1, Axis a2, Axis a3){
  if(a1 == a2 || a2 == a3 || a3 == a1) return 0;
  if( (a1 == Axis::x && a2 == Axis::y && a3 == Axis::z) ||
      (a1 == Axis::y && a2 == Axis::z && a3 == Axis::x) ||
      (a1 == Axis::z && a2 == Axis::x && a3 == Axis::y) ){
    return 1;
  }
  return -1;
}

marqu::Configuration::Configuration(int configuration, int N) : 
  N(N), orientations(new std::pair<Sign, Axis>[N]) {
  set(configuration);
}

marqu::Configuration::Configuration(std::pair<Sign, Axis> * orientations, int N) :
  N(N), orientations(new std::pair<Sign, Axis>[N]) {
  set(orientations);
}

marqu::Configuration::Configuration(const std::string & orientations) : 
  N(orientations.size()/2), orientations(new std::pair<Sign, Axis>[N]) {
  set(orientations);
}

marqu::Configuration::Configuration(const Configuration & other) :
  N(other.N), configuration(other.configuration), 
  orientations(new std::pair<Sign, Axis>[N]){
  for (int i = 0; i < N; ++i) {
    orientations[i] = other.orientations[i];
  }
}

marqu::Configuration::Configuration(Configuration && other) noexcept 
: N(other.N),
  configuration(other.configuration),
  orientations(other.orientations)
{
  other.orientations = nullptr;
}

marqu::Configuration & marqu::Configuration::operator=(const Configuration & other){
  if(this != &other) {
    delete[] orientations;

    N = other.N;
    configuration = other.configuration;
    orientations = new std::pair<Sign, Axis>[N];
    for (int i = 0; i < N; ++i) {
      orientations[i] = other.orientations[i];
    }
  }

  return *this;
}

marqu::Configuration & marqu::Configuration::operator=(Configuration && other) noexcept {
  if (this != &other) {
    delete[] orientations;

    N = other.N;
    orientations = other.orientations;
    configuration = other.configuration;

    other.orientations = nullptr;
  }
  return *this;
}

marqu::Configuration::~Configuration(){
  delete[] orientations;
}

void marqu::Configuration::set(const std::string & orientations){
  for(int i = 0; i < N; ++i){
    this->orientations[i].first = toSign(orientations[2*i]);
    this->orientations[i].second = toAxis(orientations[2*i+1]);
  }
  updateConfiguration();
}

void marqu::Configuration::set(const std::pair<Sign, Axis> * orientations){
  for(int i = 0; i < N; ++i){
    this->orientations[i].first = orientations[i].first;
    this->orientations[i].second = orientations[i].second;
  }
  updateConfiguration();
}

void marqu::Configuration::set(int configuration){
  this->configuration = configuration;
  for (int i = 0; i < N; ++i) {
    int siteConfiguration = configuration % 6;
    orientations[i].first = static_cast<Sign>(siteConfiguration % 2);
    orientations[i].second = static_cast<Axis>(siteConfiguration / 2); 
    configuration /= 6;
  }
}

void marqu::Configuration::subSet(const Configuration & other, 
    const std::vector<std::size_t> & sites){
  for(std::size_t i = 0; i < sites.size(); i++){
    orientations[sites[i]] = other[i];
  }
  updateConfiguration();
}

std::string marqu::Configuration::toString() const{
  std::string res = "";
  for(int i = 0; i < N; ++i){
    res += toChar(this->orientations[i].first);
    res += toChar(this->orientations[i].second);
  }
  return res;
}

int marqu::Configuration::subFlattened(const std::vector<std::size_t> & sites) const{
  int res = 0;
  int multiplier = 1;
  for(int i : sites){
    res += multiplier * toInt(orientations[i].first, orientations[i].second);
    multiplier *= 6;
    //std::cout << "site " << i << " " << res << std::endl;
  }
  return res;
}

void marqu::Configuration::updateConfiguration(){
  configuration = 0;
  int multiplier = 1;
  for(int i = 0; i < N; ++i){
    configuration += multiplier * 
      toInt(orientations[i].first, orientations[i].second);
    multiplier *= 6;
  }
}

int marqu::flattenConfiguration(const std::string & orientations){
  Configuration configuration(orientations);
  return configuration.flattened();
}

