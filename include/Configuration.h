#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <string>
#include <cstdint>
#include <string>

namespace marqu{

  enum class Axis : std::uint8_t {x = 0, y = 1, z = 2};
  enum class Sign : bool {plus = false, minus = true};

  Axis toAxis(char axis);
  char toChar(Axis axis);
  Sign toSign(char sign);
  char toChar(Sign sign);

  int leviCivita(Axis, Axis, Axis);

  template <typename T> T operator*(Sign sign, T value) {
    return (sign == Sign::plus) ? value : -value;
  }
  template <typename T> T operator*(T value, Sign sign) {
    return sign * value; 
  }
  inline int operator*(Sign sign, bool value) {
    return (sign == Sign::plus) ? static_cast<int>(value) : -static_cast<int>(value);
  }
  inline int operator*(bool value, Sign sign) {
    return sign * value;
  }
  inline Sign operator-(Sign sign){
    return (sign == Sign::minus) ? Sign::plus : Sign::minus;
  }

  class Configuration{
    public:
      Configuration(int configuration, int N);
      Configuration(std::pair<Sign, Axis> * orientations, int N);
      Configuration(const std::string & orientations); //More user friendly
      Configuration(const Configuration & other); 
      Configuration(Configuration && other) noexcept; 
      Configuration & operator=(const Configuration & other);
      Configuration & operator=(Configuration&& other) noexcept;
      ~Configuration();

      void set(const std::string & orientations);
      void set(std::pair<Sign, Axis> * orientations); 
      void set(int configuration);

      Sign sign(int site) const {return orientations[site].first;}; 
      Axis axis(int site) const {return orientations[site].second;};
      std::string toString() const;
      int flattened() const {return configuration;}; 
      int getN() const {return N;};

    protected:
      void updateConfiguration();

      std::pair<Sign, Axis> * orientations;
      int configuration; //Flattened configuration index
      int N; //Number of spins
  }; 

  int flattenConfiguration(const std::string & orientations);
}



#endif
