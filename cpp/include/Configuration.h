#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <string>
#include <vector>
#include <cstdint>

namespace marqu{

  enum class Axis : std::uint8_t {x = 0, y = 1, z = 2};
  enum class Sign : bool {plus = false, minus = true};

  Axis toAxis(char axis);
  char toChar(Axis axis);
  Sign toSign(char sign);
  char toChar(Sign sign);
  int toInt(Sign sign, Axis axis);
  int toInt(const std::pair<Sign, Axis> & orientation);

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
  inline Sign operator*(Sign a, Sign b) {
    return static_cast<Sign>(static_cast<bool>(a) ^ static_cast<bool>(b));
  }
  inline Sign operator-(Sign sign){
    return static_cast<Sign>(!static_cast<bool>(sign));
    //return (sign == Sign::minus) ? Sign::plus : Sign::minus;
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

      bool operator==(const Configuration & other) const;
      bool operator<(const Configuration & other) const;
      bool operator<=(const Configuration & other) const;
      const std::pair<Sign, Axis> & operator[](std::size_t index) const{
	return orientations[index];};

      void set(const std::string & orientations);
      void set(const std::pair<Sign, Axis> * orientations); 
      void set(int configuration);
      void subSet(const Configuration & other, const std::vector<std::size_t> & sites);

      Sign sign(std::size_t site) const {return orientations[site].first;}; 
      Axis axis(std::size_t site) const {return orientations[site].second;};
      std::string toString() const;
      int flattened() const {return configuration;}; 
      int subFlattened(const std::vector<std::size_t> & sites) const;
      int getN() const {return N;};

    protected:
      void updateConfiguration();

      int N; //Number of spins
      int configuration; //Flattened configuration index, might overflow!
      std::pair<Sign, Axis> * orientations;
  }; 

  int flattenConfiguration(const std::string & orientations);
}



#endif
