#ifndef TREAPPARTICLESIMULATOR_H
#define TREAPPARTICLESIMULATOR_H

#include "BaseParticleSimulator.h"
#include "Particle.h"
#include <cstdint>

namespace marqu{
  class TreapParticleSimulator : public BaseParticleSimulator{
    public:
      TreapParticleSimulator() : int_dist(1, INT32_MAX) {};
      ~TreapParticleSimulator() {clear(root);};

      void discreteTimeStep(double dt);
      double gillespieTimeStep();

      void displayParticles() const;
      double compressionRate() {return (double)root->cumulParticles/root->size;};

      bool annihilateParticles = true;

    protected:
      void clearParticles();
      void addParticle(Particle && particle);
      void addParticle(const Particle & particle);
      //Not in the base class
      void removeParticle(const Particle & particle);
      void moveParticle(const Particle & out, const Particle & in);
      void moveParticle(const Particle & out, Particle && in);
      virtual void markovStep(const Particle & particle); 

      std::uniform_int_distribution<int> int_dist;
      bool warnDiscrete = true;

      // Nested treap type
      class Node {
	public:
	  Node(const Particle & particle, int prio);
	  Node(Particle && particle, int prio) noexcept;

	  //The tree is sorted is ascending order of configuration flat indices
	  Configuration configuration;
	  int nParticles;
	  int nAntiParticles;
	  double rate; // Assume event rate depends only on the configuration
	  int size = 1; // subtree size
	  double cumulRate; // subtree summed event rates
	  int cumulParticles = 1; // subtree particle number
	  int priority;
	  Node *left, *right;
      };

      // Functions for treap manipulation
      int getSize(Node * node){return node ? node->size : 0;}
      int getCumulParticles(Node * node){return node ? node->cumulParticles : 0;}
      double getCumulRate(Node * node){return node ? node->cumulRate : 0.0;}
      void recalc(Node * node); 
      void split(Node * tree, const Configuration & config, 
	  Node *& left, Node *& right); 
      Node * merge(Node * left, Node * right);
      Node * findConfig(Node * tree, const Configuration & config);
      void remove(Node *& tree, const Configuration & config);
      // Returns true if the node gets cleared
      bool remove(Node * tree, const Particle &); 
      // Returns 0: added particle, 1: config not found, 2: cleared node 
      int add(Node * tree, const Particle &); 
      void clear(Node *& tree);
      void displayParticles(Node * tree) const;
      void discreteTimeStep(Node * tree, double dt, int count);
      void gillespieTimeStep(Node * tree, double cumulRate);
      Node * root = nullptr;
  };
}

#endif
