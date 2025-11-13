#include "TreapParticleSimulator.h"
#include <iostream>
#include <algorithm>


void marqu::TreapParticleSimulator::discreteTimeStep(double dt){
  if(warnDiscrete){
    std::cout << "Warning: using discrete timesteps in the TreapParticleSimulator might be ineficient. The time step complexity is O(n' log n') as opposed to O(n) for VectorParticleSimulator, for n particles distributed by n' configurations" << std::endl;
    warnDiscrete = false;
    //TODO: This can possibly be fixed by copying the tree into a vector, 
    //looping through the vector but acting on the tree and deleting the vector
    //at the end
  }
  for(int i = 0; i < particleNumber; i++){
    discreteTimeStep(root, dt, i); 
  }
}

double marqu::TreapParticleSimulator::gillespieTimeStep(){
  //std::cout << "Gillespie: " << particleNumber << "
  if(totalRate <= 0) {
    throw std::runtime_error("All particles are stationary");
  }
  std::exponential_distribution<> exp_dist(totalRate);
  double dt = exp_dist(gen);
  double randRate = uni_dist(gen)*totalRate;
  gillespieTimeStep(root, randRate);
  return dt;
}


void marqu::TreapParticleSimulator::clearParticles(){
  clear(root);
  trackClear();
}
void marqu::TreapParticleSimulator::addParticle(const Particle & particle){
  //std::cout << "Adding : " << particle.type << " " << particle.configuration.toString() << std::endl;
  int res = add(root, particle);
  //std::cout << "Res: " << res << std::endl;
  if(res == 1){
    Node *left, *right;
    split(root, particle.configuration.flattened(), left, right);
    root = merge(merge(left, new Node(particle, int_dist(gen))), right);
  }
  else if(res == 2) remove(root, particle.configuration.flattened());
  trackAdd(particle); 
}

void marqu::TreapParticleSimulator::addParticle(Particle && particle){
  int res = add(root, particle);
  if(res == 1){
    Node *left, *right;
    split(root, particle.configuration.flattened(), left, right);
    root = merge(merge(left, new Node(std::move(particle), int_dist(gen))), right);
  }
  else if(res == 2) remove(root, particle.configuration.flattened());
  trackAdd(particle);
}

void marqu::TreapParticleSimulator::removeParticle(const Particle & particle){
  if(remove(root, particle)) remove(root, particle.configuration.flattened());
  trackRemove(particle);
}

void marqu::TreapParticleSimulator::moveParticle(const Particle & out, const Particle & in){
  removeParticle(out);
  addParticle(in);
}

void marqu::TreapParticleSimulator::moveParticle(const Particle & out, Particle && in){
  removeParticle(out);
  addParticle(std::move(in));
}

// Simplified for unitary evolution (TODO: incorporate open systems later)
void marqu::TreapParticleSimulator::markovStep(const Particle & particle){
  std::pair<Configuration, Sign> res = randomEvent(particle);
  bool type = (res.second == Sign::plus) ? particle.type : ! particle.type;
  Node * config = findConfig(root, res.first.flattened());
  double rate = config ? config->rate : eventRate(res.first);

  Particle newParticle(std::move(res.first), type, rate);
  addParticle(newParticle);
}

marqu::TreapParticleSimulator::Node::Node(const Particle & particle, int prio):
  configuration(particle.configuration), cumulRate(particle.eventRate), 
  rate(particle.eventRate), priority(prio), left(nullptr), right(nullptr),
  nParticles(particle.type ? 1 : 0), nAntiParticles(particle.type ? 0 : 1)
{
}

marqu::TreapParticleSimulator::Node::Node(Particle && particle, int prio) noexcept:
configuration(std::move(particle.configuration)), priority(prio),
  rate(particle.eventRate), cumulRate(particle.eventRate), 
  left(nullptr), right(nullptr),
  nParticles(particle.type ? 1 : 0), nAntiParticles(particle.type ? 0 : 1)
{
}

void marqu::TreapParticleSimulator::recalc(Node* node){
  if (!node) return;
  node->cumulRate = (node->nParticles + node->nAntiParticles) * node->rate + 
    getCumulRate(node->left) + getCumulRate(node->right);
  node->cumulParticles = node->nParticles + node->nAntiParticles +
    getCumulParticles(node->left) + getCumulParticles(node->right);
  node->size = 1 + getSize(node->left) + getSize(node->right);
}

// Node with config goes to the right tree
void marqu::TreapParticleSimulator::split(Node * tree, int config, Node *& left, Node *& right) {
  if(!tree) return void(left = right = nullptr);
  if(config <= tree->configuration.flattened()){
    split(tree->left, config, left, tree->left);
    right = tree;
  }
  else{
    split(tree->right, config, tree->right, right);
    left = tree;
  }
  recalc(tree);
}

marqu::TreapParticleSimulator::Node * marqu::TreapParticleSimulator::merge(Node * left, Node * right){
  if(!left || !right) return left ? left : right;
  if(left->priority > right->priority){
    left->right = merge(left->right, right);
    recalc(left);
    return left;
  } 
  else{
    right->left = merge(left, right->left);
    recalc(right);
    return right;
  }
}

marqu::TreapParticleSimulator::Node * marqu::TreapParticleSimulator::findConfig(Node * tree, int config){
  if(!tree) return nullptr;
  if(tree->configuration.flattened() == config) return tree;
  else if(config < tree->configuration.flattened()){
    return findConfig(tree->left, config);
  }
  else{
    return findConfig(tree->right, config);
  }
}

void marqu::TreapParticleSimulator::remove(Node *& tree, int config){
  Node * left, * rm, * right;
  split(tree, config, left, right);
  split(right, config + 1, rm, right); 
  delete rm;
  tree = merge(left, right);
}

bool marqu::TreapParticleSimulator::remove(Node * tree, const Particle & particle){
  bool emptyNode = false;
  if(!tree) throw std::runtime_error("Removing a particle that doesn't exist");
  if(tree->configuration.flattened() == particle.configuration.flattened()){
    particle.type ? tree->nParticles-- : tree->nAntiParticles--;
    if(tree->nParticles == 0 && tree->nAntiParticles == 0) emptyNode = true;
  }
  else if(particle.configuration.flattened() < tree->configuration.flattened()){
    emptyNode = remove(tree->left, particle);
  }
  else{
    emptyNode = remove(tree->right, particle);
  }
  recalc(tree);
  return emptyNode;
}

int marqu::TreapParticleSimulator::add(Node * tree, const Particle & particle){
  int returnVal = 0;
  if(!tree) return 1;
  if(tree->configuration.flattened() == particle.configuration.flattened()){
    particle.type ? tree->nParticles++ : tree->nAntiParticles++;
    //std::cout << "Adding: " << particle.type << " " << particle.configuration.toString() << std::endl;
    //std::cout << "Check: " << tree->configuration.flattened() << " " << particle.configuration.flattened() << " " << tree->configuration.toString() << " " << particle.configuration.toString() << std::endl;
    if(annihilateParticles){
      int nAnnihilate = std::min(tree->nParticles, tree->nAntiParticles);
      //std::cout << "anni: " << tree->configuration.toString() << " " << tree->nParticles << " " << tree->nAntiParticles << std::endl;
      if(nAnnihilate > 0){
	tree->nParticles -= nAnnihilate;
	tree->nAntiParticles -= nAnnihilate;
	trackAnnihilate(nAnnihilate, tree->rate);
	if(tree->nParticles == tree->nAntiParticles) returnVal = 2;
	//std::cout << "Changed: " << returnVal << std::endl;
      }
    }
  }
  else if(particle.configuration.flattened() < tree->configuration.flattened()){
    returnVal = add(tree->left, particle);
  }
  else{
    returnVal = add(tree->right, particle);
  }
  recalc(tree);
  //std::cout << "Returning: " << returnVal << std::endl;
  return returnVal;
}

void marqu::TreapParticleSimulator::clear(Node *& tree){
  if (!tree) return;
  clear(tree->left);
  clear(tree->right);
  delete tree;
  tree = nullptr;  
}

void marqu::TreapParticleSimulator::displayParticles(Node * tree) const{
  if(!tree) return;
  displayParticles(tree->left);
  std::cout << "Configuration " << tree->configuration.toString() << ": " 
    << tree->nParticles << " particles, " << tree->nAntiParticles 
    << " antiparticles" << std::endl;
  displayParticles(tree->right);
}

void marqu::TreapParticleSimulator::discreteTimeStep(Node * tree, double dt, int count){
  if(!tree) return;

  int nLeft = getCumulParticles(tree->left);
  int nNode = tree->nParticles + tree->nAntiParticles;
  if(count <= nLeft) return discreteTimeStep(tree->left, dt, count);
  else if(nLeft + nNode < count) 
    return discreteTimeStep(tree->right, dt, count - nLeft - nNode);

  count -= nLeft;
  Particle particle(tree->configuration, true, tree->rate);
  if(tree->nParticles < count) particle.type = false;
  if(uni_dist(gen) < tree->rate*dt) markovStep(particle);
}

void marqu::TreapParticleSimulator::gillespieTimeStep(Node * tree, double cumulRate){
  if (!tree) return;

  double rateLeft = getCumulRate(tree->left);
  double rateNode = (tree->nParticles + tree->nAntiParticles) * tree->rate;
  if (cumulRate < rateLeft) return gillespieTimeStep(tree->left, cumulRate);
  else if(rateLeft + rateNode < cumulRate) 
    return gillespieTimeStep(tree->right, cumulRate - rateLeft - rateNode);

  cumulRate -= rateLeft;
  Particle particle(tree->configuration, true, tree->rate);
  if(tree->nParticles * tree->rate < cumulRate) particle.type = false;
  markovStep(particle);
}
