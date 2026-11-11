#include "TreapParticleSimulatorKill.h"
#include <stdexcept>

marqu::TreapParticleSimulatorKill::TreapParticleSimulatorKill(double rate)
    : killRate(rate / 2) {
  totalRate = (root ? root->cumulRate : 0.0) +
              (particleNumber > 2 ? killRate * particleNumber : 0.0);
}

void marqu::TreapParticleSimulatorKill::setKillRate(double rate) {
  killRate = rate / 2;
  totalRate = (root ? root->cumulRate : 0.0) +
              (particleNumber > 2 ? killRate * particleNumber : 0.0);
}

marqu::Particle
marqu::TreapParticleSimulatorKill::getParticleByIndex(Node *tree, int count) {
  if (!tree) {
    throw std::runtime_error("Treap particle index out of bounds");
  }

  int nLeft = getCumulParticles(tree->left);
  if (count < nLeft) {
    return getParticleByIndex(tree->left, count);
  }

  int nNode = tree->nParticles + tree->nAntiParticles;
  if (count < nLeft + nNode) {
    count -= nLeft;
    bool type = (count < tree->nParticles);
    return Particle(tree->configuration, type, tree->rate);
  }

  return getParticleByIndex(tree->right, count - nLeft - nNode);
}

void marqu::TreapParticleSimulatorKill::deleteRandomParticle() {
  std::uniform_int_distribution<int> dist(0, particleNumber - 1);
  int idx = dist(gen);
  Particle p = getParticleByIndex(root, idx);
  removeParticle(p);
}

double marqu::TreapParticleSimulatorKill::gillespieTimeStep() {
  if (totalRate <= 0.0) {
    throw std::runtime_error("All particles are stationary");
  }

  std::exponential_distribution<> exp_dist(totalRate);
  double dt = exp_dist(gen);

  double randRate = uni_dist(gen) * totalRate;

  if (randRate < root->cumulRate) {
    TreapParticleSimulator::gillespieTimeStep(root, randRate);
  } else {
    deleteRandomParticle();
    deleteRandomParticle();
  }

  totalRate = (root ? root->cumulRate : 0.0) +
              (particleNumber > 2 ? killRate * particleNumber : 0.0);
  return dt;
}
