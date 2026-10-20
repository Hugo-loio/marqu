#include "TreapCappedParticleSimulator2.h"
#include <stdexcept>

marqu::Particle marqu::TreapCappedParticleSimulator2::getParticleByIndex(Node *tree, int count) {
  if (!tree) {
    throw std::runtime_error("Treap index out of bounds");
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

void marqu::TreapCappedParticleSimulator2::deleteRandomParticle() {
  if (particleNumber <= 0 || !root) return;
  std::uniform_int_distribution<int> dist(0, particleNumber - 1);
  int idx = dist(gen);
  Particle p = getParticleByIndex(root, idx);
  removeParticle(p);
}

void marqu::TreapCappedParticleSimulator2::markovStep(const Particle &particle) {
  TreapParticleSimulator::markovStep(particle);
  while (particleNumber > maxParticleNumber) {
    deleteRandomParticle();
  }
}
