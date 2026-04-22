# MarQu

A library for simulating the Schrödinger equation of spin chains based on a classical stochastic interpretation as a continuous-time Markov chain. 

Details of the underlying mathematical formulation and algorithm implementation can be found in:

**TBA**

If you use this code in your research please cite the work above.


## Features

MarQu splits the problem of quantum simulation into two distinct steps: a Python-based preprocessing phase for mathematical optimization, and a high-performance C++ engine for the Monte Carlo execution.

* **Rate Matrix Generation (Python):** Construct transition rate matrices from arbitrary local Hamiltonian and Lindbladian/noise terms.
* **Optimal Gauge Transformation (Python):** Utilizes Linear Programming (`scipy.optimize.linprog`) to find the optimal gauge transformation that minimizes the sum of negative rates, effectively mitigating the sign problem and particle proliferation.
* **High-Performance Simulation Engine (C++):** 
    * `TreapParticleSimulator`: Handles "negative" Markov chains by simulating particles and antiparticles.
Uses a dynamic Treap (Tree + Heap) data structure to maintain efficient O(log N) insertion, deletion, and annihilation operations.
    * `ClassicalParticleSimulator`: Optimized simulator for regimes where sufficient noise has driven the system into a purely classical state (all positive rates).
* **Initial State Sampling:** Built-in support for sampling Product States and volume-law entangled non-local Bell Pair states.

## User guide

Install with `./install.sh` and uninstall with `./unistall.sh`.

Examples on how to use the package provided in the `examples` directory.

## Architecture Diagrams

To help you understand the high-level design of MarQu, we provide class diagrams in the `diagrams/` directory written in [PlantUML](https://plantuml.com/), which outline the primary user-facing APIs for both the CPP and Python packages.
