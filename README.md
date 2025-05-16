# GPS Correlator Simulation
This repo is intended to serve as a useful implementation of the algorithms described in a paper titled "Derivation and Validation of a Higher-Fidelity GPS Correlator Model" from the ION/IEEE PLANS 2025 conference proceedings. Hopefully, this can either be a "black box" that readers may easily use, or a means to validate the behavior of custom implementations.

## CPP Implementation
For now, only a c++ implementation is provided, compatible with C++11 or higher. It is header-only, so it may easily be included in your project. The only dependency is Eigen. 

### Using with CMake
If you wish to use gcs in your CMake project, simply add the include directory to your desired target. As an example, if you cloned this repo such that it is contained in a subdirectory called `gcs`, then add the following to your CMakeLists file, altering `target_name_here` to your desired target:
```
target_include_directories(target_name_here ${CMAKE_CURRENT_SOURCE_DIR}/gcs/cpp/include)
```
Then, in your code, add the include:
```
#include <gcs/correlator_sim.hpp>
```

### Example Usage
TODO

## Contact
For any questions on the contents of the paper or this repository, email bzb0120@auburn.edu
