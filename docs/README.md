# libcloudph++ documentation

🌦️ Welcome to the libcloudph++ documentation. This is a work in progress.





### Microphysical Schemes


####  💧 **Lagrangian Scheme** (`src/`)
- Particle-based microphysics using Super-Droplet Method, based on Shima et al. (2009) and Shima et al. (2020).
- **Available processes**:
  - Growth of cloud droplets with the Maxwell-Mason equation
  - Collision-coalescence
  - Advection and sedimentation
  - Ice processes: work in progress
- **Type**: Header-only library

####  ☁️ **1-Moment Bulk Scheme** (`include/libcloudph++/blk_1m/`)
- Single-moment bulk scheme based on Kessler (1995), with ice microphysics parametrization based on Grabowski (1999).
Only the total mass of water per category (cloud / rain / iceA / iceB) is considered.
  Transport equations for the
  cloud water mixing ratio rc, the rain water mixing ratio
  rr, ice A mixing ratio ria and ice B mixing ration rib are solved in addition to the state variables θ and rv representing heat and moisture content.
It is a simplistic approach. It does not contain information about the shape of particle size
  distribution.
- **Available processes**:
    - Condensation and evaporation
    - Autoconversion and collection 
    - Sedimentation
    - Ice processes (A and B types), including nucleation, growth by deposition and riming
- **Type**: Header-only library

####  🌧️ **2-Moment Bulk Scheme** (`include/libcloudph++/blk_2m/`)
- The double-moment scheme implemented in libcloudph++ was introduced by Morrison and Grabowski
  (2007). Similarly to the singlemoment approach, the double-moment warm-rain scheme
  assumes that condensed water is divided into two categories:
  cloud water and rain water. In addition to the total mass of
  water in both categories, concentrations of droplets and drops
  are also predicted. As a result, the scheme considers two moments of particle size distribution. In the Eulerian framework, four transport equations for cloud droplet
  concentration, cloud water mixing ratio, rain drop concentration and rain water mixing ratio are solved.
- **Available processes**:
    - Condensation and evaporation
    - Autoconversion and collection
    - Sedimentation
- **Type**: Header-only library


<br><br><br>



###  Common Components

#### **1. Shared Utilities** (`include/libcloudph++/common/`)
- **Thermodynamics**: Saturation vapor pressure, latent heats
- **Terminal velocity**: Particle fall speed calculations
- **Physical constants**: Atmospheric and water properties
- **Mathematical utilities**: Numerical methods and helpers



#### 2. Build System

The project uses **CMake** as its build system:
- `CMakeLists.txt`: Main build configuration
- `cmake/`: Additional CMake modules and utilities
- Header-only schemes require no compilation
- Lagrangian scheme produces linkable libraries

#### 3. Language Bindings

The `bindings/` directory provides interfaces for:
- **Python** NumPy 
- **Fortran**

#### 4. Testing

The `tests/` directory contains:
- **Unit tests**: individual component validation
- **Integration tests**: full scheme testing
- **Benchmark cases**: performance and accuracy validation
- **Inter-scheme comparisons**:

---
