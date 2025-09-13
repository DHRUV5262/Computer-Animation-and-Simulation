# Real-Time Graphics Playground

A personal exploration of **advanced computer graphics** techniques—including soft-body physics and motion-capture animation—implemented in modern C++ and OpenGL.

---

## ✨ Projects

### 1. Physics-Based *Jello Cube* Simulator
Deforms and jiggles in real time using a **mass-spring lattice** with structural, shear and bend springs.

* Stable numerical integration (explicit Euler, semi-implicit Euler, RK4)
* Self-contained collision detection & response (planes + friction)
* Adjustable elasticity, damping and gravity in the GUI
* Wireframe / shaded rendering toggle for debugging
* FPS-independent time-stepping for consistent behaviour

### 2. Motion-Capture Interpolation Engine
Loads **ASF / AMC** skeletons and generates new motion by blending sparse keyframes.

* Supports Euler, quaternion and cubic Bézier interpolation
* Interactive scrubber + playback controls at up to 120 fps
* Side-by-side comparison view for different interpolation modes
* Automatic keyframe extraction with configurable interval
* Screenshot & video recording utilities

---

## 🚀 Quick Start

clone
git clone 
cd 

build (Linux / macOS)
cd build && cmake .. && make -j

run demos
./bin/JelloCubeSim
./bin/MocapPlayer


Windows users: open the provided **Visual Studio 2022** solution and build the `ALL_BUILD` target.

---

## 🛠️ Tech Stack
| Layer            | Tools / Libraries           |
|------------------|-----------------------------|
| Language         | C++17                       |
| Graphics         | OpenGL 4.6, GLSL            |
| GUI              | FLTK 1.3                    |
| Math / Helpers   | GLM, Eigen                  |
| Build System     | CMake + Make / Visual Studio|

---

## 📊 Results
* **Soft-body demo:** elastic deformation, realistic damping, stable under large time steps.
* **Mocap demo:** quaternion Slerp produces visibly smoother rotations than Euler Lerp; Bézier curves give the most natural ease-in/ease-out.

---

## ✍️ What I Learned
* Designing numerically stable mass-spring systems
* Implementing multiple integration schemes and comparing accuracy vs. performance
* Parsing and rendering hierarchical skeletons
* Mastering quaternion algebra for rotation interpolation
* Building cross-platform OpenGL applications with CMake


