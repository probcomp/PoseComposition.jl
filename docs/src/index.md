# PoseComposition.jl

```@meta
CurrentModule = PoseComposition
```

---

PoseComposition.jl is a library for representing and manipulating poses.  It
places particular emphasis on:

* **Algebraic structure.**  This package implements two different group
  structures on poses.  In the first group structure, composition corresponds
  to change of coordinate frame.  This enables easy manipulation of chains of
  coordinate frames which occur commonly in operations on scene graphs:
  ```math
  \text{pose}_{\text{floor} \to \text{napkin}} = 
  \text{pose}_{\text{floor} \to \text{table}} *
  \text{pose}_{\text{table} \to \text{plate}} *
  \text{pose}_{\text{plate} \to \text{napkin}}
  ```
  In the second group structure, position and orientation are treated
  independently; this models, e.g., inertial motion, where an object has
  separate translational and rotational velocities:
  ```math
  \begin{aligned}
  \mathbf{x} &= \mathbf{x}_0 + \mathbf{v}_{\text{trans}} \Delta t \\
  \boldsymbol{\omega} &= \boldsymbol{\omega}_0 \cdot \mathbf{v}_{\text{rot}}^{\Delta t}
  \end{aligned}
  ```
  Inertial continuous-time dynamics are implemented concisely as
  [`interp`](@ref), which is itself a special case of the [exponential
  map](https://en.wikipedia.org/wiki/Exponential_map_(Lie_theory)) on a Lie
  algebra:
  ```math
  \text{pose}_t = \texttt{interp}(\text{pose}_0,\, \mathbf{v}_{\text{pose}},\, t)
  ```

* **Clear documentation.**  The geometric meaning of the data structures and
  operations is clearly documented, so as to minimize confusion from the many
  ambiguities that typically plague geometry code.

* **Simple equations, ergonomic code.**  Because change of coordinate frame is
  the basic algebraic operation, large equations involving many coordinate
  frames can be solved via simple algebra to express the relative pose between
  any two frames as a function of the others.  For example, in "[Pose Algebra
  In
  Action](https://web.mit.edu/bzinberg/www/GenSceneGraphs-docs/dev/poses/#Example:-Pose-Algebra-In-Action),"
  we solve the equation
  ```
  pose1 * getContactPlane(getShape(g, :obj1), :top)
        * planarContactTo6DOF(pc)
  == pose2 * getContactPlane(getShape(g, :obj2), :curved, 0)
  ```
  for `planarContactTo6DOF(pc)` by simple division in the pose group:
  ```
  planarContactTo6DOF(pc) ==
      (pose1 * getContactPlane(getShape(g, :obj1), :top))
      \ (pose2 * getContactPlane(getShape(g, :obj2), :curved, 0))
  ```
  By contrast, computing and expressing the position and orientation separately
  would require a much larger expression with nested operations and
  harder-to-understand code.


## Quick Start

### Pose composition example

Define `p1`, the pose of a batter (relative to some arbitrary world coordinate
frame) in a game of baseball:
```julia
p1 = Pose([10, 10, 0], UnitQuaternion(1, 0, 0, 0))
```
Define `p1_2`, the relative pose of the baseball in the batter's coordinate
frame, and compute `p2`, the pose of the baseball in the world coordinate
frame:
```julia
p1_2 = Pose([0, 90, 5], RotZYX(0.3, 0.4, 0.5))
p2 = p1 * p1_2
```
Suppose the ball leaves the pitcher's hand at time $t=0$.  If the pitch is a
fastball and there is negligible wind, then the ball's pose at time `t` can be
approximated as inertial:
```julia
# Translational and rotational velocity of the ball relative to the batter
v = Pose([0, -115, -0.2], RotZYX(0, 0, 10))
# Pose of the ball relative to the batter as a function of time
p1_2(t) = p2 * interp(v, t)
```

### Point cloud example

Suppose a robot views the world through a camera that has pose `pose_cam` in
the world coordinate frame, and perceives a point cloud (represented as a ``3
\times N`` matrix) `ptcloud_cam` in the camera's coordinate frame.  Then the
same point cloud expressed in the world coordinate frame is
```julia
ptcloud_world = pose_cam * ptcloud_cam
```
Conversely, if we know of a point cloud `ptcloud_world` expressed in the world
coordinate frame and want to re-express it in the camera's coordinate frame:
```julia
ptcloud_cam = pose_cam \ ptcloud_world
```
