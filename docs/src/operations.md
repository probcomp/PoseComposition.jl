# Operations on Poses

```@meta
CurrentModule = PoseComposition
```

---

A [`Pose`](@ref) can represent either a relative pose (relative to some parent
coordinate frame which must be supplied by context) or an absolute pose.  An
absolute pose is, by definition, a pose relative to the world coordinate frame.

## First group structure: change of coordinate frame

A central operation on poses is change of coordinate frame:

```
(absolute pose of frame2) = (absolute pose of frame1) * (relative pose of frame2 relative to frame1)
```

This operation gives the set of poses a group structure -- call it the "pose
group" ``G``.  This ``G`` is the
[opposite group](https://en.wikipedia.org/wiki/Opposite_group) of the group of
[rigid transformations](https://en.wikipedia.org/wiki/Special_Euclidean_group)

``SE_3 = \mathbb{R}^3 \ltimes SO_3(\mathbb{R})``

(technical note[^1]), meaning that relative poses are "applied" on the right,
rather than written as function applications on the left.

[^1]:
    In code, elements of ``G`` whose orientation is represented as a
    `Rotations.UnitQuaternion` (as opposed to some other subtype of
    `Rotations.Rotation{3}`) remember their sign (`q` and `-q` are considered
    equal by `Base.:==`, but their fields are not equal).

The group operation is implemented via `Pose`-specific methods of the standard
Julia functions in `Base`:

* Multiplication: `pose1 * pose2`
* Inversion: `inv(pose)`, such that
  `pose * inv(pose) == inv(pose) * pose == IDENTITY_POSE`
* Right-division: `pose1 / pose2 := pose1 * inv(pose2)`
* Left-division: `pose1 \ pose2 := inv(pose1) * pose2`

#### Convention: First translate, then rotate

Note that we define the group operation on the pose group by doing the
translation first, then the rotation.  So the coordinate frame `pose1 * pose2`
has an origin whose coordinates in the frame of `pose1` are `pose2.pos` (the
translational component of `pose2`).

Note also that the orientation `orn := pose.orientation` represents the linear
operator

    v ↦ Rotations.RotMatrix{3}(orn) * v   =:   orn * v

In other words, we are using the convention that matrix-vector multiplication
has the vector on the right.  The x, y and z axes of an object with pose `pose`
are `pose.orientation * [1, 0, 0]`, `pose.orientation * [0, 1, 0]` and
`pose.orientation * [0, 0, 1]` respectively.

##### Example: Which order do I multiply things in?

Suppose we start with coordinate frame `p1`, and we want to translate its
origin by `p2.pos` and then rotate its coordinate axes via the linear map
`v ↦ p2.orientation * v`.  Would the pose corresponding to the new coordinate
frame be written as `p1 * p2` or `p2 * p1`?

Answer: `p1 * p2`.

#### Mathematical details: contravariance, coordinate frames, and rigid motions

There are two ways to think about poses:

1. A pose is a coordinate frame, described with respect to some other base
   coordinate frame.  The `pos` says where the frame's origin is, and the
   `orientation` says what directions its right-handed orthonormal axes point
   in.
2. A pose is a rigid motion.  The translational component is `pos` and the
   orientation component is `orientation`.

When we talk about the [group
action](https://en.wikipedia.org/wiki/Group_action_(mathematics)) of ``G`` on
itself by multiplication, we want to say that a rigid transformation (way 2
above) acts on a coordinate frame (way 1 above), by moving the origin a
translational offset of `pos` and matrix-multiplying its rotation by
`orientation`.  The thing to note is that, even though function application is
usually written on the left, this group action *right*-associative, so it must
be represented by right-multiplication, not left-multiplication (or else the
associative law breaks).

### Points and point clouds

Above we explained how a pose can act on another pose.  It is also possible for
a pose to act on a point in ``\mathbb{R}^3``.  To see how, note that for poses
`a` and `b`, the position of the product, `(a * b).pos`, does not depend on
`b.orientation`.  Therefore it is valid to define, for any `a::Pose` and any
``\texttt{bpos} ∈ \mathbb{R}^3``,

    a * bpos := (a * Pose(bpos, orn)).pos

where `orn` is any orientation (and its value does not affect the result).

In words, `a * bpos` is the coordinates in the world frame of the vector whose
coordinates in the frame represented by `a` are `bpos`.

From this it's straightforward to see that the usual associativity law holds:

    a * (b * cpos) = (a * b) * cpos

Therefore pose–point multiplication is indeed a group action of ``G`` on
``\mathbb{R}^3``.

This package implements both `Pose * point` and a vectorized version, `Pose *
pointcloud`, where `pointcloud` is a ``3 × N`` matrix in which each column
represents a point:

```julia-repl
julia> Pose([1, 2, 3], RotZYX(0.4, 0.5, 0.6)) * @SVector([7, 8, 9])

3-element StaticArrays.SArray{Tuple{3},Float64,1,3} with indices SOneTo(3):
 11.340627918648092
  8.023198113213361
 10.126885626769848

julia> Pose([1, 2, 3], RotZYX(0.4, 0.5, 0.6)) * [7  10  13  16;
                                                 8  11  14  17;
                                                 9  12  15  18]

3×4 Array{Float64,2}:
 11.3406  15.3024  19.2641  23.2258
  8.0232  10.5473  13.0714  15.5955
 10.1269  12.3481  14.5693  16.7904
```


## Second group structure: direct product of position and orientation

The composition operation described above and denoted by `*` is useful when
manipulating chains of relative poses, which commonly occur in scene graphs.
We now define a second group structure on the set of `Pose`s: the direct
product ``\mathbb{R}^3 × SO_3(\mathbb{R})``.  We denote the group operation in
this second structure by `⊗`, and it means "add the positions and (separately)
compose the orientations."  This operation models inertial motion, where an
object has separate translational and rotational velocities:

```math
\begin{aligned}
\mathbf{x} &= \mathbf{x}_0 + Δ\mathbf{x} \\
\boldsymbol{ω} &= \boldsymbol{ω}_0 \cdot Δ\boldsymbol{ω}
\end{aligned}
```

or simply

```julia
pose1 = pose0 ⊗ Δpose
```

This second group operation is implemented via the following infix operators:

* Multiplication: `pose1 ⊗ pose2` (``\LaTeX``: `\otimes`)
* Right-division: `pose1 ⊘ pose2` (``\LaTeX``: `\oslash`)
* Left-division: `pose1 ⦸ pose2` (``\LaTeX``: `\obslash`)

The direct product structure is useful for performing computations that treat
translational and rotational components separately, such as computing the pose
of a moving object by integrating its velocity.  The translational and
rotational components of the velocity are integrated separately -- in
differential geometry terms, motion is governed by the [exponential
map](https://en.wikipedia.org/wiki/Exponential_map_(Lie_theory)) on the Lie
algebra of ``\mathbb{R}^3 × SO_3(\mathbb{R})``, not the Lie algebra of
``SE_3``.

### Inter/Extrapolation, a.k.a. stepping time forward a fractional number of times

One nice property of the direct product (and its physical application of
integrating a velocity) is that there is a straightforward analogue of the
physics 101 equation

```math
\mathbf{x} = \mathbf{x}_0 + t \mathbf{v}
```

where ``\mathbf{x}_0`` is the initial displacement of an object, ``t`` is time,
and ``\mathbf{v}`` is the object's velocity, which is assumed to stay constant.
That analogue is

    pose = pose0 ⊗ interp(v, t)

where `pose`, `pose0` and `v` are `Pose`s, and `t` is a `Real`.  The function
[`interp`](@ref) defines how to interpolate between doing nothing and
translating by `v.pos` (namely, scalar multiplication) and how to interpolate
between doing nothing and applying the rotation `v.orientation` (namely,
[SLERP](https://en.wikipedia.org/wiki/Slerp#Quaternion_Slerp)).

There is also a variant of `interp` that takes two pose arguments and
interpolates between them:

    pose = interp(pose0, pose1, t)

which is equivalent to

    pose = pose0 ⊗ interp(pose0 ⦸ pose1, t)

In particular, we have `pose = pose0` when `t = 0` and `pose = pose1` when `t =
1`.
