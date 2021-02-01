module PoseComposition

import Base: @kwdef
import Rotations: Rotation, UnitQuaternion, RotZYX
import StaticArrays: StaticVector, SVector, @SVector

include("docstring_extensions.jl")

@doc """Fields:
$(TYPEDFIELDS)$(raw"""

---

Struct representing the pose (position and orientation) of an object.

Can represent either a relative pose (relative to some parent coordinate frame
which must be supplied by context) or an absolute pose.  An absolute pose is,
by definition, a pose relative to the world coordinate frame.

# First Group Structure: Composition of Poses

There is a natural way to combine two poses, which gives the set of poses a
group structure.  This group (call it the "pose group" ``G =
SE_3^{\text{op}}``) is the [opposite
group](https://en.wikipedia.org/wiki/Opposite_group) of the group of rigid
transformations

``SE_3 = \mathbb{R}^3 \ltimes SO_3(\mathbb{R})``

except that in code, elements of ``G`` whose orientation is represented as a
`Rotations.UnitQuaternion` (as opposed to some other subtype of
`Rotations.Rotation{3}`) remember their sign (`q` and `-q` are considered equal
by `Base.:==`, but their fields are not equal).

This group operation is implemented via `Pose`-specific methods of the standard
Julia functions in `Base`:

* Multiplication: `pose1 * pose2`
* Inversion: `inv(pose)`, such that
  `pose * inv(pose) == inv(pose) * pose == IDENTITY_POSE`
* Right-division: `pose1 / pose2 := pose1 * inv(pose2)`
* Left-division: `pose1 \ pose2 := inv(pose1) * pose2`

The right action ``G`` on itself (by right-multiplication) is equivalent to the
left action of ``SE_3`` on itself (by left-composition).

## First translate, then rotate

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

## Contravariance: Left action of ``SE_3`` = Right action of ``G``

### TLDR

"Which order do I multiply things in?"

Say we want to start with pose `p1` and translate its origin by `p2` and
matrix-multiply its orientation by `p2`.  Would that new pose be written as
`p1 * p2` or `p2 * p1`?

Answer: `p1 * p2`.

## Mathematical details

There are two ways to think of poses:

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

## `Pose * point`

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


# Second group structure: Direct Product of Position and Orientation

The composition operation described above and denoted by `*` is useful when
manipulating scene graphs.  We now define a second group structure on the set
of `Pose`s: the direct product ``\mathbb{R}^3 × SO_3(\mathbb{R})``.  We denote
the group operation in this second structure by `⊗`, and it means "add the
positions and compose the orientations."

This second group operation is implemented via infix operators defined in the
`GenSceneGraphs` module:

* Multiplication: `pose1 ⊗ pose2` (``\LaTeX``: `\otimes`)
* Right-division: `pose1 ⊘ pose2` (``\LaTeX``: `\oslash`)
* Left-division: `pose1 ⦸ pose2` (``\LaTeX``: `\obslash`)

The direct product structure is useful for performing computations that treat
translational and rotational components separately, such as computing the pose
of a moving object by integrating its velocity.  The translational and
rotational components of the velocity are integrated separately -- in
differential geometry terms, motion is governed by the exponential map on the
Lie algebra of ``\mathbb{R}^3 \times SO_3(\mathbb{R})``, not the Lie algebra of
``SE_3``.

## Inter/Extrapolation, a.k.a. stepping time forward a fractional number of times

One nice property of the direct product (and its physical application of
integrating a velocity) is that there is a straightforward analogue of the
physics 101 equation

`` \vec{x} = \vec{x}_0 + t \vec{v} ``

where ``\vec{x}_0`` is the initial displacement of an object, ``t`` is time,
and ``\vec{v}`` is the object's velocity, which is assumed to stay constant.
That analogue is

    pose = pose0 ⊗ interp(v, t)

where `pose`, `pose0` and `v` are `Pose`s, and `t` is a `Real`.  The function
[`interp`](@ref) defines how to interpolate between doing nothing and
translating by `v.pos` (namely, scalar multiplication) and how to interpolate
between doing nothing and applying the rotation `v.orientation` (namely,
[SLERP](https://en.wikipedia.org/wiki/Slerp#Quaternion_Slerp)).
""")""" Pose
@kwdef struct Pose
  """Origin of the object's coordinate frame."""
  pos::StaticVector{3, <:Real}
  """Orientation of the object's coordinate frame."""
  orientation::Rotation{3}
end


### Constructors ###

function Pose(x::Real, y::Real, z::Real, orientation::Rotation{3})::Pose
  Pose(@SVector([x, y, z]), orientation)
end

function Pose(pos::AbstractVector{<:Real}, orientation::Rotation{3})::Pose
  (x, y, z) = pos
  Pose(@SVector([x, y, z]), orientation)
end

function Pose(pos::AbstractVector{<:Real}, ypr::NamedTuple{(:yaw, :pitch, :roll)})::Pose
  Pose(pos, RotZYX(ypr.yaw, ypr.pitch, ypr.roll))
end

function Pose(x::Real, y::Real, z::Real, ypr::NamedTuple{(:yaw, :pitch, :roll)})::Pose
  Pose(x, y, z, RotZYX(ypr.yaw, ypr.pitch, ypr.roll))
end

function Base.isapprox(a::Pose, b::Pose; kwargs...)::Bool
  (isapprox(a.pos, b.pos; kwargs...) &&
   isapprox(a.orientation, b.orientation; kwargs...))
end


### Constants ###

"""
The identity quaternion, representing the identity orientation.
"""
IDENTITY_ORN = UnitQuaternion(1, 0, 0, 0)

"""
Identity pose, a.k.a. the relative pose of any coordinate frame relative to itself.

This is the identity element for the pose group (the identity element is the same in both
group structures).
"""
const IDENTITY_POSE = Pose([0, 0, 0], IDENTITY_ORN)


### Pretty-printing ###

function Base.show(io::IO, pose::Pose)
  # Not sure whether most people will want (yaw, pitch, roll) or quaternion
  # components here.
  print(io, strWithQuat(pose))
end

function strWithYPR(pose::Pose)::String
  (yaw, pitch, roll) = _yawPitchRoll(pose.orientation)
  "Pose⟨pos=$(pose.pos), orientation=(yaw=$yaw, pitch=$pitch, roll=$roll)⟩"
end

function strWithQuat(pose::Pose)::String
  q = UnitQuaternion(pose.orientation)
  "Pose⟨pos=$(pose.pos), orientation=(w=$(q.w), x=$(q.x), y=$(q.y), z=$(q.z))⟩"
end


function _yawPitchRoll(orn::Rotation{3})
  ypr = RotZYX(orn)
  (yaw=ypr.theta1, pitch=ypr.theta2, roll=ypr.theta3)
end

componentsWXYZ(q::UnitQuaternion) = @SVector([q.w, q.x, q.y, q.z])

"""
Like `isapprox`, but does not consider a quaternion to be equivalent to its
negative (even though they correspond to the same rotation matrix).  Note that
this is stricter than `Base.isapprox`, since for a `Rotations.UnitQuaternion`
`q`, we have `-q ≈ q` and in fact `-q == q`.
"""
function isapproxIncludingQuaternionSign(a::Pose, b::Pose; kwargs)::Bool
  (isapprox(a.pos, b.pos; kwargs...) &&
   isapprox(componentsWXYZ(UnitQuaternion(a.orientation)), 
            componentsWXYZ(UnitQuaternion(b.orientation));
            kwargs...))
end



# Operations for combining poses.
Base.:(*)(a::Pose, b::Pose)::Pose = Pose(
    a.pos + a.orientation * b.pos,
    a.orientation * b.orientation)

Base.:(/)(a::Pose, b::Pose)::Pose = Pose(
    a.pos - (a.orientation / b.orientation) * b.pos,
    a.orientation / b.orientation)

Base.:(\)(a::Pose, b::Pose)::Pose = Pose(
    a.orientation \ (-a.pos + b.pos),
    a.orientation \ b.orientation)

Base.:(^)(pose::Pose, t::Real) = Pose(t * pose.pos,
                                      # quaternion exponentiation = SLERP
                                      quatPow(UnitQuaternion(pose.orientation),
                                              t))

function Base.inv(a::Pose)::Pose
  Pose(a.orientation \ -a.pos,
       inv(a.orientation))
end

# Action of a pose on a point.
function Base.:(*)(a::Pose, bpos::StaticVector{3, <:Real})
  a.pos + a.orientation * bpos
end

function Base.:(\)(a::Pose, bpos::StaticVector{3, <:Real})
  a.orientation \ (-a.pos + bpos)
end

# Convenience wrappers when the user supplies a Vector instead of a StaticVector
Base.:(*)(a::Pose, bpos::AbstractVector{<:Real}) = a * SVector{3}(bpos)
Base.:(\)(a::Pose, bpos::AbstractVector{<:Real}) = a \ SVector{3}(bpos)


"""
Vectorized pose–point multiplication.  Returns the matrix whose `i`th column is
`a * bpoints[:, i]`.

The matrix `bpoints` must have 3 rows, as each column represents a point in 3D
space.
"""
function Base.:(*)(a::Pose, bpoints::AbstractMatrix{<:Real})
  size(bpoints, 1) == 3 || error(
      "Must pass a 3×N matrix (one column per point)")
  a.pos .+ a.orientation * bpoints
end

"""
Vectorized version of pose–point left division.  Returns the matrix whose `i`th
column is `a \\ bpoints[:, i]`.

The matrix `bpoints` must have 3 rows, as each column represents a point in 3D
space.
"""
function Base.:(\)(a::Pose, bpoints::AbstractMatrix{<:Real})
  size(bpoints, 1) == 3 || error(
      "Must pass a 3×N matrix (one column per point)")
  a.orientation \ (-a.pos .+ bpoints)
end


"""
!!! note "TODO"
    This code mostly duplicates `GenDirectionalStats.hopf`.  The two should
    probably be consolidated into one.

Returns a rotation that carries the z-axis to `newZ`, with the remaining degree
of freedom determined by `planarAngle` as described below.

Start with the case `planarAngle = 0`.  In that case, the returned rotation is
the unique (except at singularities) rotation that carries `[0, 0, 1]` to
`newZ` "along a great circle" (more precisely: the unique rotation that carries
`[0, 0, 1]` to `newZ` and whose equator contains `[0, 0, 1]` and `newZ`;
"equator" means the unique great circle that is fixed setwise by the rotation).

Next consider the general case.  This works the same as the above special case,
except that we precede that rotation with a rotation by angle `planarAngle`
around `[0, 0, 1]` (or equivalently, we follow that rotation with a rotation by
`planarAngle` around `newZ`).

The name of this function comes from the fact that we are using geodesics
(great circles) to define a coordinate chart on the fiber over `newZ`  in the
Hopf fibration.

See also: [`invGeodesicHopf`](@ref)
"""
function geodesicHopf(newZ::StaticVector{3, <:Real}, planarAngle::Real)
  @assert norm(newZ) ≈ 1
  zUnit = @SVector([0, 0, 1])
  if newZ ≈ -zUnit
    @warn "Singularity: anti-parallel z-axis, rotation has an undetermined degree of freedom"
    axis = @SVector([1, 0, 0])
    geodesicAngle = π
  elseif newZ ≈ zUnit
    # Choice of axis doesn't matter here as long as it's nonzero
    axis = @SVector([1, 0, 0])
    geodesicAngle = 0
  else
    axis = cross(zUnit, newZ)
    @assert !(axis ≈ zero(axis)) || newZ ≈ zUnit
    geodesicAngle = let θ = asin(clamp(norm(axis), -1, 1))
      dot(zUnit, newZ) > 0 ? θ : π - θ
    end
  end
  return (AngleAxis(geodesicAngle, axis...) *
          AngleAxis(planarAngle, zUnit...))
end

geodesicHopf(newZ::AbstractVector{<:Real}, planarAngle::Real) = geodesicHopf(
    SVector{3}(newZ), planarAngle)


"""
Inverse function of [`geodesicHopf`](@ref).

Satisfies the round-trip conditions

    geodesicHopf(invGeodesicHopf(r)...) == r

and

    invGeodesicHopf(geodesicHopf(newZ, planarAngle))
    == (newZ=newZ, planarAngle=planarAngle)
"""
function invGeodesicHopf(r::Rotation{3})::NamedTuple{(:newZ, :planarAngle)}
  zUnit = @SVector([0, 0, 1])
  newZ = r * zUnit
  if newZ ≈ -zUnit
    @warn "Singularity: anti-parallel z-axis, planarAngle is undetermined"
    planarAngle = 0
  else
    # Solve `planarRot == AngleAxis(planarAngle, zUnit...)` for `planarAngle`
    planarRot = AngleAxis(geodesicHopf(newZ, 0) \ r)
    axis = @SVector([planarRot.axis_x, planarRot.axis_y, planarRot.axis_z])
    # `axis` is either `zUnit` or `-zUnit`, and we need to ensure that it's
    # `zUnit`.  (Exception: the degenerate case `planarAngle == 0`)
    if axis[3] < 0
      axis = -axis
      planarAngle = -planarRot.theta
    else
      planarAngle = planarRot.theta
    end
    atol = 1e-14
    @assert isapprox(axis, zUnit; atol=atol) ||
            abs(rem2pi(planarAngle, RoundNearest)) < atol
  end
  return (newZ=newZ, planarAngle=planarAngle)
end


# Second group structure on poses: Direct product ``ℝ^3 × SO_3(ℝ)``.
⊗(a::Pose, b::Pose) = Pose(a.pos + b.pos, a.orientation * b.orientation)
⊘(a::Pose, b::Pose) = Pose(a.pos - b.pos, a.orientation / b.orientation)
⦸(a::Pose, b::Pose) = Pose(b.pos - a.pos, a.orientation \ b.orientation)


"""
Interpolates between the identity pose and `b`.

Namely, `interp(b, 0) == IDENTITY_POSE` and `interp(b, 1) == b`.  The position
is interpolated linearly and the orientation is interpolated by quaternion
SLERP.  That is, this interpolation treats position and orientation
independently, as in the [`⊗`](@ref) operation (not the [`*`](@ref
Base.:*(::Pose, ::Pose)) operation).
"""
interp(b::Pose, t::Real) = Pose(t * b.pos,
                                quatPow(UnitQuaternion(b.orientation), t))

"""
Like [`interp`](@ref interp(::Pose, ::Real)), but interpolates between two
given poses rather than always starting at the identity.  That is,

    interp(a, b, 0) == a
    interp(a, b, 1) == b

and as a special case, we have

    interp(b, t) == interp(IDENTITY_POSE, b, t)
"""
interp(a::Pose, b::Pose, t::Real) = a * interp(a \ b, t)


function quatPow(q::UnitQuaternion, t::Real)
  # TODO: Once https://github.com/JuliaGeometry/Rotations.jl/issues/126 is
  # fixed, this special case won't be necessary
  if t == 0 || q == one(UnitQuaternion) || q == -one(UnitQuaternion)
    return one(UnitQuaternion)
  end
  return exp(t * log(q))
end


end  # module PoseComposition
