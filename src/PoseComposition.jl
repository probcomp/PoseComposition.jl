module PoseComposition

import Base: @kwdef
import LinearAlgebra: dot, norm, cross
import Rotations: AngleAxis, Rotation, QuatRotation, RotZYX, RotMatrix
import StaticArrays: StaticVector, SVector, @SVector

include("docstring_extensions.jl")

"""
Struct representing the pose (position and orientation) of an object.

Can represent either a relative pose (relative to some parent coordinate frame
which must be supplied by context) or an absolute pose.  An absolute pose is,
by definition, a pose relative to the world coordinate frame.

See [Operations on Poses](@ref) for documentation of this data structure and
the operations defined on it.
"""
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
IDENTITY_ORN = one(QuatRotation)

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
  q = QuatRotation(pose.orientation)
  "Pose⟨pos=$(pose.pos), orientation=(w=$(q.w), x=$(q.x), y=$(q.y), z=$(q.z))⟩"
end


function _yawPitchRoll(orn::Rotation{3})
  ypr = RotZYX(orn)
  (yaw=ypr.theta1, pitch=ypr.theta2, roll=ypr.theta3)
end

componentsWXYZ(q::QuatRotation) = @SVector([q.q.s, q.q.v1, q.q.v2, q.q.v3])

"""
Like `isapprox`, but does not consider a quaternion to be equivalent to its
negative (even though they correspond to the same rotation matrix).  Note that
this is stricter than `Base.isapprox`, since for a `Rotations.QuatRotation`
`q`, we have `-q ≈ q` and in fact `-q == q`.
"""
function isapproxIncludingQuaternionSign(a::Pose, b::Pose; kwargs)::Bool
  (isapprox(a.pos, b.pos; kwargs...) &&
   isapprox(componentsWXYZ(QuatRotation(a.orientation)), 
            componentsWXYZ(QuatRotation(b.orientation));
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
                                      quatPow(QuatRotation(pose.orientation),
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
                                quatPow(QuatRotation(b.orientation), t))

"""
Like [`interp`](@ref interp(::Pose, ::Real)), but interpolates between two
given poses rather than always starting at the identity.  That is,

    interp(a, b, 0) == a
    interp(a, b, 1) == b

and as a special case, we have

    interp(b, t) == interp(IDENTITY_POSE, b, t)
"""
interp(a::Pose, b::Pose, t::Real) = a * interp(a ⦸ b, t)


function quatPow(q::QuatRotation, t::Real)
  # TODO: Once https://github.com/JuliaGeometry/Rotations.jl/issues/126 is
  # fixed, this special case won't be necessary
  if t == 0 || q == one(QuatRotation) || q == -one(QuatRotation)
    return one(QuatRotation)
  end
  return RotMatrix{3}(exp(t * log(q)))
end


end  # module PoseComposition
