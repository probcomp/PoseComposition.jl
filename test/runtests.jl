import PoseComposition: Pose, IDENTITY_POSE, IDENTITY_ORN, interp, ⊗, ⊘
import Rotations: QuatRotation, RotZYX
import StaticArrays: SVector, @SVector
import Test: @test, @testset

@testset "Pose group operations" begin
  p1 = Pose(1.0, 1.1, 1.2, QuatRotation(√(0.1), √(0.2), √(0.3), √(0.4)))
  p2 = Pose(2.0, 2.1, 2.2, QuatRotation(√(0.2), √(0.2), √(0.3), √(0.3)))
  p3 = Pose(3.0, 3.1, 3.2, QuatRotation(√(0.3), √(0.3), √(0.3), √(0.1)))

  function testPosesApprox(a, b; kwargs...)
    @test isapprox(a.pos, b.pos; kwargs...)
    @test isapprox(a.orientation, b.orientation; kwargs...)
  end

  testPosesApprox((p1 * p2) * p3, p1 * (p2 * p3))
  testPosesApprox((p1 * p2) / p2, p1)
  testPosesApprox((p1 / p2) * p2, p1)
  testPosesApprox(p1 * (p1 \ p2), p2)
  testPosesApprox(inv(p1) * p1, IDENTITY_POSE; atol=1e-14)
  testPosesApprox(p1 * inv(p1), IDENTITY_POSE; atol=1e-14)
end


@testset "Commutative diagram: a * b.pos == (a * b).pos" begin
  a = Pose([1, 2, 3], RotZYX(0.4, 0.5, 0.6))
  bpos = @SVector([7, 8, 9])
  b = Pose(bpos, RotZYX(1.0, 1.1, 1.2))
  @test a * bpos ≈ (a * b).pos
end

@testset "Associativity: (a * b) * cpos == a * (b * cpos)" begin
  a = Pose([1, 2, 3], RotZYX(0.4, 0.5, 0.6))
  b = Pose([7, 8, 9], RotZYX(1.0, 1.1, 1.2))
  cpos = @SVector([13, 14, 15])
  @test (a * b) * cpos ≈ a * (b * cpos)
end

@testset "Left division: a::Pose \\ bpos::StaticVector == inv(a) * bpos" begin
  a = Pose([1, 2, 3], RotZYX(0.4, 0.5, 0.6))
  bpos = @SVector([7, 8, 9])
  @test a \ bpos ≈ inv(a) * bpos
end

@testset "Generic vector (non-`StaticVector`) wrappers, sanity check" begin
  a = Pose([1, 2, 3], RotZYX(0.4, 0.5, 0.6))
  bpos = @SVector([7, 8, 9])
  bpos_ = [7, 8, 9]
  @test a * bpos_ ≈ a * bpos
  @test a \ bpos_ ≈ a \ bpos
end

@testset "Equivalence of `Pose * Matrix` with `Pose * StaticVector` on each column" begin
  a = Pose([1, 2, 3], RotZYX(0.4, 0.5, 0.6))
  bpoints = reshape(7:21, (3, 5))
  @test a * bpoints ≈ reduce(hcat,
      a * SVector{3}(bpoint) for bpoint in eachcol(bpoints))
end

@testset "Edge case: Pose * (3-by-0 Matrix)" begin
  a = Pose([1, 2, 3], RotZYX(0.4, 0.5, 0.6))
  bpoints = zeros(3, 0)
  @test a * bpoints == zeros(3, 0)
end

@testset "Left division: a::Pose \\ bpoints::Matrix == inv(a) * bpoints" begin
  a = Pose([1, 2, 3], RotZYX(0.4, 0.5, 0.6))
  bpoints = reshape(7:21, (3, 5))
  @test a \ bpoints ≈ inv(a) * bpoints
end


@testset "`interp` endpoint conditions" begin
  a = Pose([10, 11, 12], IDENTITY_ORN)
  b = Pose([1, 2, 3], RotZYX(0.4, 0.5, 0.6))
  @test interp(a, b, 0) ≈ a
  @test interp(a, b, 1) ≈ b
end

@testset "`t ↦ interp(b, t)` is a homomorphism from (ℝ, +) to (poses, ⊗)" begin
  b = Pose([1, 2, 3], RotZYX(0.4, 0.5, 0.6))
  @test interp(b, 0.5) ⊗ interp(b, 0.7) ≈ interp(b, 1.2)
end

@testset "`interp(b, t) == b ⊗ ... ⊗ b` (`t` times) when `t` is an integer" begin
  b = Pose([1, 2, 3], RotZYX(0.4, 0.5, 0.6))
  @test interp(b, -1) ≈ IDENTITY_POSE ⊘ b
  @test interp(b, 2) ≈ b ⊗ b
  @test interp(b, 0.5) ⊗ interp(b, 0.5) ≈ b
end
