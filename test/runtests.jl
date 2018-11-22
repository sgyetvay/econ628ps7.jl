using econ627ps7
using Test
g(x)=(x-1.0)^2
gprime(x)=2*(x-1)

@testset "Basic Newtons Method Tests" begin
    @test newton(g,gprime,10.0) ≈ 1.0  atol=0.00001
    @test newton(g,10)≈ 1.0 atol=0.00001
end;

@testset "Convergence" begin
    @test (newton(x->x^2+1,x₀=10)==nothing)
    @test (newton(x->x^2+1,x->2*x,x₀=10)==nothing)
end;

@testset "Max Iter" begin
    @test (newton(g,x₀=10,maxiter=5)==nothing)
    @test (newton(g,gprime,x₀=10,maxiter=5)==nothing)
end;

@testset "Tolerance" begin
    x1=(newton(g,x₀=10,tol=1E-6))
    x2=(newton(g,x₀=10,tol=1E-5))
    @test(abs(x1-1)<abs(x2-1))
end;

@testset "BigFloat" begin
    @test (newton(g,gprime,x₀=BigFloat(10))≈1.0)
    @test (newton(g,x₀=BigFloat(10))≈1.0)
end;
