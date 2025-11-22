using Base.FastMath: mul_fast, add_fast, div_fast
using TropicalNumbers: ∧

function test_semiring(a::T, b::T, c::T) where {T <: AbstractSemiring}
    # 1 is the multiplicative identity
    @test one(T) * a ≈ a * one(T) ≈ a

    # 0 is the additive identity
    @test zero(T) + a ≈ a + zero(T) ≈ a

    # 0 is absorbing
    @test zero(T) * a ≈ a * zero(T) ≈ zero(T)

    # multiplication associates
    @test (a * b) * c ≈ a * (b * c)

    # addition associates
    @test (a + b) + c ≈ a + (b + c)    

    # addition commutes
    @test a + b ≈ b + a

    # multiplication right-distributes over addition
    @test (a + b) * c ≈ a * c + b * c

    # multiplication left-distributes over addition
    @test a * (b + c) ≈ a * b + a * c
end

function test_quantale(a::T, b::T, c::T) where {T <: AbstractSemiring}
    test_semiring(a, b, c)

    # ⊤ is the identity for infimum
    @test typemax(T) ∧ a ≈ a ∧ typemax(T) ≈ a

    # addition and infimum are connected by the absorbtion law
    @test a + (a ∧ b) ≈ a ∧ (a + b) ≈ a

    # lattice is residuated
    @test (a * b <= c) == (b <= a \ c) == (a <= c / b)

    # partial order is reflexive
    @test a <= a

    # strict order is irreflexive
    @test !(a < a)

    # partial order is transitive
    @test (a <= b && b <= c) <= (a <= c)

    # strict order is transitive
    @test (a < b && b < c) <= (a < c)

    # partial ordering agrees with lattice structure
    @test (a <= b) == (a ≈ a ∧ b)
    @test (a <= b) == (b ≈ a + b)

    # strict ordering agrees with partial ordering
    @test (a < b) == (a != b && a <= b)

    # lattice is a Heyting algebra
    @test imp(a, a) == typemax(T)
    @test a ∧ imp(a, b) == a ∧ b
    @test b ∧ imp(a, b) == b
    @test imp(a, b ∧ c) == imp(a, b) ∧ imp(a, c)

    # complement agrees with implication
    @test not(a) == imp(a, zero(T))
end

function test_fast(a::T, b::T, c::T) where {T <: AbstractSemiring}
    @test a + b ≈ add_fast(a, b)
    @test a * b ≈ mul_fast(a, b)
    @test a ∧ b ≈ inf_fast(a, b)
    @test a / b ≈ div_fast(a, b)
    @test a \ b ≈ ldiv_fast(a, b)
    @test imp(a, b) ≈ imp_fast(a, b)

    @test (a * b) + c ≈ fma(a, b, c)
    @test (a ∧ b) + c ≈ fia(a, b, c)
    @test (a \ b) ∧ c ≈ fli(a, b, c)
    @test (a / b) ∧ c ≈ fri(a, b, c)
    @test imp(a, b) ∧ c ≈ fii(a, b, c) 
end

function test_isless(a::T, b::T) where {T <: AbstractSemiring}
    @test isless(a, b) == (a < b)
end

@testset "printing" begin
    types = (
        TropicalMinPlus,
        TropicalMaxPlus,
        TropicalMaxMul,
        TropicalBitwise,
        TropicalMaxMin,
        TropicalMinPlusF64,
        TropicalMaxPlusF64,
        TropicalMaxMulF64,
        TropicalAndOr,
        TropicalBitwiseI64,
        TropicalMaxMinF64,
    )

    for T in types
        @test repr(T) isa String
    end
end

@testset "interface" begin
    types = (
        TropicalMinPlusF64,
        TropicalMaxPlusF64,
        TropicalMaxMulF64,
        TropicalAndOr,
        TropicalMaxMinF64,
    )

    for T in types
        a = rand(T)
        b = rand(T)

        test_isless(a, b)
        test_isless(a, a)
        test_isless(a, zero(T))
        test_isless(a, typemax(T))
    end

    for T in types
        a = rand(T)
        b = rand(T)
        c = rand(T)

        test_fast(a, b, c)
        test_quantale(a, b, c)
        test_quantale(zero(T), b, c)
    end

    types = (
        TropicalMinPlusF64,
        TropicalMaxPlusF64,
        TropicalMaxMulF64,
        TropicalAndOr,
        TropicalBitwiseI64,
        TropicalMaxMinF64,
    )

    for T in types
        a = rand(T)
        b = rand(T)
        c = rand(T)

        test_fast(a, b, c)
        test_quantale(a, b, c)
        test_quantale(zero(T), b, c)
    end

    types = (
        TropicalMinPlus{Rational{Int}},
        TropicalMaxPlus{Rational{Int}},
        TropicalMaxMul{Rational{Int}},
    )

    a = 1 // 2
    b = 2 // 3
    c = 3 // 4

    for T in types
        test_fast(T(a), T(b), T(c))
        test_quantale(T(a), T(b), T(c))
        test_quantale(zero(T), T(b), T(c))
        test_quantale(T(a), zero(T), zero(T))
        test_quantale(typemax(T), T(b), T(c))
        test_quantale(typemax(T), zero(T), T(c))
        test_quantale(typemax(T), T(b), typemax(T))
    end
end
