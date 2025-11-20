"""
    AbstractSemiring{T} <: Number

A [`semiring`](https://en.wikipedia.org/wiki/Semiring) is a quintuple (R, +, ×, 0, 1), where

* (R, +, 0) is a commutative monoid
* (R, ×, 1) is a monoid
* multiplication (×) distributes over addition (+)
* zero (0) is absorbing

This package implements the following semirings.

* TropicalAndOr, ({0, 1}, or, and, 0, 1);
* TropicalMaxPlus, ([-∞, +∞], max, +, -∞, 0);
* TropicalMinPlus, ([-∞, +∞], min, +, +∞, 0);
* TropicalMaxMin, ([-∞, +∞], max, min, -∞, +∞)
* TropicalMaxMul, ([0, +∞], max, ×, 0, 1).

Fast semiring matrix multiplication is implemented in the following libraries.

* [`TropicalGEMM`](https://github.com/TensorBFS/TropicalGEMM.jl/)
* [`CuTropicalGEMM`](https://github.com/ArrogantGao/CuTropicalGEMM.jl/)

"""
abstract type AbstractSemiring <: Number end

function Base.zero(::T) where {T <: AbstractSemiring}
    return zero(T)
end

function Base.one(::T) where {T <: AbstractSemiring}
    return one(T)
end

function Base.:*(a::T, b::Bool) where {T <: AbstractSemiring}
    return b ? a : zero(T)
end

function Base.:*(b::Bool, a::T) where {T <: AbstractSemiring}
    return b ? a : zero(T)
end

struct Semiring{A <: AbstractSemiringAlgebra, T} <: AbstractSemiring
    n::T

    function Semiring{A, T}(n) where {A <: AbstractSemiringAlgebra, T}
        return new{A, T}(n)
    end
end

function Semiring{A}(n::T) where {A <: AbstractSemiringAlgebra, T}
    return Semiring{A, T}(n)
end

function Semiring{A}(a::Semiring{A}) where {A <: AbstractSemiringAlgebra}
    return Semiring{A}(content(a))
end

function Semiring{A, T}(a::Semiring{A}) where {A <: AbstractSemiringAlgebra, T}
    return Semiring{A, T}(content(a))
end

function content(a::Semiring)
    return a.n
end

function content(::Type{Semiring{A, T}}) where {A <: AbstractSemiringAlgebra, T}
    return T
end

function Base.show(io::IO, a::Semiring)
    print(io, content(a))
    return
end

function Base.isnan(a::Semiring)
    return isnan(content(a))
end

function Base.isinf(a::Semiring)
    return isinf(content(a))
end

function Base.rand(rng::AbstractRNG, sampler::SamplerType{Semiring{A, T}}) where {A <: AbstractSemiringAlgebra, T}
    n =  rand(rng, T)
    return Semiring{A}(n)
end

function Base.:(==)(a::Semiring{A}, b::Semiring{A}) where {A <: AbstractSemiringAlgebra}
    return content(a) == content(b)
end

function Base.isapprox(a::Semiring{A}, b::Semiring{A}; kw...) where {A <: AbstractSemiringAlgebra}
    return isapprox(content(a), content(b); kw...)
end

function Base.isapprox(a::AbstractArray{<:Semiring{A}}, b::AbstractArray{<:Semiring{A}}; kw...) where {A <: AbstractSemiringAlgebra}
    return all(isapprox.(a, b; kw...))
end

function Base.promote_rule(::Type{Semiring{A, T}}, ::Type{Semiring{A, U}}) where {A <: AbstractSemiringAlgebra, T, U}
    V = promote_type(T, U)
    return Semiring{A, V}
end

# --------- #
# Semirings #
# --------- #

function Base.zero(::Type{Semiring{A, T}}) where {A <: AbstractSemiringAlgebra, T}
    n = zero_alg(A, T)
    return Semiring{A}(n)
end

function Base.one(::Type{Semiring{A, T}}) where {A <: AbstractSemiringAlgebra, T}
    n = one_alg(A, T)
    return Semiring{A}(n)
end

function Base.:+(a::Semiring{A}, b::Semiring{A}) where {A <: AbstractSemiringAlgebra}
    n = add_alg(A, content(a), content(b))
    return Semiring{A}(n)
end

function Base.FastMath.add_fast(a::Semiring{A}, b::Semiring{A}) where {A <: AbstractSemiringAlgebra}
    n = add_fast_alg(A, content(a), content(b))
    return Semiring{A}(n)
end

function Base.:*(a::Semiring{A}, b::Semiring{A}) where {A <: AbstractSemiringAlgebra}
    n = mul_alg(A, content(a), content(b))
    return Semiring{A}(n)
end

function Base.FastMath.mul_fast(a::Semiring{A}, b::Semiring{A}) where {A <: AbstractSemiringAlgebra}
    n = mul_fast_alg(A, content(a), content(b))
    return Semiring{A}(n)
end

function Base.fma(a::Semiring{A}, b::Semiring{A}, c::Semiring{A}) where {A <: AbstractSemiringAlgebra}
    n = mul_add_alg(A, content(a), content(b), content(c))
    return Semiring{A}(n)
end

# --------- #
# Quantales #
# --------- #

function Base.typemin(::Type{T}) where {T <: Semiring{<:AbstractQuantaleAlgebra}}
    return zero(T)
end

function Base.typemin(::T) where {T <: Semiring{<:AbstractQuantaleAlgebra}}
    return zero(T)
end

function Base.typemax(::Type{Semiring{A, T}}) where {A <: AbstractQuantaleAlgebra, T}
    n = typemax_alg(A, T)
    return Semiring{A}(n)
end

function Base.typemax(::T) where {T <: Semiring{<:AbstractQuantaleAlgebra}}
    return typemax(T)
end

function inf(a::Semiring{A}, b::Semiring{A}) where {A <: AbstractQuantaleAlgebra}
    n = inf_alg(A, content(a), content(b))
    return Semiring{A}(n)
end

function inf_fast(a, b)
    return inf(a, b)
end

function inf_fast(a::Semiring{A}, b::Semiring{A}) where {A <: AbstractQuantaleAlgebra}
    n = inf_fast_alg(A, content(a), content(b))
    return Semiring{A}(n)
end

function Base.:\(a::Semiring{A}, b::Semiring{A}) where {A <: AbstractQuantaleAlgebra}
    n = ldiv_alg(A, content(a), content(b))
    return Semiring{A}(n)
end

function ldiv_fast(a, b)
    return a \ b
end

function ldiv_fast(a::Semiring{A}, b::Semiring{A}) where {A <: AbstractQuantaleAlgebra}
    n = ldiv_fast_alg(A, content(a), content(b))
    return Semiring{A}(n)
end

function fli(a, b, c)
    return (a \ b) ∧ c
end

function fli(a::Semiring{A}, b::Semiring{A}, c::Semiring{A}) where {A <: AbstractQuantaleAlgebra}
    n = inf_ldiv_alg(A, content(a), content(b), content(c))
    return Semiring{A}(n)
end

function Base.:/(b::Semiring{A}, a::Semiring{A}) where {A <: AbstractQuantaleAlgebra}
    n = rdiv_alg(A, content(b), content(a))
    return Semiring{A}(n)
end

function Base.:/(b::T, a::Bool) where {T <: Semiring{<:AbstractQuantaleAlgebra}}
     return a ? b : typemax(T)
end

function Base.:/(b::Bool, a::T) where {T <: Semiring{<:AbstractQuantaleAlgebra}}
    return b ? one(T) / a : zero(T)
end

function Base.FastMath.div_fast(b::Semiring{A}, a::Semiring{A}) where {A <: AbstractQuantaleAlgebra}
    n = rdiv_fast_alg(A, content(b), content(a))
    return Semiring{A}(n)
end

function fri(b, a, c)
    return (b / a) ∧ c
end

function fri(b::Semiring{A}, a::Semiring{A}, c::Semiring{A}) where {A <: AbstractQuantaleAlgebra}
    n = inf_rdiv_alg(A, content(b), content(a), content(c))
    return Semiring{A}(n)
end

function Base.:<=(a::Semiring{A}, b::Semiring{A}) where {A <: AbstractQuantaleAlgebra}
    return leq_alg(A, content(a), content(b))
end

function Base.:<(a::Semiring{A}, b::Semiring{A}) where {A <: AbstractQuantaleAlgebra}
    return lt_alg(A, content(a), content(b))
end

# -------- #
# Tropical #
# -------- #

function Base.:^(a::Semiring{A}, b::Number) where {A <: AbstractTropicalAlgebra}
    n = exp_alg(A, content(a), b)
    return Semiring{A}(n)
end

function Base.:^(a::Semiring{A}, b::Integer) where {A <: AbstractTropicalAlgebra}
    n = exp_alg(A, content(a), b)
    return Semiring{A}(n)
end

function Base.inv(a::Semiring{A}) where {A <: AbstractTropicalAlgebra}
    n = inv_alg(A, content(a))
    return Semiring{A}(n)
end
