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
abstract type AbstractSemiring{T} <: Number end

function content(::Type{<:AbstractSemiring{T}}) where {T}
    return T
end

function Base.show(io::IO, a::AbstractSemiring)
    n = content(a)
    print(io, "$(n)s")
    return
end

function Base.show(io::IO, ::MIME"text/plain", a::AbstractSemiring)
    return show(io, a)
end

function Base.isnan(x::AbstractSemiring)
    return isnan(content(a))
end

function Base.isapprox(a::T, b::T; kw...) where {T <: AbstractSemiring}
    return isapprox(content(a), content(b); kw...)
end

function Base.isapprox(A::AbstractArray{<:AbstractSemiring}, B::AbstractArray{<:AbstractSemiring}; kw...)
    return all(isapprox.(A, B; kw...))
end

function Base.typemin(::T) where {T <: AbstractSemiring}
    return typemin(T)
end

function Base.typemax(::T) where {T <: AbstractSemiring}
    return typemax(T)
end

function Base.zero(::T) where {T <: AbstractSemiring}
    return zero(T)
end

function Base.one(::T) where {T <: AbstractSemiring}
    return one(T)
end

function Base.inv(a::T) where {T <: AbstractSemiring}
    return one(T) / a
end

function invint(a::T) where {T <: AbstractSemiring}
    return div(one(T), a)
end

function Base.:*(a::T, b::Bool) where {T <: AbstractSemiring}
    return b ? a : zero(T)
end

function Base.:*(b::Bool, a::T) where {T <: AbstractSemiring}
    return b ? a : zero(T)
end

function Base.:/(a::T, b::Bool) where {T <: AbstractSemiring}
    return b ? a : a / zero(T)
end

function Base.div(a::T, b::Bool) where {T <: AbstractSemiring}
    return b ? a : div(a, zero(T))
end

function Base.:/(b::Bool, a::T) where {T <: AbstractSemiring}
    return b ? inv(a) : zero(T)
end

function Base.div(b::Bool, a::T) where {T <: AbstractSemiring}
    return b ? invint(a) : zero(T)
end

function Base.:(==)(a::T, b::T) where {T <: AbstractSemiring}
    return content(a) == content(b)
end

function Base.isless(a::T, b::T) where {T <: AbstractSemiring}
    return a < b || !(b >= a) || isless(content(a), content(b))
end

"""
    AbstractIdempotentSemiring{T} <: AbstractSemiring{T}

A semiring (R, +, ×, 0, 1) is additively idempotent if

* x + x = x

for all x ∈ R.
"""
abstract type AbstractIdempotentSemiring{T} <: AbstractSemiring{T} end

function Base.zero(::Type{T}) where {T <: AbstractIdempotentSemiring}
    return typemin(T)
end

function Base.:+(a::T, b::T) where {T <: AbstractIdempotentSemiring}
    return sup(a, b)
end

function Base.:/(a::T, b::Bool) where {T <: AbstractIdempotentSemiring}
    return b ? a : typemax(a)
end

function Base.:div(a::T, b::Bool) where {T <: AbstractIdempotentSemiring}
    return b ? a : typemax(a)
end

function Base.:<=(a::T, b::T) where {T <: AbstractIdempotentSemiring}
    return sup(a, b) == b
end

function Base.:>=(a::T, b::T) where {T <: AbstractIdempotentSemiring}
    return b <= a
end

function Base.:<(a::T, b::T) where {T <: AbstractIdempotentSemiring}
    return a != b && a <= b
end

function Base.:>(a::T, b::T) where {T <: AbstractIdempotentSemiring}
    return b < a
end

"""
    AbstractSimpleSemiring{T} <: AbstractIdemptentSemiring{T}

A additively idempotent semiring (R, +, ×, 0, 1) is simple if

* a × b = inf(a, b)

for all a, b ∈ R.
"""
abstract type AbstractSimpleSemiring{T} <: AbstractIdempotentSemiring{T} end

function Base.:*(a::T, b::T) where {T <: AbstractSimpleSemiring}
    return inf(a, b)
end

function Base.one(::Type{T}) where {T <: AbstractSimpleSemiring}
    return typemax(T)
end

function Base.inv(a::T) where {T <: AbstractSimpleSemiring}
    return typemax(T)
end

function invint(a::T) where {T <: AbstractSimpleSemiring}
    return typemax(T)
end
