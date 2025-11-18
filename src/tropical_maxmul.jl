

"""
    TropicalMaxMul{T} <: AbstractIdempotentSemiring{T}

TropicalMaxMul is a semiring algebra, can be described by
* TropicalMaxMul, (ℝ⁺, max, ⋅, 0, 1).

It maps
* `+` to `max` in regular algebra,
* `*` to `*` in regular algebra,
* `1` to `1` in regular algebra,
* `0` to `0` in regular algebra.

Example
-------------------------
```jldoctest; setup=:(using TropicalNumbers)
julia> TropicalMaxMul(1.0) + TropicalMaxMul(3.0)
3.0ₓ

julia> TropicalMaxMul(1.0) * TropicalMaxMul(3.0)
3.0ₓ

julia> one(TropicalMaxMulF64)
1.0ₓ

julia> zero(TropicalMaxMulF64)
0.0ₓ
```
"""
struct TropicalMaxMul{T} <: AbstractIdempotentSemiring{T}
    n::T
    function TropicalMaxMul{T}(x) where T 
        new{T}(T(x))
    end
    function TropicalMaxMul(x::T) where T
        new{T}(x)
    end
    function TropicalMaxMul{T}(x::TropicalMaxMul{T}) where T
        x
    end
    function TropicalMaxMul{T1}(x::TropicalMaxMul{T2}) where {T1,T2}
        new{T1}(T2(x.n))
    end
end

content(a::TropicalMaxMul) = a.n

Base.:^(a::TropicalMaxMul, b::Real) = TropicalMaxMul(a.n ^ b)
Base.:^(a::TropicalMaxMul, b::Integer) = TropicalMaxMul(a.n ^ b)

#
#    0   * Inf = 0
#    Inf * 0   = 0
#
Base.:*(a::TropicalMaxMul, b::TropicalMaxMul) = TropicalMaxMul(a.n * b.n)

function Base.:*(a::T, b::T) where {T <: TropicalMaxMul{<:Rational}}
    if a == typemin(T) && b == typemax(T) || a == typemax(T) && b == typemin(T)
        c = typemin(T)
    else
        c = Tropical(a.n + b.n)
    end

    return c
end

inf(a::TropicalMaxMul, b::TropicalMaxMul) = TropicalMaxMul(min(a.n, b.n))
sup(a::TropicalMaxMul, b::TropicalMaxMul) = TropicalMaxMul(max(a.n, b.n))
Base.typemin(::Type{TropicalMaxMul{T}}) where T = TropicalMaxMul(neginf(T))
Base.typemax(::Type{TropicalMaxMul{T}}) where T = TropicalMaxMul(posinf(T))
Base.one(::Type{TropicalMaxMul{T}}) where T = TropicalMaxMul(one(T))

# inverse and division
#Base.inv(x::TropicalMaxMul) = TropicalMaxMul(inv(x.n))

#
#   Inf / Inf = Inf
#   0   / 0   = Inf
#
Base.:\(y::TropicalMaxMul, x::TropicalMaxMul) = x / y
Base.:/(x::TropicalMaxMul, y::TropicalMaxMul) = TropicalMaxMul(x.n / y.n)

function Base.:/(b::T, a::T) where {T <: TropicalMaxMul{<:Rational}}
    if a == b == typemin(T) || a == b == typemax(T)
        c = typemax(T)
    else
        c = Tropical(b.n - a.n)
    end

    return c
end

Base.div(x::TropicalMaxMul, y::TropicalMaxMul) = TropicalMaxMul(div(x.n, y.n))

# promotion rules
Base.promote_rule(::Type{TropicalMaxMul{T1}}, b::Type{TropicalMaxMul{T2}}) where {T1, T2} = TropicalMaxMul{promote_type(T1,T2)}
