"""
    TropicalMaxPlus{T} = Tropical{T} <: AbstractIdempotentSemiring{T}

TropicalMaxPlus is a semiring algebra, can be described by
* Tropical (TropicalMaxPlus), (ℝ, max, +, -Inf, 0).

It maps
* `+` to `max` in regular algebra,
* `*` to `+` in regular algebra,
* `1` to `0` in regular algebra,
* `0` to `-Inf` in regular algebra (for integer content types, this is chosen as a small integer).

Example
-------------------------
```jldoctest; setup=:(using TropicalNumbers)
julia> TropicalMaxPlus(1.0) + TropicalMaxPlus(3.0)
3.0ₜ

julia> TropicalMaxPlus(1.0) * TropicalMaxPlus(3.0)
4.0ₜ

julia> one(TropicalMaxPlusF64)
0.0ₜ

julia> zero(TropicalMaxPlusF64)
-Infₜ
```
"""
struct Tropical{T} <: AbstractIdempotentSemiring{T}
    n::T
    Tropical{T}(x) where T = new{T}(T(x))
    function Tropical(x::T) where T
        new{T}(x)
    end
    function Tropical{T}(x::Tropical{T}) where T
        x
    end
    function Tropical{T1}(x::Tropical{T2}) where {T1,T2}
        new{T1}(T2(x.n))
    end
end

content(a::Tropical) = a.n

Base.:^(a::Tropical, b::Real) = Tropical(a.n * b)
Base.:^(a::Tropical, b::Integer) = Tropical(a.n * b)

#
#   -Inf +  Inf = -Inf
#    Inf + -Inf = -Inf
#
Base.:*(a::Tropical, b::Tropical) = Tropical(a.n + b.n)

function Base.:*(a::T, b::T) where {T <: Tropical{<:Rational}}
    if a == typemin(T) && b == typemax(T) || a == typemax(T) && b == typemin(T)
        c = typemin(T)
    else
        c = Tropical(a.n + b.n)
    end

    return c
end

inf(a::Tropical, b::Tropical) = Tropical(min(a.n, b.n))
sup(a::Tropical, b::Tropical) = Tropical(max(a.n, b.n))
Base.typemin(::Type{Tropical{T}}) where {T} = Tropical(neginf(T))
Base.typemax(::Type{Tropical{T}}) where {T} = Tropical(posinf(T))
Base.one(::Type{Tropical{T}}) where T = Tropical(zero(T))

# inverse and division
#Base.inv(x::Tropical) = Tropical(-x.n)

#
#    Inf -  Inf =  Inf
#   -Inf - -Inf =  Inf
#
Base.:\(y::Tropical, x::Tropical) = x / y
Base.:/(x::Tropical, y::Tropical) = Tropical(x.n - y.n)

function Base.:/(b::T, a::T) where {T <: Tropical{<:Rational}}
    if a == b == typemin(T) || a == b == typemax(T) 
        c = typemax(T)
    else
        c = Tropical(b.n - a.n)
    end

    return c
end

Base.div(x::Tropical, y::Tropical) = x / y

# promotion rules
Base.promote_rule(::Type{Tropical{T1}}, b::Type{Tropical{T2}}) where {T1, T2} = Tropical{promote_type(T1,T2)}
