"""
    TropicalMaxMin{T} <: AbstractSimpleSemiring{T}

TropicalMaxMin is a semiring algebra, can be described by
* TropicalMaxMin, (ℝ, max, min, -Inf, Inf).

It maps
* `+` to `max` in regular algebra,
* `*` to `min` in regular algebra,
* `0` to `-Inf` in regular algebra (for integer content types, this is a small integer).
* `1` to `Inf` in regular algebra, (for integer content types, this is a large integer)

Example
-------------------------
```jldoctest; setup=:(using TropicalNumbers)
julia> TropicalMaxMin(1.0) + TropicalMaxMin(3.0)
3.0ₛ

julia> TropicalMaxMin(1.0) * TropicalMaxMin(3.0)
1.0ₛ

julia> zero(TropicalMaxMinF64)
-Infₛ

julia> one(TropicalMaxMinF64)
Infₛ
```
"""
struct TropicalMaxMin{T} <: AbstractSimpleSemiring{T}
    n::T
end

function TropicalMaxMin(a::TropicalMaxMin)
    return TropicalMaxMin(a.n)
end

function TropicalMaxMin{T}(a::TropicalMaxMin) where {T}
    return TropicalMaxMin{T}(a.n)
end

function Base.promote_rule(::Type{TropicalMaxMin{U}}, ::Type{TropicalMaxMin{V}}) where {U, V}
    W = promote_type(U, V)
    return TropicalMaxMin{W}
end

function content(a::TropicalMaxMin)
    return a.n
end

function inf(a::TropicalMaxMin, b::TropicalMaxMin)
    n = min(a.n, b.n)
    return TropicalMaxMin(n)
end

function sup(a::TropicalMaxMin, b::TropicalMaxMin)
    n = max(a.n, b.n)
    return TropicalMaxMin(n)
end

function Base.typemin(::Type{TropicalMaxMin{T}}) where {T}
    n = neginf(T)
    return TropicalMaxMin(n)
end

function Base.typemax(::Type{TropicalMaxMin{T}}) where {T}
    n = posinf(T)
    return TropicalMaxMin(n)
end

function Base.:\(a::TropicalMaxMin, b::TropicalMaxMin)
    return b / a
end

function Base.:/(b::TropicalMaxMin, a::TropicalMaxMin)
    return a <= b ? typemax(b) : b
end

function Base.:div(b::TropicalMaxMin, a::TropicalMaxMin)
    return b / a
end
