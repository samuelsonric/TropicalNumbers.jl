struct MaxMin <: AbstractLatticeAlgebra end

"""
    TropicalMaxMin{T} <: AbstractSemiring

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
3.0

julia> TropicalMaxMin(1.0) * TropicalMaxMin(3.0)
1.0

julia> zero(TropicalMaxMinF64)
-Inf

julia> one(TropicalMaxMinF64)
Inf
```
"""
const TropicalMaxMin = Semiring{MaxMin}

add_alg(::Type{MaxMin}, a, b) = max(a, b)
mul_alg(::Type{MaxMin}, a, b) = min(a, b)

zero_alg(::Type{MaxMin}, ::Type{T}) where {T} = neginf(T)
one_alg(::Type{MaxMin}, ::Type{T}) where {T} = posinf(T)

ldiv_alg(::Type{MaxMin}, a, b::T) where {T} = ifelse(a <= b, one_alg(MaxMin, T), b)
rdiv_alg(::Type{MaxMin}, b, a) = ldiv_alg(MaxMin, a, b)

leq_alg(::Type{MaxMin}, a, b) = a <= b
lt_alg(::Type{MaxMin}, a, b) = a < b

add_fast_alg(::Type{MaxMin}, a, b) = Base.FastMath.max_fast(a, b)
mul_fast_alg(::Type{MaxMin}, a, b) = Base.FastMath.min_fast(a, b)

# --------------- #
# other operators #
# --------------- #

function Base.isless(a::TropicalMaxMin, b::TropicalMaxMin)
    return a < b
end
