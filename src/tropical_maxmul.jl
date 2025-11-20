struct MaxMul <: AbstractTropicalAlgebra end

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
3.0

julia> TropicalMaxMul(1.0) * TropicalMaxMul(3.0)
3.0

julia> one(TropicalMaxMulF64)
1.0

julia> zero(TropicalMaxMulF64)
0.0
```
"""
const TropicalMaxMul = Semiring{MaxMul}

add_alg(::Type{MaxMul}, a, b) = max(a, b)
inf_alg(::Type{MaxMul}, a, b) = min(a, b)

#
#    0   * Inf = 0
#    Inf * 0   = 0
#
mul_alg(::Type{MaxMul}, a, b) = a * b
exp_alg(::Type{MaxMul}, a, b) = a ^ b
inv_alg(::Type{MaxMul}, a) = inv(a)

function mul_alg(::Type{MaxMul}, a::T, b::T) where {T <: Rational}
    ⊤ = typemax(T)
    ⊥ = zero(T)

    if a == ⊥ && b == ⊤ || a == ⊤ && b == ⊥
        c = ⊥
    else
        c = a * b
    end

    return c
end

zero_alg(::Type{MaxMul}, ::Type{T}) where {T} = zero(T)
typemax_alg(::Type{MaxMul}, ::Type{T}) where {T} = posinf(T)
one_alg(::Type{MaxMul}, ::Type{T}) where {T} = one(T)

#
#   Inf / Inf = Inf
#   0   / 0   = Inf
#
ldiv_alg(::Type{MaxMul}, a, b) = b / a

function ldiv_alg(::Type{MaxMul}, a::T, b::T) where {T <: Rational}
    ⊤ = typemax(T)
    ⊥ = zero(T)

    if a == b == ⊥ || a == b == ⊤
        c = ⊤
    else
        c = b / a
    end

    return c
end

rdiv_alg(::Type{MaxMul}, b, a) = ldiv_alg(MaxMul, a, b)

leq_alg(::Type{MaxMul}, a, b) = a <= b
lt_alg(::Type{MaxMul}, a, b) = a < b

add_fast_alg(::Type{MaxMul}, a, b) = Base.FastMath.max_fast(a, b)
mul_fast_alg(::Type{MaxMul}, a, b) = Base.FastMath.mul_fast(a, b)
ldiv_fast_alg(::Type{MaxMul}, a, b) = Base.FastMath.div_fast(b, a)
rdiv_fast_alg(::Type{MaxMul}, b, a) = ldiv_fast_alg(MaxMul, a, b)

# --------------- #
# other operators #
# --------------- #

function Base.isless(a::TropicalMaxMul, b::TropicalMaxMul)
    return a < b
end
