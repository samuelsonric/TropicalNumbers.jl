struct MaxPlus <: AbstractTropicalAlgebra end

"""
    TropicalMaxPlus{T} = Tropical{T} <: AbstractSemiring

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
const TropicalMaxPlus = Semiring{MaxPlus}
const Tropical = TropicalMaxPlus

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
3.0ₜ

julia> TropicalMaxMin(1.0) * TropicalMaxMin(3.0)
1.0ₜ

julia> zero(TropicalMaxMinF64)
-Infₜ

julia> one(TropicalMaxMinF64)
Infₜ
```
"""
const TropicalMaxMin = Semiring{LatticeAlgebra{MaxPlus}}

add_alg(::Type{MaxPlus}, a, b) = max(a, b)
inf_alg(::Type{MaxPlus}, a, b) = min(a, b)

#
#   -Inf +  Inf = -Inf
#    Inf + -Inf = -Inf
#
mul_alg(::Type{MaxPlus}, a, b) = a + b

function mul_alg(::Type{MaxPlus}, a::T, b::T) where {T <: Rational}
    ⊤ = typemax(T)
    ⊥ = typemin(T)

    if a == ⊥ && b == ⊤ || a == ⊤ && b == ⊥
        c = ⊥
    else
        c = a + b
    end

    return c
end

exp_alg(::Type{MaxPlus}, a, b) = a * b
inv_alg(::Type{MaxPlus}, a) = -a

zero_alg(::Type{MaxPlus}, ::Type{T}) where {T} = neginf(T)
typemax_alg(::Type{MaxPlus}, ::Type{T}) where {T} = posinf(T)
one_alg(::Type{MaxPlus}, ::Type{T}) where {T} = zero(T)

#
#    Inf -  Inf =  Inf
#   -Inf - -Inf =  Inf
#
ldiv_alg(::Type{MaxPlus}, a, b) = b - a

function ldiv_alg(::Type{MaxPlus}, a::T, b::T) where {T <: Rational}
    ⊤ = typemax(T)
    ⊥ = typemin(T)

    if a == b == ⊥ || a == b == ⊤
        c = ⊤
    else
        c = b - a
    end

    return c
end

leq_alg(::Type{MaxPlus}, a, b) = a <= b
lt_alg(::Type{MaxPlus}, a, b) = a < b

add_fast_alg(::Type{MaxPlus}, a, b) = Base.FastMath.max_fast(a, b)
mul_fast_alg(::Type{MaxPlus}, a, b) = Base.FastMath.add_fast(a, b)
inf_fast_alg(::Type{MaxPlus}, a, b) = Base.FastMath.min_fast(a, b)
ldiv_fast_alg(::Type{MaxPlus}, a, b) = Base.FastMath.sub_fast(b, a)

# --------------- #
# other operators #
# --------------- #

function Base.isless(a::TropicalMaxPlus, b::TropicalMaxPlus)
    return a < b
end

function Base.isless(a::TropicalMaxMin, b::TropicalMaxMin)
    return a < b
end
