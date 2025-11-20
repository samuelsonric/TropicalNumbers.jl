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
3.0

julia> TropicalMaxPlus(1.0) * TropicalMaxPlus(3.0)
4.0

julia> one(TropicalMaxPlusF64)
0.0

julia> zero(TropicalMaxPlusF64)
-Inf
```
"""
const Tropical = Semiring{MaxPlus}
const TropicalMaxPlus = Tropical

add_alg(::Type{MaxPlus}, a, b) = max(a, b)
inf_alg(::Type{MaxPlus}, a, b) = min(a, b)

#
#   -Inf +  Inf = -Inf
#    Inf + -Inf = -Inf
#
mul_alg(::Type{MaxPlus}, a, b) = a + b
exp_alg(::Type{MaxPlus}, a, b) = a * b
inv_alg(::Type{MaxPlus}, a) = -a

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

rdiv_alg(::Type{MaxPlus}, b, a) = ldiv_alg(MaxPlus, a, b)

leq_alg(::Type{MaxPlus}, a, b) = a <= b
lt_alg(::Type{MaxPlus}, a, b) = a < b

add_fast_alg(::Type{MaxPlus}, a, b) = Base.FastMath.max_fast(a, b)
mul_fast_alg(::Type{MaxPlus}, a, b) = Base.FastMath.add_fast(a, b)
ldiv_fast_alg(::Type{MaxPlus}, a, b) = Base.FastMath.sub_fast(b, a)
rdiv_fast_alg(::Type{MaxPlus}, b, a) = ldiv_fast_alg(MaxPlus, a, b)

# --------------- #
# other operators #
# --------------- #

function Base.isless(a::Tropical, b::Tropical)
    return a < b
end
