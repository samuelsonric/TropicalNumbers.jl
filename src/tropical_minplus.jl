struct MinPlus <: AbstractTropicalAlgebra end

"""
    TropicalMinPlus{T} <: AbstractSemiring

TropicalMinPlus is a semiring algebra, can be described by
* TropicalMinPlus, (ℝ, min, +, Inf, 0).

It maps
* `+` to `min` in regular algebra,
* `*` to `+` in regular algebra,
* `1` to `0` in regular algebra,
* `0` to `Inf` in regular algebra (for integer content types, this is chosen as a large integer).

Example
-------------------------
```jldoctest; setup=:(using TropicalNumbers)
julia> TropicalMinPlus(1.0) + TropicalMinPlus(3.0)
1.0

julia> TropicalMinPlus(1.0) * TropicalMinPlus(3.0)
4.0

julia> one(TropicalMinPlusF64)
0.0

julia> zero(TropicalMinPlusF64)
Inf
```
"""
const TropicalMinPlus = Semiring{MinPlus}

add_alg(::Type{MinPlus}, a, b) = min(a, b)
inf_alg(::Type{MinPlus}, a, b) = max(a, b)

#
#   -Inf +  Inf =  Inf
#    Inf + -Inf =  Inf
#
mul_alg(::Type{MinPlus}, a, b) = a + b
exp_alg(::Type{MinPlus}, a, b) = a * b
inv_alg(::Type{MinPlus}, a) = -a

function mul_alg(::Type{MinPlus}, a::T, b::T) where {T <: Rational}
    ⊤ = typemax(T)
    ⊥ = typemin(T)

    if a == ⊥ && b == ⊤ || a == ⊤ && b == ⊥
        c = ⊤
    else
        c = a + b
    end

    return c
end

zero_alg(::Type{MinPlus}, ::Type{T}) where {T} = posinf(T)
typemax_alg(::Type{MinPlus}, ::Type{T}) where {T} = neginf(T)
one_alg(::Type{MinPlus}, ::Type{T}) where {T} = zero(T)

#
#    Inf -  Inf = -Inf
#   -Inf - -Inf = -Inf
#
ldiv_alg(::Type{MinPlus}, a, b) = b - a

function ldiv_alg(::Type{MinPlus}, a::T, b::T) where {T <: Rational}
    ⊤ = typemax(T)
    ⊥ = typemin(T)

    if a == b == ⊥ || a == b == ⊤
        c = ⊥
    else
        c = b - a
    end

    return c
end

leq_alg(::Type{MinPlus}, a, b) = a >= b
lt_alg(::Type{MinPlus}, a, b) = a > b

add_fast_alg(::Type{MinPlus}, a, b) = Base.FastMath.min_fast(a, b)
mul_fast_alg(::Type{MinPlus}, a, b) = Base.FastMath.add_fast(a, b)
inf_fast_alg(::Type{MinPlus}, a, b) = Base.FastMath.max_fast(a, b)
ldiv_fast_alg(::Type{MinPlus}, a, b) = Base.FastMath.sub_fast(b, a)

# --------------- #
# other operators #
# --------------- #

function Base.isless(a::TropicalMinPlus, b::TropicalMinPlus)
    return a < b
end
