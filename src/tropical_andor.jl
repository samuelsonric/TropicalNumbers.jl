struct AndOr <: AbstractLatticeAlgebra end

"""
    TropicalAndOr <: AbstractSemiring

TropicalAndOr is a semiring algebra, can be described by
* TropicalAndOr, ([T, F], or, and, false, true).

It maps
* `+` to `or` in regular algebra,
* `*` to `and` in regular algebra,
* `1` to `true` in regular algebra,
* `0` to `false` in regular algebra.

For the parallel bit-wise version, see [`TropicalBitwise`](@ref).

Example
-------------------------
```jldoctest; setup=:(using TropicalNumbers)
julia> TropicalAndOr(true) + TropicalAndOr(false)
true

julia> TropicalAndOr(true) * TropicalAndOr(false)
false

julia> one(TropicalAndOr)
true

julia> zero(TropicalAndOr)
false
```
"""
const TropicalAndOr = Semiring{AndOr, Bool}

add_alg(::Type{AndOr}, a::Bool, b::Bool) = a || b
mul_alg(::Type{AndOr}, a::Bool, b::Bool) = a && b

zero_alg(::Type{AndOr}, ::Type{Bool}) = false
one_alg(::Type{AndOr}, ::Type{Bool}) = true

ldiv_alg(::Type{AndOr}, a::Bool, b::Bool) = b || !a
rdiv_alg(::Type{AndOr}, b::Bool, a::Bool) = ldiv_alg(AndOr, a, b)

leq_alg(::Type{AndOr}, a::Bool, b::Bool) = a <= b
lt_alg(::Type{AndOr}, a::Bool, b::Bool) = a < b

# --------------- #
# other operators #
# --------------- #

function Base.isless(a::TropicalAndOr, b::TropicalAndOr)
    return a < b
end
