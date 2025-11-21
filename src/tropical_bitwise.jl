struct Bitwise <: AbstractLatticeAlgebra end

"""
    TropicalBitwise{T} <: AbstractSemiring

`TropicalBitwise` is a semiring algebra that parallelizes the [`TropicalAndOr`](@ref) algebra,
It can be described by
* TropicalBitwise, (ℝ, |, &, 0, ~0).

It maps
* `+` to `|`
* `*` to `&`
* `0` to `0`
* `1` to `~0`

Example
-------------------------
```jldoctest; setup=:(using TropicalNumbers)
julia> TropicalBitwise(1) + TropicalBitwise(3)
3

julia> TropicalBitwise(1) * TropicalBitwise(3)
1

julia> zero(TropicalBitwiseI64)
0

julia> one(TropicalBitwiseI64)
-1
```
"""
const TropicalBitwise = Semiring{Bitwise}
const TropicalAndOr = TropicalBitwise{Bool}

add_alg(::Type{Bitwise}, a, b) = a | b
mul_alg(::Type{Bitwise}, a, b) = a & b
not_alg(::Type{Bitwise}, a) = ~a

zero_alg(::Type{Bitwise}, ::Type{T}) where {T} = zero(T)
one_alg(::Type{Bitwise}, ::Type{T}) where {T} = ~zero(T)

ldiv_alg(::Type{Bitwise}, a, b) = b | ~a

# --------------- #
# other operators #
# --------------- #

function Base.isless(a::TropicalAndOr, b::TropicalAndOr)
    return a < b
end
