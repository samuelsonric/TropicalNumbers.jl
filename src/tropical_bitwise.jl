"""
    TropicalBitwise{T} <: AbstractSimpleSemiring{T}

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
3ₛ

julia> TropicalBitwise(1) * TropicalBitwise(3)
1ₛ

julia> zero(TropicalBitwiseI64)
0ₛ

julia> one(TropicalBitwiseI64)
-1ₛ
```
"""
struct TropicalBitwise{T} <: AbstractSimpleSemiring{T}
    n::T
end

function TropicalBitwise(a::TropicalBitwise)
    return TropicalBitwise(a.n)
end

function TropicalBitwise{T}(a::TropicalBitwise) where {T}
    return TropicalBitwise{T}(a.n)
end

function Base.promote_rule(::Type{TropicalBitwise{U}}, ::Type{TropicalBitwise{V}}) where {U, V}
    W = promote_type(U, V)
    return TropicalBitwise{W}
end

function content(a::TropicalBitwise)
    return a.n
end

function inf(a::TropicalBitwise, b::TropicalBitwise)
    n = a.n & b.n
    return TropicalBitwise(n)
end

function sup(a::TropicalBitwise, b::TropicalBitwise)
    n = a.n | b.n
    return TropicalBitwise(n)
end

function Base.typemin(::Type{TropicalBitwise{T}}) where {T}
    n = zero(T)
    return TropicalBitwise(n)
end

function Base.typemax(::Type{TropicalBitwise{T}}) where {T}
    n = ~zero(T)
    return TropicalBitwise(n)
end

function Base.:\(a::TropicalBitwise, b::TropicalBitwise)
    return b / a
end

function Base.:/(b::TropicalBitwise, a::TropicalBitwise)
    return TropicalBitwise(b.n | ~a.n)
end

function Base.div(b::TropicalBitwise, a::TropicalBitwise)
    return b / a
end
