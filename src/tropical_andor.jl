

"""
    TropicalAndOr <: AbstractSimpleSemiring{Bool}

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
trueₜ

julia> TropicalAndOr(true) * TropicalAndOr(false)
falseₜ

julia> one(TropicalAndOr)
trueₜ

julia> zero(TropicalAndOr)
falseₜ
```
"""
struct TropicalAndOr <: AbstractSimpleSemiring{Bool}
    n::Bool
    TropicalAndOr(x::T) where T <: Bool = new(x)
end

content(a::TropicalAndOr) = a.n

inf(a::TropicalAndOr, b::TropicalAndOr) = TropicalAndOr(a.n && b.n)
sup(a::TropicalAndOr, b::TropicalAndOr) = TropicalAndOr(a.n || b.n)

Base.typemin(::Type{TropicalAndOr}) = TropicalAndOr(false)
Base.typemax(::Type{TropicalAndOr}) = TropicalAndOr(true)

Base.:\(a::TropicalAndOr, b::TropicalAndOr) = b / a
Base.:/(b::TropicalAndOr, a::TropicalAndOr) = TropicalAndOr(b.n || !a.n)
Base.div(b::TropicalAndOr, a::TropicalAndOr) = b / a

# Base.inv(a::TropicalAndOr) = TropicalAndOr(!a.n)
# invint(a::TropicalAndOr) = inv(a)
