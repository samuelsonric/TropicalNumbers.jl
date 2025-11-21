"""
    CountingTropical{T, C} <: Number

Counting tropical number type is also a semiring algebra.
It is tropical algebra with one extra field for counting, it is introduced in [arXiv:2008.06888](https://arxiv.org/abs/2008.06888).

Example
-------------------------
```jldoctest; setup=:(using TropicalNumbers)
julia> CountingTropical(1.0, 5.0) + CountingTropical(3.0, 2.0)
(3.0, 2.0)

julia> CountingTropical(1.0, 5.0) * CountingTropical(3.0, 2.0)
(4.0, 10.0)

julia> one(CountingTropicalF64)
(0.0, 1.0)

julia> zero(CountingTropicalF64)
(-Inf, 0.0)
```
"""
struct CountingTropical{T,C} <: Number
    n::T
    c::C
end
CountingTropical(n::T) where {T <: Real} = CountingTropical{T}(n)
CountingTropical{T}(n) where {T} = CountingTropical{T, T}(n)
CountingTropical{T, C}(n) where {T, C} = CountingTropical{T, C}(n, one(C))
CountingTropical{T, C}(a::Tropical) where {T, C} = CountingTropical{T, C}(a.n)
CountingTropical{T, C}(a::CountingTropical{T, C}) where {T, C} = a
CountingTropical{T, C}(a::CountingTropical) where {T, C} = CountingTropical{T, C}(a.n, a.c)

Base.:*(a::CountingTropical, b::CountingTropical) = CountingTropical(a.n + b.n, a.c * b.c)
Base.:^(a::CountingTropical, b::Real) = CountingTropical(a.n * b, a.c ^ b)
Base.:^(a::CountingTropical, b::Integer) = CountingTropical(a.n * b, a.c ^ b)

function Base.:+(a::CountingTropical, b::CountingTropical)
    if a.n > b.n
        n = a.n
        c = a.c
    elseif a.n == b.n
        n = a.n
        c =  a.c + b.c
    else
        n = b.n
        c = b.c
    end

    CountingTropical(n, c)
end

Base.inv(a::CountingTropical) = CountingTropical(-a.n, a.c)
Base.typemin(::Type{CountingTropical{T, C}}) where {T, C} = CountingTropical(neginf(T), zero(C))

Base.zero(::Type{CountingTropical{T, C}}) where {T, C} = typemin(CountingTropical{T, C})
Base.zero(::Type{CountingTropical{T}}) where {T} = zero(CountingTropical{T, T})
Base.zero(::T) where {T <: CountingTropical} = zero(T)

Base.one(::Type{CountingTropical{T, C}}) where {T, C} = CountingTropical(zero(T), one(C))
Base.one(::Type{CountingTropical{T}}) where {T} = one(CountingTropical{T,T})
Base.one(::T) where {T <: CountingTropical} = one(T)

Base.isapprox(a::CountingTropical, b::CountingTropical; kw...) = isapprox(a.n, b.n; kw...) && isapprox(a.c, b.c; kw...)
Base.isapprox(x::AbstractArray{<:CountingTropical}, y::AbstractArray{<:CountingTropical}; kw...) = all(isapprox.(x, y; kw...))

Base.show(io::IO, t::CountingTropical) = Base.print(io, (t.n, t.c))
Base.promote_rule(::Type{CountingTropical{T, C}}, b::Type{CountingTropical{U, D}}) where {T, U, C, D} = CountingTropical{promote_type(T, U), promote_type(C, D)}

content(a::CountingTropical) = a.n
content(::Type{<:CountingTropical{T}}) where {T} = T
Base.isnan(a::CountingTropical) = isnan(a.n)
Base.isinf(a::CountingTropical) = isinf(a.n)
Base.:(==)(a::CountingTropical, b::CountingTropical) = a.n == b.n
Base.:<=(a::CountingTropical, b::CountingTropical) = a.n <= b.n
Base.:<(a::CountingTropical, b::CountingTropical) = a.n < b.n
Base.isless(a::CountingTropical, b::CountingTropical) = a < b
Base.:*(a::CountingTropical, b::Bool) = ifelse(b, a, zero(a))
Base.:*(a::Bool, b::CountingTropical) = b * a
Base.:/(b::CountingTropical, a::Bool) = ifelse(a, b, b / zero(a))
Base.:/(b::Bool, a::CountingTropical) = ifelse(b, inv(a), zero(a))
