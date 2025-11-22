"""
    AbstractSemiringAlgebra

A [`semiring`](https://en.wikipedia.org/wiki/Semiring) is a quintuple (R, +, ×, 0, 1), where

  - (R, +, 0) is a commutative monoid
  - (R, ×, 1) is a monoid
  - multiplication (×) distributes over addition (+)
  - zero (0) is absorbing
"""
abstract type AbstractSemiringAlgebra end

"""
    AbstractQuantaleAlgebra <: AbstractSemiringAlgebra

A [`semiring`](https://en.wikipedia.org/wiki/Semiring) (R, +, ×, 0, 1) is called a 
[`quantale`](https://en.wikipedia.org/wiki/Quantale) if it is additionally a complete
lattice whose supremum operation coincides with addition (+).
"""
abstract type AbstractQuantaleAlgebra <: AbstractSemiringAlgebra end

"""
    AbstractCommutativeQuantaleAlgebra <: AbstractQuantaleAlgebra end

A quantale (R, +, ×, 0, 1) is called commutative if the multiplication operation (×)
commutates.
"""
abstract type AbstractCommutativeQuantaleAlgebra <: AbstractQuantaleAlgebra end

"""
    AbstractLatticeAlgebra <: AbstractCommutativeQuantaleAlgebra

A quantale (R, +, ×, 0, 1) is called Cartesian, if its infimum operation
coincides with multiplcation (×).
"""
abstract type AbstractLatticeAlgebra <: AbstractCommutativeQuantaleAlgebra end

"""
    AbstractTropicalAlgebra <: AbstractCommutativeQuantaleAlgebra

The tropical semirings

  - ([-∞, +∞], ∧, +, +∞, 0)
  - ([-∞, +∞], ∨, +, -∞, 0)
  - ([0, +∞], ∨, *, 0, 0)

are commutative quantales with a well-defined exponentiation operation.
"""
abstract type AbstractTropicalAlgebra <: AbstractCommutativeQuantaleAlgebra end

"""
    LatticeAlgebra{T <: AbstractQuantaleAlgebra} <: AbstractLatticeAlgebra

Transform a quantale into a Cartesian quantale by replacing multiplication (×)
with the infimum operation.
"""
struct LatticeAlgebra{T <: AbstractQuantaleAlgebra} <: AbstractLatticeAlgebra end

# --------- #
# Semirings #
# --------- #

"""
    zero_alg(::Type{T}, A::Type) where {T <: AbstractSemiringAlgebra}

Get an additive identity 0 of type `A`.
"""
zero_alg(::Type{T}, A::Type) where {T <: AbstractSemiringAlgebra}

function zero_alg(::Type{LatticeAlgebra{T}}, ::Type{A}) where {T <: AbstractQuantaleAlgebra, A}
    return zero_alg(T, A)
end

"""
    one_alg(::Type{T}, A::Type) where {T <: AbstractSemiringAlgebra}

Get a multiplicative identity 1 of type `A`.
"""
one_alg(::Type{T}, A::Type) where {T <: AbstractSemiringAlgebra}

function one_alg(::Type{LatticeAlgebra{T}}, ::Type{A}) where {T <: AbstractQuantaleAlgebra, A}
    return typemax_alg(T, A)
end

"""
    add_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}

Compute the sum a + b.
"""
add_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}

function add_alg(::Type{LatticeAlgebra{T}}, a, b) where {T <: AbstractQuantaleAlgebra}
    return add_alg(T, a, b)
end

"""
    add_fast_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}

Compute the sum a + b.
"""
function add_fast_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}
    return add_alg(T, a, b)
end

function add_fast_alg(::Type{LatticeAlgebra{T}}, a, b) where {T <: AbstractQuantaleAlgebra}
    return add_fast_alg(T, a, b)
end

"""
    mul_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}

Compute the product a × b.
"""
mul_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}

function mul_alg(::Type{LatticeAlgebra{T}}, a, b) where {T <: AbstractQuantaleAlgebra}
    return inf_alg(T, a, b)
end

"""
    mul_fast_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}

Compute the product a × b.
"""
function mul_fast_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}
    return mul_alg(T, a, b)
end

function mul_fast_alg(::Type{LatticeAlgebra{T}}, a, b) where {T <: AbstractQuantaleAlgebra}
    return inf_fast_alg(T, a, b)
end

"""
    mul_add_alg(::Type{T}, a, b, c) where {T <: AbstractSemiringAlgebra}

Compute (a × b) + c.
"""
function mul_add_alg(::Type{T}, a, b, c) where {T <: AbstractSemiringAlgebra}
    return add_fast_alg(T, mul_fast_alg(T, a, b), c)
end

function mul_add_alg(::Type{LatticeAlgebra{T}}, a, b, c) where {T <: AbstractQuantaleAlgebra}
    return inf_add_alg(T, a, b, c)
end


# --------- #
# Quantales #
# --------- #

"""
    typemax_alg(::Type{T}, A::Type) where {T <: AbstractQuantaleAlgebra}

Get a top element ⊤ of type `A`.
"""
typemax_alg(::Type{T}, A::Type) where {T <: AbstractQuantaleAlgebra}

function typemax_alg(::Type{T}, ::Type{A}) where {T <: AbstractLatticeAlgebra, A}
    return one_alg(T, A)
end

"""
    inf_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}

Compute the infimum a ∧ b.
"""
inf_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}

function inf_alg(::Type{T}, a, b) where {T <: AbstractLatticeAlgebra}
    return mul_alg(T, a, b)
end

"""
    inf_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}

Compute infimum a ∧ b.
"""
function inf_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    return inf_alg(T, a, b)
end

function inf_fast_alg(::Type{T}, a, b) where {T <: AbstractLatticeAlgebra}
    return mul_fast_alg(T, a, b)
end

"""
    inf_add_alg(::Type{T}, a, b, c) where {T <: AbstractQuantaleAlgebra}

Compute (a ∧ b) + c.
"""
function inf_add_alg(::Type{T}, a, b, c) where {T <: AbstractQuantaleAlgebra}
    return add_fast_alg(T, inf_fast_alg(T, a, b), c)
end

"""
    ldiv_alg(::Type{T}, a, b) whre {T <: AbstractQuantaleAlgebra}

Compute the residual a \\ b.
"""
ldiv_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}

function ldiv_alg(::Type{LatticeAlgebra{T}}, a, b) where {T <: AbstractQuantaleAlgebra}
    return imp_alg(T, a, b)
end

"""
    ldiv_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}

Compute the residual a \\ b.
"""
function ldiv_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    return ldiv_alg(T, a, b)
end

function ldiv_fast_alg(::Type{LatticeAlgebra{T}}, a, b) where {T <: AbstractQuantaleAlgebra}
    return imp_fast_alg(T, a, b)
end

function inf_rdiv_alg(::Type{T}, b, a, c) where {T <: AbstractCommutativeQuantaleAlgebra}
    return inf_ldiv_alg(T, a, b, c)
end

"""
    inf_ldiv_alg(::Type{T}, a, b, c) where {T <: AbstractQuantaleAlgebra}

Compute (a \\ b) ∧ c.
"""
function inf_ldiv_alg(::Type{T}, a, b, c) where {T <: AbstractQuantaleAlgebra}
    return inf_fast_alg(T, ldiv_fast_alg(T, a, b), c)
end

function inf_ldiv_alg(::Type{LatticeAlgebra{T}}, a, b, c) where {T <: AbstractQuantaleAlgebra}
    return inf_imp_alg(T, a, b, c)
end

"""
    rdiv_alg(::Type{T}, b, a) where {T <: AbstractQuantaleAlgebra}

Compute the residual b / a.
"""
rdiv_alg(::Type{T}, b, a) where {T <: AbstractQuantaleAlgebra}

function rdiv_alg(::Type{T}, b, a) where {T <: AbstractCommutativeQuantaleAlgebra}
    return ldiv_alg(T, a, b)
end

"""
    rdiv_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}

Compute the residual b / a.
"""
function rdiv_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    return rdiv_alg(T, a, b)
end

function rdiv_fast_alg(::Type{T}, b, a) where {T <: AbstractCommutativeQuantaleAlgebra}
    return ldiv_fast_alg(T, a, b)
end

"""
    inf_rdiv_alg(::Type{T}, b, a, c) where {T <: AbstractQuantaleAlgebra}

Compute (b / a) ∧ c.
"""
function inf_rdiv_alg(::Type{T}, b, a, c) where {T <: AbstractQuantaleAlgebra}
    return inf_fast_alg(T, rdiv_fast_alg(T, b, a), c)
end

"""
    imp_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}

Compute the implication a → b.
"""
imp_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}

function imp_alg(::Type{T}, a, b) where {T <: AbstractLatticeAlgebra}
    return ldiv_alg(T, a, b)
end

function imp_alg(::Type{T}, a, b::A) where {T <: AbstractTropicalAlgebra, A}
    if leq_alg(T, a, b)
        c = typemax_alg(T, A)
    else
        c = b
    end

    return c
end

"""
    imp_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}

Compute the implication a → b.
"""
function imp_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    return imp_alg(T, a, b)
end

function imp_fast_alg(::Type{T}, a, b) where {T <: AbstractLatticeAlgebra}
    return ldiv_fast_alg(T, a, b)
end

"""
    inf_imp_alg(::Type{T}, a, b, c) where {T <: AbstractQuantaleAlgebra}

Compute (a → b) ∧ c.
"""
function inf_imp_alg(::Type{T}, a, b, c) where {T <: AbstractQuantaleAlgebra}
    return inf_fast_alg(T, imp_fast_alg(T, a, b), c)
end

function inf_imp_alg(::Type{T}, a, b, c) where {T <: AbstractLatticeAlgebra}
    return inf_ldiv_alg(T, a, b, c)
end

function inf_imp_alg(::Type{T}, a, b::A, c::A) where {T <: AbstractTropicalAlgebra, A}
    if leq_alg(T, a, b)
        d = c
    else
        d = inf_fast_alg(T, b, c)
    end

    return d
end

"""
    inv_alg(::Type{T}, a) where {T <: AbstractQuantaleAlgebra}

Compute the inverse 1 / a.
"""
function inv_alg(::Type{T}, a::A) where {T <: AbstractQuantaleAlgebra, A}
    return rdiv_alg(T, one_alg(T, A), a)
end

function inv_alg(::Type{T}, a::A) where {T <: AbstractLatticeAlgebra, A}
    return one_alg(T, A)
end

"""
    not_alg(::Type{T}, a) where {T <: AbstractQuantaleAlgebra}

Compute the psuedo-complement a → 0.
"""
function not_alg(::Type{T}, a::A) where {T <: AbstractQuantaleAlgebra, A}
    return imp_alg(T, a, zero_alg(T, A))
end

function not_alg(::Type{LatticeAlgebra{T}}, a) where {T <: AbstractQuantaleAlgebra}
    return not_alg(T, a)
end

function not_alg(::Type{T}, a::A) where {T <: AbstractTropicalAlgebra, A}
    if leq_alg(T, a, zero_alg(T, A))
        c = typemax_alg(T, A)
    else
        c = zero_alg(T, A)
    end

    return c
end

"""
    leq_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}

Evaluate a ≤ b.
"""
function leq_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    return add_alg(T, a, b) == b
end

function leq_alg(::Type{LatticeAlgebra{T}}, a, b) where {T <: AbstractQuantaleAlgebra}
    return leq_alg(T, a, b)
end

"""
    lt_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}

Evaluate a < b.
"""
function lt_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    return a != b && leq_alg(T, a, b)
end

function lt_alg(::Type{LatticeAlgebra{T}}, a, b) where {T <: AbstractQuantaleAlgebra}
    return lt_alg(T, a, b)
end

"""
    exp_alg(::Type{T}, a, b) where {T <: AbstractTropicalAlgebra}

Compute the exponent aᵇ.
"""
exp_alg(::Type{T}, a, b) where {T <: AbstractTropicalAlgebra}
