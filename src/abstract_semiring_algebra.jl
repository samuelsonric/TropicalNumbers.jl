abstract type AbstractSemiringAlgebra end

abstract type AbstractQuantaleAlgebra <: AbstractSemiringAlgebra end

abstract type AbstractCommutativeQuantaleAlgebra <: AbstractQuantaleAlgebra end

abstract type AbstractLatticeAlgebra <: AbstractCommutativeQuantaleAlgebra end

abstract type AbstractTropicalAlgebra <: AbstractCommutativeQuantaleAlgebra end

struct LatticeAlgebra{T <: AbstractQuantaleAlgebra} <: AbstractLatticeAlgebra end

# --------- #
# Semirings #
# --------- #

"""
    zero_alg(::Type{T}, A::Type) where {T <: AbstractSemiringAlgebra}
"""
zero_alg(::Type{T}, A::Type) where {T <: AbstractSemiringAlgebra}

function zero_alg(::Type{LatticeAlgebra{T}}, ::Type{A}) where {T <: AbstractQuantaleAlgebra, A}
    return zero_alg(T, A)
end

"""
    one_alg(::Type{T}, A::Type) where {T <: AbstractSemiringAlgebra}
"""
one_alg(::Type{T}, A::Type) where {T <: AbstractSemiringAlgebra}

function one_alg(::Type{LatticeAlgebra{T}}, ::Type{A}) where {T <: AbstractQuantaleAlgebra, A}
    return typemax_alg(T, A)
end

"""
    add_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}
"""
add_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}

function add_alg(::Type{LatticeAlgebra{T}}, a, b) where {T <: AbstractQuantaleAlgebra}
    return add_alg(T, a, b)
end

"""
    add_fast_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}
"""
function add_fast_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}
    return add_alg(T, a, b)
end

function add_fast_alg(::Type{LatticeAlgebra{T}}, a, b) where {T <: AbstractQuantaleAlgebra}
    return add_fast_alg(T, a, b)
end

"""
    mul_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}
"""
mul_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}

function mul_alg(::Type{LatticeAlgebra{T}}, a, b) where {T <: AbstractQuantaleAlgebra}
    return inf_alg(T, a, b)
end

"""
    mul_fast_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}
"""
function mul_fast_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}
    return mul_alg(T, a, b)
end

function mul_fast_alg(::Type{LatticeAlgebra{T}}, a, b) where {T <: AbstractQuantaleAlgebra}
    return inf_fast_alg(T, a, b)
end

"""
    mul_add_alg(::Type{T}, a, b, c) where {T <: AbstractSemiringAlgebra}
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
"""
typemax_alg(::Type{T}, A::Type) where {T <: AbstractQuantaleAlgebra}

function typemax_alg(::Type{T}, ::Type{A}) where {T <: AbstractLatticeAlgebra, A}
    return one_alg(T, A)
end

"""
    inf_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
"""
inf_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}

function inf_alg(::Type{T}, a, b) where {T <: AbstractLatticeAlgebra}
    return mul_alg(T, a, b)
end

"""
    inf_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
"""
function inf_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    return inf_alg(T, a, b)
end

function inf_fast_alg(::Type{T}, a, b) where {T <: AbstractLatticeAlgebra}
    return mul_fast_alg(T, a, b)
end

"""
    inf_add_alg(::Type{T}, a, b, c) where {T <: AbstractQuantaleAlgebra}
"""
function inf_add_alg(::Type{T}, a, b, c) where {T <: AbstractQuantaleAlgebra}
    return add_fast_alg(T, inf_fast_alg(T, a, b), c)
end

"""
    ldiv_alg(::Type{T}, a, b) whre {T <: AbstractQuantaleAlgebra}
"""
ldiv_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}

function ldiv_alg(::Type{LatticeAlgebra{T}}, a, b) where {T <: AbstractQuantaleAlgebra}
    return imp_alg(T, a, b)
end

"""
    ldiv_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
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
"""
function inf_ldiv_alg(::Type{T}, a, b, c) where {T <: AbstractQuantaleAlgebra}
    return inf_fast_alg(T, ldiv_fast_alg(T, a, b), c)
end

function inf_ldiv_alg(::Type{LatticeAlgebra{T}}, a, b, c) where {T <: AbstractQuantaleAlgebra}
    return inf_imp_alg(T, a, b, c)
end

"""
    rdiv_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
"""
rdiv_alg(::Type{T}, b, a) where {T <: AbstractQuantaleAlgebra}

function rdiv_alg(::Type{T}, b, a) where {T <: AbstractCommutativeQuantaleAlgebra}
    return ldiv_alg(T, a, b)
end

"""
    rdiv_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
"""
function rdiv_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    return rdiv_alg(T, a, b)
end

function rdiv_fast_alg(::Type{T}, b, a) where {T <: AbstractCommutativeQuantaleAlgebra}
    return ldiv_fast_alg(T, a, b)
end

"""
    inf_rdiv_alg(::Type{T}, b, a, c) where {T <: AbstractQuantaleAlgebra}
"""
function inf_rdiv_alg(::Type{T}, b, a, c) where {T <: AbstractQuantaleAlgebra}
    return inf_fast_alg(T, rdiv_fast_alg(T, b, a), c)
end

"""
    imp_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
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
"""
function imp_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    return imp_alg(T, a, b)
end

function imp_fast_alg(::Type{T}, a, b) where {T <: AbstractLatticeAlgebra}
    return ldiv_fast_alg(T, a, b)
end

"""
    inf_imp_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
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
"""
function inv_alg(::Type{T}, a::A) where {T <: AbstractQuantaleAlgebra, A}
    return div_alg(T, one_alg(T, A), a)
end

function inv_alg(::Type{T}, a::A) where {T <: AbstractLatticeAlgebra, A}
    return one_alg(T, A)
end

"""
    not_alg(::Type{T}, a) where {T <: AbstractQuantaleAlgebra}
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
"""
function leq_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    return add_alg(T, a, b) == b
end

function leq_alg(::Type{LatticeAlgebra{T}}, a, b) where {T <: AbstractQuantaleAlgebra}
    return leq_alg(T, a, b)
end

"""
    lt_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
"""
function lt_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    return a != b && leq_alg(T, a, b)
end

function lt_alg(::Type{LatticeAlgebra{T}}, a, b) where {T <: AbstractQuantaleAlgebra}
    return lt_alg(T, a, b)
end

"""
    exp_alg(::Type{T}, a, b) where {T <: AbstractTropicalAlgebra}
"""
exp_alg(::Type{T}, a, b) where {T <: AbstractTropicalAlgebra}
