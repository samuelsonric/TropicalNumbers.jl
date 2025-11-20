abstract type AbstractSemiringAlgebra end

abstract type AbstractQuantaleAlgebra <: AbstractSemiringAlgebra end

abstract type AbstractLatticeAlgebra <: AbstractQuantaleAlgebra end

abstract type AbstractTropicalAlgebra <: AbstractQuantaleAlgebra end

# --------- #
# Semirings #
# --------- #

"""
    zero_alg(::Type{T}, A::Type) where {T <: AbstractSemiringAlgebra}
"""
zero_alg(::Type{T}, A::Type) where {T <: AbstractSemiringAlgebra}

"""
    one_alg(::Type{T}, A::Type) where {T <: AbstractSemiringAlgebra}
"""
one_alg(::Type{T}, A::Type) where {T <: AbstractSemiringAlgebra}

"""
    add_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}
"""
add_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}

"""
    add_fast_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}
"""
function add_fast_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}
    c = add_alg(T, a, b)
    return c
end

"""
    mul_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}
"""
mul_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}

"""
    mul_fast_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}
"""
function mul_fast_alg(::Type{T}, a, b) where {T <: AbstractSemiringAlgebra}
    c = mul_alg(T, a, b)
    return c
end

"""
    mul_add_alg(::Type{T}, a, b, c) where {T <: AbstractSemiringAlgebra}
"""
function mul_add_alg(::Type{T}, a, b, c) where {T <: AbstractSemiringAlgebra}
    return add_fast_alg(T, mul_fast_alg(T, a, b), c)
end

# --------- #
# Quantales #
# --------- #

"""
    typemax_alg(::Type{T}, A::Type) where {T <: AbstractQuantaleAlgebra}
"""
typemax_alg(::Type{T}, A::Type) where {T <: AbstractQuantaleAlgebra}

"""
    inf_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
"""
inf_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}

"""
    inf_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
"""
function inf_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    c = inf_alg(T, a, b)
    return c
end

"""
    ldiv_alg(::Type{T}, a, b) whre {T <: AbstractQuantaleAlgebra}
"""
ldiv_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}

"""
    ldiv_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
"""
function ldiv_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    c = ldiv_alg(T, a, b)
    return c
end

"""
    inf_ldiv_alg(::Type{T}, a, b, c) where {T <: AbstractQuantaleAlgebra}
"""
function inf_ldiv_alg(::Type{T}, a, b, c) where {T <: AbstractQuantaleAlgebra}
    return inf_fast_alg(T, ldiv_fast_alg(T, a, b), c)
end

"""
    rdiv_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
"""
rdiv_alg(::Type{T}, b, a) where {T <: AbstractQuantaleAlgebra}

"""
    rdiv_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
"""
function rdiv_fast_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    c = rdiv_alg(T, a, b)
    return c
end

"""
    inf_rdiv_alg(::Type{T}, b, a, c) where {T <: AbstractQuantaleAlgebra}
"""
function inf_rdiv_alg(::Type{T}, b, a, c) where {T <: AbstractQuantaleAlgebra}
    return inf_fast_alg(T, rdiv_fast_alg(T, b, a), c)
end

"""
    leq_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
"""
function leq_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    return add_alg(T, a, b) == b
end

"""
    lt_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
"""
function lt_alg(::Type{T}, a, b) where {T <: AbstractQuantaleAlgebra}
    return a != b && leq_alg(T, a, b)
end

# -------- #
# Lattices #
# -------- #

function typemax_alg(::Type{T}, A::Type) where {T <: AbstractLatticeAlgebra}
    return one_alg(T, A)
end

function inf_alg(::Type{T}, a, b) where {T <: AbstractLatticeAlgebra}
    return mul_alg(T, a, b)
end

function inf_fast_alg(::Type{T}, a, b) where {T <: AbstractLatticeAlgebra}
    return mul_fast_alg(T, a, b)
end

# -------- #
# Tropical #
# -------- #

"""
    exp_alg(::Type{T}, a, b) where {T <: AbstractTropicalAlgebra}
"""
exp_alg(::Type{T}, a, b) where {T <: AbstractTropicalAlgebra}

"""
    inv_alg(::Type{T}, a) where {T <: AbstractTropicalAlgebra}
"""
function inv_alg(::Type{T}, a) where {T <: AbstractTropicalAlgebra}
    return div_alg(T, one_alg(T, a), a)
end
