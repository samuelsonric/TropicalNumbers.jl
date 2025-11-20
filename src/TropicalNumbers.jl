module TropicalNumbers

using Random
using Random: SamplerType

export content, neginf, posinf, inf, fli, fri, inf_fast, ldiv_fast, ∧
export TropicalTypes, AbstractSemiring

export TropicalAndOr
export Tropical, TropicalF64, TropicalF32, TropicalF16, TropicalI64, TropicalI32, TropicalI16
export TropicalMaxPlus, TropicalMaxPlusF64, TropicalMaxPlusF32, TropicalMaxPlusF16, TropicalMaxPlusI64, TropicalMaxPlusI32, TropicalMaxPlusI16
export TropicalMinPlus, TropicalMinPlusF64, TropicalMinPlusF32, TropicalMinPlusF16, TropicalMinPlusI64, TropicalMinPlusI32, TropicalMinPlusI16
export TropicalMaxMul, TropicalMaxMulF64, TropicalMaxMulF32, TropicalMaxMulF16, TropicalMaxMulI64, TropicalMaxMulI32, TropicalMaxMulI16
export TropicalMaxMin, TropicalMaxMinF64, TropicalMaxMinF32, TropicalMaxMinF16, TropicalMaxMinI64, TropicalMaxMinI32, TropicalMaxMinI16
export TropicalBitwise, TropicalBitwiseI64, TropicalBitwiseI32, TropicalBitwiseI16
export CountingTropical, CountingTropicalF16, CountingTropicalF32, CountingTropicalF64, CountingTropicalI16, CountingTropicalI32, CountingTropicalI64

neginf(::Type{T}) where {T} = typemin(T)
neginf(::Type{T}) where {T <: Integer} = T(-999999)
neginf(::Type{Int16}) = Int16(-16384)
neginf(::Type{Int8}) = Int8(-64)
posinf(::Type{T}) where {T} = -neginf(T)

include("abstract_semiring_algebra.jl")
include("semiring.jl")
include("tropical_maxplus.jl")
include("tropical_andor.jl")
include("tropical_minplus.jl")
include("tropical_maxmul.jl")
include("tropical_maxmin.jl")
include("tropical_bitwise.jl")
include("counting_tropical.jl")

const TropicalTypes{T} = Union{CountingTropical{T}, Tropical{T}}
const ∧ = inf

# alias
# defining constants like `TropicalF64`.
for NBIT in [16, 32, 64]
    @eval const $(Symbol(:Tropical, :F, NBIT)) = Tropical{$(Symbol(:Float, NBIT))}
    @eval const $(Symbol(:TropicalMaxPlus, :F, NBIT)) = TropicalMaxPlus{$(Symbol(:Float, NBIT))}
    @eval const $(Symbol(:TropicalMinPlus, :F, NBIT)) = TropicalMinPlus{$(Symbol(:Float, NBIT))}
    @eval const $(Symbol(:TropicalMaxMul, :F, NBIT)) = TropicalMaxMul{$(Symbol(:Float, NBIT))}
    @eval const $(Symbol(:TropicalMaxMin, :F, NBIT)) = TropicalMaxMin{$(Symbol(:Float, NBIT))}
    @eval const $(Symbol(:CountingTropical, :F, NBIT)) = CountingTropical{$(Symbol(:Float, NBIT)),$(Symbol(:Float, NBIT))}

    @eval const $(Symbol(:Tropical, :I, NBIT)) = Tropical{$(Symbol(:Int, NBIT))}
    @eval const $(Symbol(:TropicalMaxPlus, :I, NBIT)) = TropicalMaxPlus{$(Symbol(:Int, NBIT))}
    @eval const $(Symbol(:TropicalMinPlus, :I, NBIT)) = TropicalMinPlus{$(Symbol(:Int, NBIT))}
    @eval const $(Symbol(:TropicalMaxMul, :I, NBIT)) = TropicalMaxMul{$(Symbol(:Int, NBIT))}
    @eval const $(Symbol(:TropicalMaxMin, :I, NBIT)) = TropicalMaxMin{$(Symbol(:Int, NBIT))}
    @eval const $(Symbol(:TropicalBitwise, :I, NBIT)) = TropicalBitwise{$(Symbol(:Int, NBIT))}
    @eval const $(Symbol(:CountingTropical, :I, NBIT)) = CountingTropical{$(Symbol(:Int, NBIT)),$(Symbol(:Float, NBIT))}
end

end # module
