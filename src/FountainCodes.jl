# Copyright 2019 Albin Severinson
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

module FountainCodes

using Random, StatsBase, Statistics, Distributions
using LinearAlgebra, SparseArrays
using Primes, DataStructures

export CoefficientType, Binary, NonBinary, Code

# type system
abstract type AbstractErasureCode end

# For testing
export TestConstraints, subtract!, rand_nonzero
struct TestConstraints{Tc,Tv} <: AbstractVector{Tv}
    A::Matrix{Tc}
    Vs::Vector{Tv}
    function TestConstraints(A, Vs)
        _, n = size(A)
        length(Vs) == n || throw(DimensionMismatch("Vs has dimension $length(Vs)), but A has dimensions $(size(A))"))
        new{eltype(A),eltype(Vs)}(Matrix(A), deepcopy(Vs))
    end
end

Base.length(tc::TestConstraints) = length(tc.Vs)
Base.size(tc::TestConstraints) = (length(tc),)
Base.getindex(tc::TestConstraints, args...) = getindex(tc.Vs, args...)

function subtract!(tc::TestConstraints; coef, rpi_src, rpi_dst)
    tc.A[:, rpi_dst] .-= coef .* view(tc.A, :, rpi_src)
    tc.Vs[rpi_dst] -= coef * tc.Vs[rpi_src]
    return
end

function rand_nonzero(rng::AbstractRNG, T::Type)
    rv = rand(rng, T)
    while iszero(rv)
        rv = rand(rng, T)    
    end
    rv
end

rand_nonzero(rng::AbstractRNG, T::Type{<:AbstractFloat}) = randn(rng, T)

function rand_nonzero(rng::AbstractRNG, T::Type, dims...)
    rv = zeros(T, dims...)
    for i in 1:length(rv)
        rv[i] = rand_nonzero(rng, T)
    end
    rv
end

export dimension, constraint_matrix, generator_matrix

include("GF256.jl")
include("QMatrix.jl")
include("Decode.jl")
include("Gray.jl")
include("R10Tables.jl")
include("R10.jl")
include("LT.jl")
# include("RQTables.jl")
# include("RQ.jl")
# include("LDPC.jl")

end # module
