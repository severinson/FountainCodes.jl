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

using StatsBase, Statistics, Distributions, Random
using Primes, DataStructures, LinearAlgebra, SparseArrays
export CoefficientType, Binary, NonBinary, Code

# type system
abstract type AbstractErasureCode end

include("Numinv.jl")
include("Soliton.jl")
include("GF256.jl")
include("CodedMvNormal.jl")
include("Gray.jl")
include("QMatrix.jl")
include("DisjointSet.jl")
include("Decode.jl")
include("R10Tables.jl")
include("R10.jl")
include("LT.jl")
include("RQTables.jl")
include("RQ.jl")
include("LDPC.jl")

end # module
