export CodedMvNormal, whiten

"""

Coded MvNormal, i.e., a multivariate normal random variable with a,
possibly, rank-deficit covariance matrix.

```
struct CodedMvNormal{Cov<:AbstractMatrix,Mean<:AbstractVector} <: AbstractMvNormal
    μ::Mean
    Σ::Cov
end
```

"""
struct CodedMvNormal{T<:Real,Cov<:AbstractMatrix, Mean<:AbstractVector} <: AbstractVector{T}
    μ::Mean
    Σ::Cov
end

## constructors ##
function CodedMvNormal(μ::AbstractVector{T}, Σ::AbstractMatrix{T}) where {T<:Real}
    size(Σ, 1) == length(μ) || throw(DimensionMismatch("The dimensions of μ and Σ are inconsistent."))
    isapprox(Σ, Σ') || throw(ArgumentError("Σ must be symmetric"))
    minimum(Diagonal(Σ)) >= 0 || throw(ArgumentError("variances must be non-negative"))
    CodedMvNormal{T,typeof(Σ), typeof(μ)}(μ, Σ)
end

function CodedMvNormal(μ::AbstractVector, Σ::AbstractMatrix)
    R = Base.promote_eltype(μ, Σ)
    CodedMvNormal(convert(AbstractArray{R}, μ), convert(AbstractArray{R}, Σ))
end

Base.show(io::IO, d::CodedMvNormal) = print(io, "CodedMvNormal($(d.μ), $(d.Σ))")

## basic statistics ##
Base.length(d::CodedMvNormal) = length(d.μ)
Base.size(d::CodedMvNormal) = (length(d.μ),)
Statistics.mean(d::CodedMvNormal) = d.μ
StatsBase.params(d::CodedMvNormal) = (d.μ, d.Σ)
@inline Distributions.partype(d::CodedMvNormal{T}) where {T<:Real} = T

Statistics.var(d::CodedMvNormal) = diag(d.Σ)
Statistics.cov(d::CodedMvNormal) = Matrix(d.Σ)

Distributions.invcov(d::CodedMvNormal) = Matrix(inv(d.Σ))
Distributions.logdetcov(d::CodedMvNormal) = logdet(d.Σ)

## indexing into the mean vector ###

Base.getindex(d::CodedMvNormal, inds...) = d.μ[inds...]
Base.setindex!(d::CodedMvNormal, v, inds...) = d.μ[inds...] = v

"""
    subtract!(d::CodedMvNormal, rpi::Integer, rpj::Integer, coef::Real)

Compute μ[rpj] -= coef*μ[rpi] in-place and update the covariance
accordingly.

"""
function subtract!(d::CodedMvNormal, rpi::Integer, rpj::Integer, coef::Real)
    if iszero(d[rpi]) return d[rpj] end
    d.μ[rpj] -= coef*d.μ[rpi]
    d.Σ[rpj, :] .-= coef.*d.Σ[rpi, :]
    d.Σ[:, rpj] .-= coef.*d.Σ[:, rpi]
    return d[rpj]
end

"""
    whiten(d::CodedMvNormal, is)

Return a linear transform M that when applied to d.μ[is] has the
effect of removing correlation and normalising variance, i.e.,
M*d.μ[is] has unit covariance.

"""
function whiten(d::CodedMvNormal, is)
    # covariance of the values with indices in rpis
    Σ = d.Σ[is, is]

    # compute the whitening transform
    Ω = svd!(Σ)
    r = searchsortedfirst(Ω.S, length(Ω.S) * eps(Float64)*Ω.S[1], rev=true) - 1
    M = Diagonal(1 ./ sqrt.(Ω.S))[1:r, :] * Ω.U'
    return M
end

"""
    whiten(d::CodedMvNormal)

Returns whiten(d, 1:length(d))

"""
function whiten(d::CodedMvNormal)
    return whiten(d, 1:length(d))
end
