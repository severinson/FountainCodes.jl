
export Soliton

"Soliton distribution."
struct Soliton
    K::Int64 # number of input symbols
    mode::Int64 # robust component spike location
    delta::Float64 # interpreted as the decoder failure probability
    R::Float64
    c::Float64
    beta::Float64 # normalization constant
    function Soliton(K::Int64, mode::Int64, delta::Float64)
        if !(0 < delta < 1)
            error("delta must be 0 < delta < 1")
        end
        if !(1 <= mode <= K)
            error("mode must be 1 <= mode <= K")
        end
        R = K / mode
        c = R / log(K / delta) / sqrt(K)
        beta = sum(
            tau(K, mode, delta, R, i) + rho(K, i)
            for i=1:K
        )
        soliton = new(
            K,
            mode,
            delta,
            R,
            c,
            beta,
        )
    end
end

Base.repr(s::Soliton) = "Soliton($(s.K), $(s.mode), $(s.delta))"

"Robust component of the Soliton distribution."
function tau(K::Int64, mode::Int64, delta::Float64, R::Float64, i::Int64)
    if i < mode
        return 1 / (i * mode)
    elseif i == mode
        return log(R / delta) / mode
    elseif i <= K
        return 0
    else
        error("i must be less than or equal to K")
    end
end

"Robust component of the Soliton distribution."
function tau(f::Soliton, i::Int64)
    if i < f.mode
        return 1 / (i * f.mode)
    elseif i == f.mode
        return log(f.R / f.delta) / f.mode
    elseif i <= f.K
        return 0
    else
        error("i must be less than or equal to f.K")
    end
end

"Ideal component of the Soliton distribution."
function rho(K::Int64, i::Int64)
    if i == 1
        return 1 / K
    elseif i <= K
        return 1 / (i * (i - 1))
    else
        error("i must be less than or equal to K")
    end
end

"Ideal component of the Soliton distribution."
function rho(f::Soliton, i::Int64)
    if i == 1
        return 1 / f.K
    elseif i <= f.K
        return 1 / (i * (i - 1))
    else
        error("i must be less than or equal to f.K")
    end
end

"Soliton distribution probability density function."
function pdf(f::Soliton, i::Int64)
    return (tau(f, i) + rho(f, i)) / f.beta
end

"Soliton distribution cumulative density function."
function cdf(f::Soliton, i::Int64)
    return sum(pdf(f, j) for j=1:i)
end

"Soliton distribution inverse cdf."
function icdf(f::Soliton, v::Float64)
    return numinv(x->cdf(f, x), v, 1.0, Float64(f.K))
end

"Draw a random sample from the distribution."
function sample(f::Soliton)
    return icdf(f, rand())
end
