"Soliton distribution."
struct Soliton
    num_symbols::Int64 # number of input symbols
    mode::Int64 # robust component spike location
    delta::Float64 # interpreted as the decoder failure probability
    R::Float64
    c::Float64
    beta::Float64 # normalization constant
    function Soliton(num_symbols::Int64, mode::Int64, delta::Float64)
        if !(0 < delta < 1)
            error("delta must be 0 < delta < 1")
        end
        if !(1 <= mode <= num_symbols)
            error("mode must be 1 <= mode <= num_symbols")
        end
        R = num_symbols / mode
        c = R / log(num_symbols / delta) / sqrt(num_symbols)
        beta = sum(
            tau(num_symbols, mode, delta, R, i) + rho(num_symbols, i)
            for i=1:num_symbols
        )
        soliton = new(
            num_symbols,
            mode,
            delta,
            R,
            c,
            beta,
        )
    end
end

"Robust component of the Soliton distribution."
function tau(num_symbols::Int64, mode::Int64, delta::Float64, R::Float64, i::Int64)
    if i < mode
        return 1 / (i * mode)
    elseif i == mode
        return log(R / delta) / mode
    elseif i <= num_symbols
        return 0
    else
        error("i must be less than or equal to num_symbols")
    end
end

"Robust component of the Soliton distribution."
function tau(f::Soliton, i::Int64)
    if i < f.mode
        return 1 / (i * f.mode)
    elseif i == f.mode
        return log(f.R / f.delta) / f.mode
    elseif i <= f.num_symbols
        return 0
    else
        error("i must be less than or equal to f.num_symbols")
    end
end

"Ideal component of the Soliton distribution."
function rho(num_symbols::Int64, i::Int64)
    if i == 1
        return 1 / num_symbols
    elseif i <= num_symbols
        return 1 / (i * (i - 1))
    else
        error("i must be less than or equal to num_symbols")
    end
end

"Ideal component of the Soliton distribution."
function rho(f::Soliton, i::Int64)
    if i == 1
        return 1 / f.num_symbols
    elseif i <= f.num_symbols
        return 1 / (i * (i - 1))
    else
        error("i must be less than or equal to f.num_symbols")
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
    return numinv(x->cdf(f, x), v, 1.0, Float64(f.num_symbols))
end

"Draw a random sample from the distribution."
function sample(f::Soliton)
    return icdf(f, rand())
end
