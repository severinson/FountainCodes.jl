doc"For a monotonically increasing function f, find a value, lower <= x <=
upper, such that f(x) == target."
function numinv(f, target::Float64, lower::Float64=0.0, upper::Float64=Inf)
    f_wrapped = x -> bounds_wrap(x, f, lower, upper)
    lower2, upper2 = find_limits(f_wrapped, target)
    lower = max(lower, lower2)
    upper = min(upper, upper2)
    while (upper - lower) > 1
        x = Int64(floor(lower + (upper-lower)/2))
        v = f(x)
        if v > target
            upper = x
        elseif v < target
            lower = x
        else
            return x
        end
    end
    return lower
end

function find_limits(f, target::Float64)
    upper = 1
    lower = 0
    while f(upper) < target
        lower = upper
        upper *= 2
    end
    return lower, upper
end

function bounds_wrap(x::Int64, f, lower::Float64, upper=Float64)
    if x < lower
        return -Inf
    end
    if x > upper
        return Inf
    end
    return f(x)
end
