using Statistics


function onepoint_extrema(x)
    maxima = Int64[]
    sizehint!(maxima, 32)
    minima = Int64[]
    sizehint!(minima, 32)

    for i in 2:(length(x)-1)
        if x[i+1] < x[i] > x[i-1]
            push!(maxima, i)
        elseif x[i+1] > x[i] < x[i-1]
            push!(minima, i)
        end
    end
    maxima, minima
end

function extract_baseline!(signal, baseline)
    can_be_decomposed_further = nothing
    maxima, minima = onepoint_extrema(signal)
    print("Max and min lengths", " ", length(maxima), " ", length(minima))
    argextrema = zeros(Int64, length(maxima) + length(minima) + 2)
    argextrema[begin] = 1
    argextrema[end] = length(signal)

    if length(maxima) > 0 && length(minima) > 0
        can_be_decomposed_further = true
        argextrema[2:end-1] = [maxima minima]'[:]
        if maxima[1] < minima[1]
            argextrema[2:2:end-2] = maxima
            argextrema[3:2:end-1] = minima
        else
            argextrema[2:2:end-2] = minima
            argextrema[3:2:end-1] = maxima
        end
    else
        can_be_decomposed_further = false
        return can_be_decomposed_further
    end

    baseline_knots = similar(argextrema, Float64)
    baseline_knots[1] = mean(signal[1:2])
    baseline_knots[end] = mean(signal[end-1:end])

    for k in range(2, length(argextrema) - 1)
        baseline_knots[k] = 0.5 * (signal[argextrema[k-1]] + 
                             (argextrema[k] - argextrema[k - 1]) / (argextrema[k + 1] - argextrema[k - 1]) * 
                             (signal[argextrema[k + 1]] - signal[argextrema[k-1]])) + 
        0.5 * signal[argextrema[k]]
    end

    for k in range(1, length(argextrema) - 1)
        @. baseline[argextrema[k]:argextrema[k + 1]] = baseline_knots[k] + 
                        (baseline_knots[k + 1] - baseline_knots[k]) / (signal[argextrema[k + 1]] - signal[argextrema[k]]) * 
                            (signal[argextrema[k]:argextrema[k + 1]] - signal[argextrema[k]])
    end
    can_be_decomposed_further
end

function itd(signal, num_components)
    signalcopy = similar(signal)
    signalcopy[:] = signal
    can_be_decomposed_further = true
    baselines = []
    proper_rotations = []
    baseline = similar(signal)
    for i in range(1, num_components)
        println(i)
        can_be_decomposed_further = extract_baseline!(signalcopy, baseline)
        push!(baselines, copy(baseline))
        push!(proper_rotations, copy(signalcopy - baseline))
        signalcopy = copy(baseline)
        println(can_be_decomposed_further)
        if !can_be_decomposed_further
            break
        end
    end

    baselines, proper_rotations
end