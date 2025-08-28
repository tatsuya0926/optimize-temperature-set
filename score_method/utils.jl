module Utils

function set_temperature_ladder(
    β_min, 
    β_max; 
    M=15, 
    method=:inverse_linear
)
    if method == :geometric
        if β_min == 0
            β_min = 1e-8  # Avoid division by zero
        end 
        R = (β_max / β_min)^(1/(M - 1))
        temperatures = [1/β_max * R^(i - 1) for i in 1:M]
        return sort(1 ./ temperatures)
    elseif method == :inverse_linear
        return [β_min + (β_max - β_min) * (i - 1)/(M - 1) for i in 1:M]
    else
        @error "Unknown method for temperature ladder"
    end
end
    export set_temperature_ladder
end