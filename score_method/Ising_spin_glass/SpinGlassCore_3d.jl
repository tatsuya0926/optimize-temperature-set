module SpinGlassCore

using Random, LinearAlgebra, Statistics, StatsBase, Printf, Distributions

seed = 42
Random.seed!(seed)

function initial_state(N)
    state = ones(N, N, N)
    return state
end

function metropolis_sampler(config, β, N, Jh, Jv, Jz)
    sites = shuffle!([(i, j, k) for i in 1:N, j in 1:N, k in 1:N][:])
    for (i, j, k) in sites
        s = config[i, j, k]
        s_up = config[mod1(i - 1, N), j, k]
        s_down = config[mod1(i + 1, N), j, k]
        s_left = config[i, mod1(j - 1, N), k]
        s_right = config[i, mod1(j + 1, N), k]
        s_back = config[i, j, mod1(k - 1, N)]
        s_forward = config[i, j, mod1(k + 1, N)]
        
        local_field = (
            Jv[mod1(i-1, N), j, k] * s_up + Jv[i, j, k] * s_down + 
            Jh[i, mod1(j-1, N), k] * s_left + Jh[i, j, k] * s_right +
            Jz[i, j, mod1(k-1, N)] * s_back + Jz[i, j, k] * s_forward
        )
        
        cost = 2 * s * local_field
        if rand() < min(1.0, exp(-cost * β))
            config[i, j, k] *= -1
        end
    end
    
    return config
end

function calc_energy(config, N, Jh, Jv, Jz)
    energy = 0.0
    for i in 1:N, j in 1:N, k in 1:N
        s = config[i, j, k]
        energy -= s * (Jh[i, j, k] * config[i, mod1(j + 1, N), k] + Jv[i, j, k] * config[mod1(i + 1, N), j, k] + Jz[i, j, k] * config[i, j, mod1(k + 1, N)])
    end
    return energy
end

function score_method(
    N,
    β_replicas, 
    Jh, 
    Jv,
    Jz; 
    mcSteps=10^4, 
    eqSteps=10^3,
    exchange_interval=1,
    target=nothing
)
    M = length(β_replicas)
    config = initial_state(N)
    configs = [copy(config) for _ in 1:M]

    for _ in 1:eqSteps
        for r in 1:M
            configs[r] = metropolis_sampler(configs[r], β_replicas[r], N, Jh, Jv, Jz)
        end
    end

    exchange_probs_sum = zeros(M-1)
    exchange_counts = zeros(M-1)
    forward_energy_acceptance = zeros(M-1) # H(x(β_{i-1}))A_{i-1, i}
    backward_energy_acceptance = zeros(M-1) # H(x(β_{i+1}))A_{i, i+1}
    energy_sum = zeros(M)
    grad_estimates = zeros(M)
    step_count = 0

    for step in 1:mcSteps
        for r in 1:M
            configs[r] = metropolis_sampler(configs[r], β_replicas[r], N, Jh, Jv, Jz)
            energy_sum[r] += calc_energy(configs[r], N, Jh, Jv, Jz)
        end

        if step % exchange_interval == 0
            start = isodd(step_count) ? 2 : 1

            for r in start:2:(M - 1)
                E1 = calc_energy(configs[r], N, Jh, Jv, Jz)
                E2 = calc_energy(configs[r + 1], N, Jh, Jv, Jz)
                Δ = (β_replicas[r+1] - β_replicas[r]) * (E2 - E1)
                p = min(1.0, exp(Δ))
                if rand() < p
                    configs[r], configs[r+1] = configs[r+1], configs[r]
                end
                exchange_probs_sum[r] += p
                exchange_counts[r] += 1
                forward_energy_acceptance[r] += E1 * p
                backward_energy_acceptance[r] += E2 * p
            end
            step_count += 1
        end
    end
    energy_means = energy_sum ./ mcSteps
    exchange_probs_means = exchange_probs_sum ./ exchange_counts
    forward_energy_acceptance_mean = forward_energy_acceptance ./ exchange_counts
    backward_energy_acceptance_mean = backward_energy_acceptance ./ exchange_counts
    if target !== nothing
        A_bar = target
    else
        A_bar = mean(exchange_probs_means)
    end

    for r in 1:M
        if r == 1
            grad_estimates[1] = 2/(M-1) * (exchange_probs_means[1] - A_bar) * (- backward_energy_acceptance_mean[1] + energy_means[1] * exchange_probs_means[1])
        elseif r == M
            grad_estimates[M] = 2/(M-1) * (exchange_probs_means[M-1] - A_bar) * (- forward_energy_acceptance_mean[M-1] + energy_means[M] * exchange_probs_means[M-1])
        else
            grad_estimates[r] = 2/(M-1) * (
                (exchange_probs_means[r-1] - A_bar) * 
                    (-forward_energy_acceptance_mean[r-1] + energy_means[r] * exchange_probs_means[r-1]) + 
                (exchange_probs_means[r] - A_bar) * 
                    (-backward_energy_acceptance_mean[r] + energy_means[r] * exchange_probs_means[r])
                )
        end
    end

    return grad_estimates, exchange_probs_means
end

function calc_acceptance_and_rtt(
    N,
    β_replicas, 
    Jh, 
    Jv,
    Jz;
    mcSteps=10^4
)
    M = length(β_replicas)
    configs = [copy(initial_state(N)) for _ in 1:M]
    replica_indices = collect(1:M)
    temperature_indices = collect(1:M)

    exchange_count = 0
    exchange_probs_sum = zeros(M - 1)
    exchange_attempts = zeros(M - 1)

    replica_states = zeros(Int, M)
    round_trip_start_steps = zeros(Int, M)
    round_trip_times = [[] for _ in 1:M]

    for step in 1:mcSteps
        for i in 1:M
            current_replica_idx = replica_indices[i]
            configs[current_replica_idx] = metropolis_sampler(configs[current_replica_idx], β_replicas[i], N, Jh, Jv, Jz)
        end

        start_idx = isodd(exchange_count) ? 2 : 1
        for r in start_idx:2:(M - 1)
            rep_idx1 = replica_indices[r]
            rep_idx2 = replica_indices[r + 1]

            E1 = calc_energy(configs[rep_idx1], N, Jh, Jv, Jz)
            E2 = calc_energy(configs[rep_idx2], N, Jh, Jv, Jz)
            
            Δ = (β_replicas[r] - β_replicas[r + 1]) * (E1 - E2)
            p = min(1.0, exp(Δ))

            if rand() < p
                replica_indices[r], replica_indices[r + 1] = replica_indices[r + 1], replica_indices[r]
                temperature_indices[rep_idx1] = r + 1
                temperature_indices[rep_idx2] = r
            end
            
            exchange_probs_sum[r] += p
            exchange_attempts[r] += 1
        end
        exchange_count += 1


        for physical_replica_id in 1:M
            current_temp_idx = temperature_indices[physical_replica_id]
            current_state = replica_states[physical_replica_id]

            if current_state == 0 # WAITING
                if current_temp_idx == M
                    replica_states[physical_replica_id] = 1
                    round_trip_start_steps[physical_replica_id] = step
                end
            elseif current_state == 1 # GOING_UP
                if current_temp_idx == 1
                    replica_states[physical_replica_id] = 2
                end
            elseif current_state == 2 # GOING_DOWN
                if current_temp_idx == M
                    rtt = step - round_trip_start_steps[physical_replica_id]
                    push!(round_trip_times[physical_replica_id], rtt)

                    replica_states[physical_replica_id] = 1 # GOING_UP
                    round_trip_start_steps[physical_replica_id] = step
                end
            end
        end
    end

    all_rtts = vcat(round_trip_times...)

    exchange_prob_means = zeros(M - 1)
    for i in 1:(M - 1)
        if exchange_attempts[i] > 0
            exchange_prob_means[i] = exchange_probs_sum[i] / exchange_attempts[i]
        end
    end

    return exchange_prob_means, all_rtts
end

export initial_state, metropolis_sampler, calc_energy, score_method, calc_acceptance_and_rtt

end