module XYCore

using Random, LinearAlgebra, Statistics, StatsBase, Printf, Distributions

# seed = 42
# Random.seed!(seed)

function initial_state(N)
    state = 2π .* rand(N, N)
    return state
end

function metropolis_sampler(config, β, N)
    range_i = shuffle!(collect(1:N))
    range_j = shuffle!(collect(1:N))
    for i in range_i
        for j in range_j
            s = config[i, j]
            nb_cos_sum = cos(s - config[mod1(i + 1, N), j]) + cos(s - config[i, mod1(j + 1, N)]) +
                         cos(s - config[mod1(i - 1, N), j]) + cos(s - config[i, mod1(j - 1, N)])
            s_new = 2π * rand()
            
            nb_cos_sum_new = cos(s_new - config[mod1(i + 1, N), j]) + cos(s_new - config[i, mod1(j + 1, N)]) +
                             cos(s_new - config[mod1(i - 1, N), j]) + cos(s_new - config[i, mod1(j - 1, N)])
            cost = -(nb_cos_sum_new - nb_cos_sum)
            acceptance_prob = min(1.0, exp(-cost * β))
            if rand() < acceptance_prob
                config[i, j] = s_new
            end
        end
    end
    
    return config
end

function wolff_sampler(config, β, N)
    r_angle = 2π * rand()
    r_vec = [cos(r_angle), sin(r_angle)]
    seed_i, seed_j = rand(1:N), rand(1:N)
    cluster_stack = [(seed_i, seed_j)]
    in_cluster = falses(N, N)
    in_cluster[seed_i, seed_j] = true
    
    head = 1
    while head <= length(cluster_stack)
        (curr_i, curr_j) = cluster_stack[head]
        head += 1

        s_curr_angle = config[curr_i, curr_j]
        s_curr_vec = [cos(s_curr_angle), sin(s_curr_angle)]

        neighbors = [
            (mod1(curr_i - 1, N), curr_j), (mod1(curr_i + 1, N), curr_j),
            (curr_i, mod1(curr_j - 1, N)), (curr_i, mod1(curr_j + 1, N))
        ]

        for (next_i, next_j) in neighbors
            if !in_cluster[next_i, next_j]
                s_next_angle = config[next_i, next_j]
                s_next_vec = [cos(s_next_angle), sin(s_next_angle)]

                p_add = 1 - exp(min(0, 2 * β * dot(s_curr_vec, r_vec) * dot(s_next_vec, r_vec)))

                if rand() < p_add
                    in_cluster[next_i, next_j] = true
                    push!(cluster_stack, (next_i, next_j))
                end
            end
        end
    end

    for (i, j) in cluster_stack
        s_angle = config[i, j]
        s_vec = [cos(s_angle), sin(s_angle)]

        # s_new = s - 2 * (s・r) * r
        s_reflected_vec = s_vec - 2 * dot(s_vec, r_vec) * r_vec

        config[i, j] = atan(s_reflected_vec[2], s_reflected_vec[1])
    end
    
    return config
end

function calc_energy(config, N)
    energy = 0.0

    for i in 1:N
        for j in 1:N
            S = config[i, j]
            nb = cos(S - config[mod1(i + 1, N), j]) + cos(S - config[i, mod1(j + 1, N)])
            energy -= nb
        end
    end
    return energy
end

function calc_magnetization(config, N)
    mx = 0.0
    my = 0.0
    for i in 1:N, j in 1:N
        mx += cos(config[i, j])
        my += sin(config[i, j])
    end
    return sqrt(mx^2 + my^2) / (N * N)
end

function calc_helicity_terms(config, N)
    E_x = 0.0
    J_x = 0.0
    for i in 1:N, j in 1:N
        delta_theta = config[i, j] - config[mod1(i + 1, N), j]
        E_x += cos(delta_theta)
        J_x += sin(delta_theta)
    end
    return E_x, J_x
end

function calc_correlation_function(config, N, max_r)
    correlations = zeros(max_r)
    counts = zeros(Int, max_r)
    center_i, center_j = div(N, 2), div(N, 2)
    s_ref = config[center_i, center_j]

    for i in 1:N, j in 1:N
        dx = abs(i - center_i)
        dy = abs(j - center_j)
        dx = min(dx, N - dx)
        dy = min(dy, N - dy)
        r = round(Int, sqrt(dx^2 + dy^2))

        if 1 <= r <= max_r
            correlations[r] += cos(s_ref - config[i, j])
            counts[r] += 1
        end
    end

    for r in 1:max_r
        if counts[r] > 0
            correlations[r] /= counts[r]
        end
    end
    return correlations
end

function score_method(
    N,
    β_replicas; 
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
            configs[r] = metropolis_sampler(configs[r], β_replicas[r], N)
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
            configs[r] = metropolis_sampler(configs[r], β_replicas[r], N)
            energy_sum[r] += calc_energy(configs[r], N)
        end

        if step % exchange_interval == 0
            start = isodd(step_count) ? 2 : 1

            for r in start:2:(M - 1)
                E1 = calc_energy(configs[r], N)
                E2 = calc_energy(configs[r + 1], N)
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
    config,
    β_replicas;
    mcSteps=10^4
)
    M = length(β_replicas)
    configs = [copy(config) for _ in 1:M]
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
            configs[current_replica_idx] = metropolis_sampler(configs[current_replica_idx], β_replicas[i], N)
        end

        start_idx = isodd(exchange_count) ? 2 : 1
        for r in start_idx:2:(M - 1)
            rep_idx1 = replica_indices[r]
            rep_idx2 = replica_indices[r + 1]

            E1 = calc_energy(configs[rep_idx1], N)
            E2 = calc_energy(configs[rep_idx2], N)
            
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

export initial_state, metropolis_sampler, wolff_sampler, calc_energy, calc_magnetization, score_method, calc_acceptance_and_rtt

end
