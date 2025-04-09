@kwdef mutable struct Parameters
    const Δt::Float64
    const maxpass::Int64
    const σ::Float64
    const β::Float64
    const q::Float64
    const L_ict::Float64
    const L_out::Float64
    const L_buff::Float64
    const m::Float64
    const L::Float64
    const T::Float64
    const kr::Float64  # spring constant for returning force
    const γ::Float64  # coeff. of viscocity outside the field
    const kf::Float64
    const Lf::Float64
    const ke::Float64
    const Le::Float64
    const τ::Float64
    xPtp::MVector{2, Float64} = MVector{2, Float64}(undef)
    xRtp::MVector{2, Float64} = MVector{2, Float64}(undef)
    xRgoal::MVector{2, Float64} = MVector{2, Float64}(undef)
    if_chase::Bool = true
end

@kwdef struct Result
    t::Vector{Float64}
    u::Matrix{Float64}
    message::Symbol
    positionbyid::Matrix{Float64}
    passcount::Int64
    endtime::Float64
    endid::Int64
    passtime::Vector{Float64}
    pass2wide::Vector{Bool}
    passangle::Vector{Float64}
end

# social forces for Passer and Mover
"""
    returnforce(x, kr, L)

Passer and Mover cannot move away from the field outside.
"""
function returnforce(x, kr, L, L_buff)
    return @. ifelse(
        abs(x) < (L - L_buff) / 2,
        0,
        kr * (sign(x) * (L - L_buff) - 2x) / L_buff
    )
end

""" Nonlinear damping around the border """
function damping(x, v, γ, L, L_buff)
    return @. ifelse(
        L - L_buff <= 2 * abs(x) <= L + L_buff,
        -γ * v^3,
        0
    )
end

"""
    evadeforce(x, xD, ke, Le)

Passer and Mover are connected to DF by a spring whose 
natural length and coefficient are ke and Le, respectively.
"""
function evadeforce(x, xD, ke, Le)
    distance = norm(x - xD)
    return ke * (Le - distance) * (x - xD) / distance
end

"""
    followforce(x, xR, kf, Lf)

Passer and Mover are connected to Receiver by a spring
whose natural length and coefficient are kf and Lf, respectively.
"""
function followforce(x, xR, kf, Lf)
    distance = norm(xR - x)
    return kf * (distance - Lf) * (xR - x) / distance
end

"""
    getincenter(x1, x2, x3)

Return the position vector of the incenter of the triangle
given by `x1`, `x2`, and `x3`.
"""
function getincenter(x1, x2, x3)
    l12 = norm(x1 - x2)
    l13 = norm(x1 - x3)
    l23 = norm(x2 - x3)
    return @. (l23 * x1 + l13 * x2 + l12 * x3) / (l12 + l13 + l23)
end

"""
    equationofmotion!(du, u, p, t)

Equations of motion during a pass.

# Arguments
- `du`
- `u`: vector of state variables, which is the transpose of
    `[vP vM xP xM xR xD xB]` where each element is 
    a two dimensional vector.
- `p::Parameters`: Parameter values.
- `t`: time
"""
function equationofmotion!(du, u, p::Parameters, t)
    # Ball
    du[13:14] = (p.xRgoal - p.xPtp) / p.T
    # Receiver
    du[9:10] = (p.xRgoal - p.xRtp) / p.T
    # DF
    if p.if_chase
        # DF runs towards the pass course
        du[11:12] = (
            (u[9:10] + u[13:14]) / 2 - u[11:12]  # (xR + xB) / 2 - xD
        ) / p.τ
    else
        # DF runs towards the incenter
        du[11:12] = (
            getincenter(u[5:6], u[7:8], u[9:10]) - u[11:12]
        ) / p.τ
    end
    # Passer and Mover
    # - positions
    du[5:6] = u[1:2]  # dxP = vP
    du[7:8] = u[3:4]  # dxM = vM
    # - velocities
    du[1:2] = (
        returnforce(u[5:6], p.kr, p.L, p.L_buff) +
        damping(u[5:6], u[1:2], p.γ, p.L, p.L_buff) +
        followforce(u[5:6], u[9:10], p.kf, p.Lf) +
        evadeforce(u[5:6], u[11:12], p.ke, p.Le)
    ) / p.m  # vP
    du[3:4] = (
        returnforce(u[7:8], p.kr, p.L, p.L_buff) +
        damping(u[7:8], u[3:4], p.γ, p.L, p.L_buff) +
        followforce(u[7:8], u[9:10], p.kf, p.Lf) +
        evadeforce(u[7:8], u[11:12], p.ke, p.Le)
    ) / p.m  # vM
    return nothing
end

function get_xRgoal(xRtp, xPtp, θ_p, L)
    # Ignore the field constraint at first
    rotationmatrix(θ) = @MMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]
    return(xPtp + rotationmatrix(θ_p) * (xRtp - xPtp))
end

function truncate_sol(sol, id, extra_samples=20)
    last = minimum([length(sol.t), id + extra_samples])
    return (
        sol.t[1:last],
        Array(sol[:, 1:last])
    )
end

function record_int(sol, icttime, id_ict)
    return(
        :intercept,
        id_ict,
        icttime,
        truncate_sol(sol, id_ict)...,
    )
end

function record_bout(sol, bouttime, id_bout)
    return(
        :ball_out,
        id_bout,
        bouttime,
        truncate_sol(sol, id_bout)...,
    )
end

"""
    difference_in_angle(ang1, ang2)

Calculate `ang2 - ang1` which is in `[-pi, pi]`.
"""
function difference_in_angle(ang1, ang2)
    ϕ = ang2 - ang1
    return abs(ϕ) > π ? ϕ - 2π * sign(ϕ) : ϕ
end

"""
    simulate(rng, p)

Run a simulation of the ball possession model.
"""
function simulate(rng::Random.AbstractRNG, p::Parameters)
    # initialize
    t = 0
    passcount = 0
    message = :continue
    if_chase = false  # DF remains centered at the first pass
    #if_chase = rand(rng, Bool)  # random initial strategy for DF

    # generate initial condition
    # u is the transpose of [vP vM xP xM xR xD xB]
    u_old = MVector{14, Float64}(undef)
    # initial velocities
    u_old[1:4] .= 0  # no initial velocities
    # initial position of DF
    u_old[11:12] .= 0  # xD = [0, 0]
    # initial position of OFs: equiangular formation
    # - Passer
    xPtp = MVector{2, Float64}(p.L / 2 * [1, 0])
    u_old[5:6] = xPtp
    # - Receiver
    xRtp = MVector{2, Float64}(p.L / 2 * [cos(-2π / 3), sin(-2π / 3)])
    u_old[9:10] = xRtp
    # - Mover
    u_old[7:8] = MVector{2, Float64}(p.L / 2 * [cos(2π / 3), sin(2π / 3)])
    # Ball: Passer has the ball
    u_old[13:14] = xPtp

    # prepare to store results
    res_t = [t]
    res_u = deepcopy(u_old)

    endtime = 0.0
    endid = 1  # current length of res_t
    passtimelist = Float64[]  # empty vector
    pass2widelist = [true]
    passanglelist = Float64[]  # empty vector

    # position of OFs by ID, not role
    role2id = Dict(:P => :OF1, :R => :OF2, :M => :OF3)
    positionbyid = Dict(
        :OF1 => Matrix{Float64}(undef, 2, 1),
        :OF2 => Matrix{Float64}(undef, 2, 1),
        :OF3 => Matrix{Float64}(undef, 2, 1)
    )
    positionbyid[:OF1][:, 1] = xPtp
    positionbyid[:OF2][:, 1] = xRtp
    positionbyid[:OF3][:, 1] = u_old[7:8]

    # run simulation
    while message == :continue
        @debug println("Pass is made at t = ", round(t, digits=2))
        push!(passtimelist, t)
        n_step::Int64 = round(p.T / p.Δt)
        rounded_T = n_step * p.Δt

        θ_p = randn(rng) * p.σ * π / 180
        push!(passanglelist, θ_p)
        xRgoal = get_xRgoal(xRtp, xPtp, θ_p, p.L)

        # DF probabilistically chases the ball.
        # `if_chase`; `if_continue` -> `if_chase` in the next pass
        # true; true -> true
        # true; false -> false
        # false; true -> false
        # false; false -> true
        if_continue = rand(rng) < p.q
        if_chase = if_chase == if_continue

        # solve ODEs for one pass
        p.xPtp = xPtp
        p.xRtp = xRtp
        p.xRgoal = xRgoal
        p.if_chase = if_chase
        sol = solve(
            ODEProblem(equationofmotion!, u_old, (t, t + rounded_T), p), 
            saveat=p.Δt, 
            save_start=false,
            abstol=1e-8, 
            reltol=1e-8
        )

        # termination conditions
        # check intercept
        if_ict_array = (
            # norm along dimension 1
            sqrt.(dropdims(sum(
                abs2, 
                sol[11:12, :] - sol[13:14, :],
                dims=1
            ), dims=1))
            .< p.L_ict  # is less than L_ict
        )
        if_ict = any(if_ict_array)
        if if_ict
            # pass was ict by DF
            # extract the time of pass ict
            id_ict = collect(1:length(sol.t))[if_ict_array][1]
            icttime = sol.t[id_ict]
            @debug "intercept at time $(round(icttime, digits=3))!"
        end
        # check ball_out
        if_bout_array = (
            # distance between the edge and the ball
            dropdims(maximum(abs.(sol[13:14, :]), dims=1), dims=1) .- p.L / 2
            .> p.L_out  # is greater than L_out
        )
        if_bout = any(if_bout_array)
        if if_bout
            # ball went out of the field
            # extract the time of ball_out
            id_bout = collect(1:size(sol.t, 1))[if_bout_array][1]
            bouttime = sol.t[id_bout]
            @debug "ball_out at time $(round(bouttime, digits=3))!"
        end
        
        # truncate solution if necessary
        message, endid_n, endtime, arr_t, arr_u = begin
            if if_bout && if_ict
                if icttime <= bouttime
                    # intercept took place first
                    record_int(sol, icttime, id_ict)
                else
                    # ball_out took place first
                    record_bout(sol, bouttime, id_bout)
                end
            elseif if_ict
                # only intercept took place
                record_int(sol, icttime, id_ict)
            elseif if_bout
                # only out_of_field took place
                record_bout(sol, bouttime, id_bout)
            else
                (:continue, length(sol.t), sol.t[end], sol.t, Array(sol))
                # Array(sol)'s shape is 14 x (the num. of samples)
            end
        end
        endid += endid_n  # update endid

        # record state variables by role
        res_t = vcat(res_t, arr_t)
        res_u = hcat(res_u, arr_u)

        # record by player id
        # - current Passer
        positionbyid[role2id[:P]] = hcat(
            positionbyid[role2id[:P]], arr_u[5:6, :]
        )
        # - current Receiver
        positionbyid[role2id[:R]] = hcat(
            positionbyid[role2id[:R]], arr_u[9:10, :]
        )
        # - current Mover
        positionbyid[role2id[:M]] = hcat(
            positionbyid[role2id[:M]], arr_u[7:8, :]
        )

        # update pass-count
        if message == :continue
            # no intercept & ball remained within the field
            passcount += 1
        end
        # check if pass-count reached the limit
        if passcount >= p.maxpass
            message = :max_pass
        end

        # terminate when neccesary
        if message != :continue
            break
        end

        # Simulation continues! Prepare for the next epoch
        t += rounded_T

        # switch OFs' roles
        role2id_old = Dict(role2id)
        u_last = arr_u[:, end]  # make a copy because it loops along dim 1
        # determine if Mover becomes Receiver
        # - calculate pass courses (angles) from R to M & P
        get_angle(vec) = angle(complex(vec[1], vec[2]))
        ang_RD = get_angle(u_last[11:12] - u_last[9:10])
        # take difference
        ang_MRD = abs(difference_in_angle(
            get_angle(u_last[7:8] - u_last[9:10]),  # ang_RM
            ang_RD
        ))
        ang_PRD = abs(difference_in_angle(
            get_angle(u_last[5:6] - u_last[9:10]),  # ang_RP
            ang_RD
        ))
        # determine the pass direction
        sameM = (
            rand(rng) < (
                # probability that Receiver makes a pass to Passer
                exp(p.β * ang_PRD) /
                (exp(p.β * ang_PRD) + exp(p.β * ang_MRD))
            )
        )

        # record if the pass goes to the player with a wider pass course
        if ang_PRD == ang_MRD
            push!(pass2widelist, true)
        elseif ang_PRD > ang_MRD
            push!(pass2widelist, sameM)
        else
            push!(pass2widelist, !sameM)
        end

        # update u_old: switch roles among OFs
        # velocities of P & M: reset to zero
        u_old[1:4] .= 0
        # Receiver always becomes Passer
        role2id[:P] = role2id_old[:R]
        xPtp[:] = deepcopy(u_last[9:10])
        u_old[5:6] = xPtp
        # next Mover & Receiver
        if sameM
            # Mover remains Mover
            u_old[7:8] = u_last[7:8]
            # Passer becomes Receiver
            role2id[:R] = role2id_old[:P]
            xRtp[:] = u_last[5:6]
            u_old[9:10] = xRtp
        else
            # Mover becomes Receiver
            role2id[:R] = role2id_old[:M]
            xRtp[:] = u_last[7:8]
            u_old[9:10] = xRtp
            # Passer becomes Mover
            role2id[:M] = role2id_old[:P]
            u_old[7:8] = u_last[5:6]
        end
        # no change regarding the Ball & DF
        u_old[11:14] = u_last[11:14]
    end

    # return resuts as a BPMResult object
    return Result(
        t=res_t,
        u=res_u,
        message=message,
        positionbyid=vcat(
            positionbyid[:OF1],
            positionbyid[:OF2],
            positionbyid[:OF3]
        ),
        passcount=passcount,
        endtime=endtime,
        endid=endid,
        passtime=passtimelist,
        pass2wide=pass2widelist,
        passangle=passanglelist,
    )
end

function param2strlist(p::Parameters, seed::Int)
    d_label = Dict(
        :kf => L"k_\mathrm{f}",
        :Lf => L"L_\mathrm{f}",
        :ke => L"k_\mathrm{e}",
        :Le => L"L_\mathrm{e}",
        :β => L"\beta",
        :σ => L"\sigma",
        :T => L"T",
        :q => L"q",
        :τ => L"\tau",
        :L_ict => L"L_\mathrm{ict}",
        :m => L"m",
        :L => L"L",
        :L_out => L"L_\mathrm{out}",
        :kr => L"k_\mathrm{r}",
        :γ => L"\gamma",
        :maxpass => L"n_\mathrm{max}",
        :Δt => L"\Delta t",
    )
    l_labvar = String[]
    for (key, lab) in d_label
        val = getproperty(p, key)
        valstr = typeof(val) <: Int ? "$(val)" : "$(round(val, digits=2))"
        push!(
            l_labvar,
            lab[1:end-1] * " = $valstr \$"
        )
    end
    push!(l_labvar, L"\mathrm{seed} = %$(seed)")
    return l_labvar
end
