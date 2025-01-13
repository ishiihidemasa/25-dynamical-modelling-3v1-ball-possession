# Generate dataset with randomly sampled sets of parameter values
# NOTE: This script was run on HPC
using Distributed
using CSV
using Dates
using TimeZones

# print start time
starttime = Dates.now(tz"Europe/Berlin")
println("start: ", starttime)
# Change directory if needed
wd = "scripts"
if !occursin(wd, pwd())
    cd(wd)
end
println("I am at ", pwd())

##### simulation settings #####
seed_offset = 0  # Seed starts from seed_offset
numsample = 10  # Number of sets of parameter values
numtrial = 500  # Number of trials for each set of param values
# the number of workers
# CAUTION: num_workers + 1 threads will be occupied by this script!
num_workers = 10  # main process is not included

# fixed parameters
Δt = 1e-2
maxpass = 20
L = 6.0
L_out = 0.5
L_buff = 0.5
# value ranges for other parameters
d_range = Dict(
    :m => (55, 75),
    :T => (0.5, 1.2),
    :σ => (3, 9),
    :β => (0, 25),
    :kr => (0, 1000),
    :γ => (0, 1000),
    :kf => (0, 200),
    :Lf => (3, 8),
    :ke => (0, 200),
    :Le => (1, 6),
    :q => (0.2, 0.8),
    :τ => (0.5, 1.5),
    :L_ict => (0.4, 0.8)
)

# Path for input / output files
dir_name = "sim-dataset"
out_file = "$dir_name-numpass-data.csv"
###############################

# add worker processes
# - reduce the number of workers if there are not enough threads
addprocs(min(num_workers, max(1, Sys.CPU_THREADS - 4)))
# include the custom module in all processes
Distributed.remotecall_eval(
    Main, procs(), :(include("../src/BallPossessionModel.jl"))
)

# define a task for all workers
@everywhere begin
    # packages must be loaded on every process
    import .BallPossessionModel as BPM

    using DataFrames
    using Random
    import Statistics: mean, var

    """ 
    Sample a value within d_range[key], record it in d_output, 
    and return the value.
    """
    function record_and_return!(d_output, rng, d_range, key)
        val = (
            d_range[key][1] + 
            (d_range[key][2] - d_range[key][1]) * rand(rng)
        )
        d_output[key] = val
        return val
    end

    function run_seed(
        seed, d_range, numtrial, Δt, L, maxpass, L_out, L_buff, dir_name, v_keys
    )
        try
            println("working on seed $seed...")
            rng = Xoshiro(seed)

            # Generate and record parameter values
            d_output = Dict(
                :seed => seed,
                :numtrial => numtrial,
                :Δt => Δt,
                :L => L,
                :maxpass => maxpass,
                :L_out => L_out,
                :L_buff => L_buff,
            )

            p = BPM.Parameters(
                Δt=Δt,
                maxpass=maxpass,
                σ=record_and_return!(d_output, rng, d_range, :σ),
                β=record_and_return!(d_output, rng, d_range, :β),
                q=record_and_return!(d_output, rng, d_range, :q),
                L_ict=record_and_return!(d_output, rng, d_range, :L_ict),
                L_out=L_out,
                L_buff=L_buff,
                m=record_and_return!(d_output, rng, d_range, :m),
                L=L,
                T=record_and_return!(d_output, rng, d_range, :T),
                kr=record_and_return!(d_output, rng, d_range, :kr),
                γ=record_and_return!(d_output, rng, d_range, :γ),
                kf=record_and_return!(d_output, rng, d_range, :kf),
                Lf=record_and_return!(d_output, rng, d_range, :Lf),
                ke=record_and_return!(d_output, rng, d_range, :ke),
                Le=record_and_return!(d_output, rng, d_range, :Le),
                τ=record_and_return!(d_output, rng, d_range, :τ),
            )

            # Arrays to store results
            a_passcount = Vector{Int64}(undef, numtrial)
            a_passwide = Vector{Int64}(undef, numtrial)
            a_avearea = Vector{Float64}(undef, numtrial)
            a_vararea = Vector{Float64}(undef, numtrial)

            # Record the number of appearances of each message directly into d_output
            d_output[:intercept] = 0
            d_output[:max_pass] = 0
            d_output[:ball_out] = 0

            # Simulate
            for i in 1:numtrial
                res = BPM.simulate(rng, p)
                # Record results
                d_output[res.message] += 1
                a_passcount[i] = res.passcount
                a_passwide[i] = sum(res.pass2wide[1:end-1])
                area_ts = BPM.get_area_timeseries(res)[1:res.endid]
                a_avearea[i] = mean(area_ts)
                a_vararea[i] = var(area_ts)
            end

            d_output[:passcount_mean] = mean(a_passcount)
            d_output[:passcount_var] = var(a_passcount)
            d_output[:passwide_freq] = sum(a_passwide) / sum(a_passcount)
            d_output[:avearea_mean] = mean(a_avearea)
            d_output[:vararea_mean] = mean(a_vararea)

            @debug println("seed $seed: done")
            return DataFrame([[d_output[key]] for key in v_keys], v_keys)

        catch e
            # When an exception was raised, write it out
            logfile = "$dir_name/so$seed-error.txt"
            error_msg = sprint(showerror, e)
            st = sprint(
                (io,v) -> show(io, "text/plain", v), 
                stacktrace(catch_backtrace())
            )
            open(logfile, "w") do f
                write(f, "$(error_msg)\n$(st)")
            end
            println("seed $seed: exception was raised")
            return
        end
    end
end

# setup a channel for communication
const results_channel = RemoteChannel(() -> Channel{Any}(64))

# function to write results to CSV
function writer_task(filename)
    println("writer_task started.")
    while true
        result = take!(results_channel)
        if result === :done
            break
        end
        CSV.write(filename, result, append=true)
    end
    return "writer_task terminates."
end

# Load header for the output CSV file
header_file = "$dir_name/$dir_name-header.csv"
t_header = CSV.read(header_file, DataFrame) |> names

if !isfile("$dir_name/$out_file")
    # Create the output file if it does not exist
    CSV.write(
        "$dir_name/$out_file", 
        DataFrame([[] for _ in 1:length(t_header)], t_header)
    )
end

# start the writer process
wfunc() = writer_task("$dir_name/$out_file")
wtask = Task(wfunc)
schedule(wtask)  # start the writer

# distribute tasks to workers: wait until the loop finishes
pmap(seed_offset + 1:seed_offset + numsample; on_error=identity) do seed
    df_results = run_seed(
        seed, Dict(d_range), numtrial, Δt, L, maxpass, L_out, L_buff, 
        dir_name, [Symbol(key) for key in t_header]
    )
    put!(results_channel, df_results)
    println("seed $seed: sent results!")
end

# tell the writer_task to finish
put!(results_channel, :done)
# wait until writer_task() terminates
println(fetch(wtask))

# print end time and duration
endtime = Dates.now(tz"Europe/Berlin")
println("end: ", endtime)
elapsedsec = Dates.value(endtime - starttime) / 1000
println(
    "elapsed: $(Int(elapsedsec ÷ 3600)) h $(Int(elapsedsec % (3600) ÷ 60)) min $(Int(round(elapsedsec % 360) % 60)) sec"
)

