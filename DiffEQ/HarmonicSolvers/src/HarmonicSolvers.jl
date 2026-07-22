module HarmonicSolvers

export step!, integrate,harmonic_oscillator!, energy, exact_undamped, SHO_check, state_matrix, input_vector, control, pendulum!, harmonic_oscillator, harmonic_oscillator_free!, run_freq_sweep, measure_amp, find_max_dt,AbstractODESolver, ForwardEuler, MyODEProblem, ODESolution, SimpleHarmonicOscillator, SimplePendulum, UnstableOscillator, SymplecticEuler, BackwardEuler, RungeKutta4

using LinearAlgebra
using DifferentialEquations
abstract type AbstractODESolver end
struct ForwardEuler <: AbstractODESolver end
struct RungeKutta4 <: AbstractODESolver end
struct BackwardEuler <: AbstractODESolver end
struct SymplecticEuler <: AbstractODESolver end

struct SimpleHarmonicOscillator
    m::Float64
    c::Float64
    k::Float64
    F::Float64
    ω::Float64
end

struct MyODEProblem{F,U,TS,P}
    f!::F
    u0::U
    tspan::TS
    p::P
end

struct ODESolution{T,U}
    t::T
    u::U
end

struct SimplePendulum
    m::Float64
    c::Float64
    L::Float64
    g::Float64
    F::Float64
    ω::Float64
end

struct UnstableOscillator
    m::Any
    c::Float64
    k::Float64
    F::Float64
    ω::Float64
end
#______________________________
# SYSTEMS
function SHO_check(m,c,k,F,ω)
   if m >0 && c>=0 && k>=0
        return SimpleHarmonicOscillator(m,c,k,F,ω)
    else
        println("Invalid variable definitions")
    end
end

function state_matrix(p::SimpleHarmonicOscillator)
    return A = [0.0 1.0; -p.k/p.m -p.c/p.m]
end

function input_vector(p::SimpleHarmonicOscillator)
    return B = [0.0; p.F/p.m]
end

function control(t, p::SimpleHarmonicOscillator)
    return p.F*cos(p.ω*t)
end

function harmonic_oscillator!(dx, x, p, t)
    A= state_matrix(p)
    B= input_vector(p)
    f = control(t,p)
    dx .= A * x + B * f
    return nothing
end


function exact_undamped(ts, x0, p)
    q0, v0 = x0
    ω0 = sqrt(p.k/p.m)
    q = q0 .* cos.(ω0 .* ts) .+ (v0 / ω0) .* sin.(ω0 .* ts)
    v = -q0 * ω0 .* sin.(ω0 .* ts) .+ v0 .* cos.(ω0 .* ts)
    return [q'; v']
end 

function energy(x, p::SimpleHarmonicOscillator)
    q,v = x
    ω0 = sqrt(p.k/p.m)
    return 0.5*(ω0^2 * q^2 + v^2)
end

function pendulum!(dx,x,p::SimplePendulum)
    θ = x[1]
    dθ = x[2]

    f = p.F * cos(p.ω*t)

    dx[1] = dθ
    dx[2] = -(p.c/p.m)* dθ - (p.g/p.L)*sin(θ)+(f/p.m)
    return nothing
end

function harmonic_oscillator(x,p,t)
    A = [0.0 1.0; -p.k/p.m -p.c/p.m]
    B = [0.0; p.F/p.m]
    f = p.F*cos(p.ω*t)
    return A*x + B*f
end

function harmonic_oscillator_free!(dx,x,p,t)
    pos = x[1]
    vel = x[2]

    f = p.F*cos(p.ω*t)

    dx[1] = vel
    dx[2] = -(p.k/p.m)*pos -(p.c/p.m) * vel + (f/p.m)
    return nothing
end


function theoretical_amp(ωd,p)
    return p.F/ sqrt((p.k-p.m*ωd^2)^2 + (p.c*ωd)^2)
end 

function measure_amp(sol_t, sol_x; delay=50.0, window_size=20.0)
    window_indices = findall(t-> (t>= delay) && (t <= delay + window_size), sol_t)

    if isempty(window_indices)
        error("Time window outside solution bounds")
    end

    window_x = sol_x[window_indices]

    return 0.5*(maximum(window_x) - minumum(window_x))
end

function run_freq_sweep(frequencies, p, solver)
    measured_amps = Float64[]
    theory_amps = Float64[]

    for ωd in frequencies
        p_current = (m=p.m, c=p.c, k=p.k, f =p.F, ω=ωd)
        t,x = solver(p_current, t_end=150.0)
        push!(measured_amps, measure_amp(t,x; delay=80.0, window_size=40.0))
        push!(theory_amps, theoretical_amp(ωd,p_current))
    end
    return theory_amps, measured_amps
end 

#________________________________
# SOVLERS
# FORWARD EULER

solver = ForwardEuler()

function step!(
    ::ForwardEuler, 
    f!, 
    xnext, 
    x,
    p,
    t,
    dt)
    dx = similar(x)
    f!(dx,x,p,t)
    xnext .= x .+ dx .* dt
    return xnext
end
#_____________________
# RUNGE KUTTA 4
function step!(
    ::RungeKutta4,
    f!,
    xnext,
    x,
    p,
    t,
    dt,
)

    k1 = similar(x)
    k2 = similar(x)
    k3 = similar(x)
    k4 = similar(x)
    temporary = similar(x) 
    f!(temporary, x, p, t) 
    k1 .= temporary
    f!(temporary, x .+ 0.5 .* dt .* k1, p, t + 0.5 * dt)
    k2 .= temporary
    f!(temporary, x .+ 0.5 .* dt .* k2, p, t + 0.5 * dt)
    k3 .= temporary
    f!(temporary, x .+ dt .* k3, p, t + dt)
    k4 .= temporary
    xnext .= x .+ (dt/6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    return xnext
end
#_____________________
# BACKWARD EULER
function step!(
    ::BackwardEuler,
    f!,
    xnext,
    x,
    p::SimpleHarmonicOscillator,
    t,
    dt
)
    A = state_matrix(p)
    B = input_vector(p)
    matrix = I - dt * A
    rhs = x .+ dt .* (B .* control(t + dt, p))
    xnext .= matrix \ rhs
    return xnext
end
#_____________________
# SYMPLECTIC EULER
function step!(
    ::SymplecticEuler,
    f!,
    xnext,
    x,
    p::SimpleHarmonicOscillator,
    t,
    dt
)
    q,v = x
    dx = similar(x)

    f!(dx,x,p,t)
    a = dx[2]
 
    v_next = v + dt * a
    q_next = q + dt * v_next

    xnext[1] = q_next
    xnext[2] = v_next
    return xnext
end
#_____________________

function integrate(problem::MyODEProblem, solver::AbstractODESolver;dt)
    t0, tf = problem.tspan
    @assert tf > t0 "Final time must be greater than start time"
    @assert dt > 0 "Time step must be positive"
    nsteps = Int(ceil((tf-t0)/dt)) #if the div isnt perfect you still a small step to fill the gap
    times = range(t0, tf, length=nsteps + 1) 
    states = [
        similar(problem.u0)
        for _ in eachindex(times)
    ]
    states[1] .= problem.u0
    # TODO: initialize states[1]
    for n in 1:nsteps
        # TODO: call step!
        step!(solver, problem.f!, states[n+1], states[n], problem.p, times[n], dt)
        #next state = current + (roc *dt)
        #next state = current +(ODE*dt)
        #solver , f, u, t, dt, p(physical constants of the system)
    end
    return ODESolution(times, states)
end


end




