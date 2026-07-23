#plot fiedelity vs gate time 

gate_times = 10.0:5.0:50.0 
fidelities_standard = Float64[]
fidelities_leakage = Float64[]


const N = 100 

for T in gate_times
    times = range(0.0, T, length=N)

    a = annihilate(3)
    g1 = 1/50.0
    g_phi = 1/70.0
    c_ops = [ sqrt(g1) * a, sqrt(2*g_phi) * (a'*a) ]
    
    sys = TransmonSystem(drive_bounds = U_BOUND, levels = LEVELS, δ = DELTA, dissipation = c_ops)
    op = EmbeddedOperator(:X, sys)
    
    initial_controls = 0.1 * randn(2, N)
    pulse0 = ZeroOrderPulse(initial_controls, times)
    qtraj = UnitaryTrajectory(sys, pulse0, op)
    

    qcp = SmoothPulseProblem(qtraj, N; ddu_bound = DDU_BOUND, Q = 100.0, R = 1e-2)
    Piccolo.solve!(qcp; max_iter = 50, eval_hessian=false)
    push!(fidelities_standard, fidelity(qcp))
    
    qcp_leakage = SmoothPulseProblem(
        qtraj, N; 
        ddu_bound = DDU_BOUND, Q = 100.0, R = 1e-2, 
        piccolo_options = PiccoloOptions(
            leakage_constraint = true, 
            leakage_constraint_value = 1e-2, 
            leakage_cost = 1e-2,
        ),
    )
    Piccolo.solve!(qcp_leakage; max_iter = 50, eval_hessian = false)
    push!(fidelities_leakage, fidelity(qcp_leakage))
end

plot(gate_times, fidelities_standard, label="Standard", marker=:circle, linewidth=2)
plot!(gate_times, fidelities_leakage, label="Leakage Constrained", marker=:square, linewidth=2)

xlabel!("Gate Time (ns)")
ylabel!("Infidelity" or "Fidelity") # Depending on if you want 1 - fidelity
title!("Piccolo.jl: Gate Time vs Fidelity Sweep")
