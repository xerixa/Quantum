
using QuantumToolbox
using GLMakie

b = Bloch()
fig, _ = render(b)
point = [ 1/sqrt(3), 1/sqrt(3), 1 / sqrt(3)]
add_points!(b,point)
fig, _ = render(b)
fig

state1 = basis(2,0)
state2 = sigmax() * state1
state3 = normalize(basis(2,0) + basis(2,1))
add_states!(b, [state1, state2, state3])
fig, _ = render(b)
fig



#add_line!(b, state1, state 2)
#add_arc!(b, state1, state2)