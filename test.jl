using CollectiveSpins

# Define geometry of system
N = 5     # Number of spins
a = 0.3   # spin-spin distance
geometry = CollectiveSpins.geometry.chain(a, N)

# Create system consisting of N spins in the defined geometry
e = [0,0,1]   # Quantization axis
system = CollectiveSpins.SpinCollection(geometry, e)

# Initial quantum state
phi = 0.
theta = pi/2
Ψ0 = CollectiveSpins.mpc.blochstate(phi, theta, N)

# Perform time evolution according to master equation
T = [0:0.05:5.;]
tout, ρt = CollectiveSpins.mpc.timeevolution(T, system, Ψ0)
tout2, ρt2 = CollectiveSpins.mpc.timeevolution(T, system, ρt[end])