using cascadeddecay
using PyCall
@pyimport matplotlib.pyplot as plt

const A = [0.01:0.01:5]
const γ = 1.
#const Θ = pi/2

plt.figure(1, figsize=(6,4))
plt.figure(2, figsize=(6,4))

for Θ=range(0, stop=pi/2, length=11)
    Ω = [Omega(a, Θ, γ) for a=A]
    Γ = [Gamma(a, Θ, γ) for a=A]

    plt.figure(1)
    plt.plot(A, Ω, "b")
    plt.ylim(-0.5,0.5)

    plt.figure(2)
    plt.plot(A, Γ, "g")
    plt.ylim(-0.5,0.5)
end

plt.figure(1)
plt.ylim(-0.5,0.5)
plt.xlabel("a")
plt.ylabel("\$\\Omega(a,\\Theta)\$")
plt.savefig("images/collective_interaction_omega.pdf")

plt.figure(2)
plt.ylim(-0.5,0.5)
plt.xlabel("a")
plt.ylabel("\$\\Gamma(a,\\Theta)\$")
plt.savefig("images/collective_interaction_gamma.pdf")

plt.show()