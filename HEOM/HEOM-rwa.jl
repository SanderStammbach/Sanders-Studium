import QuantumToolbox
using HierarchicalEOM
import Plots
using LaTeXStrings
using LinearSolve # to change the solver for better GPU performance
using CUDA
CUDA.allowscalar(false) # Avoid unexpected scalar indexing


W=2

#ϵ = -5
#U = 10
σm = sigmam() ## σ-
σz = sigmaz() ## σz
II = qeye(2)  ## identity matrix
λ = 1.0          # Reorganisationsenergie
γ = 0.1          # Zerfallsrate  
N =2     
kT=0.1
β = 1/kT ## Inverse Temperatur (1/kT)
# construct the annihilation operator for both spin-up and spin-down
# (utilize Jordan–Wigner transformation)
ω0 = 1
σz = sigmaz()
σx = sigmax() 
H0 = 0.5 * ω0 * σz
Td=0.2

# Define the operator that measures the 0, 1 element of density matrix
ρ01 = Qobj([0 1; 0 0])

ψ0 = (basis(2, 0) + basis(2, 1)) / √1;

P00 = ket2dm(basis(2, 0))
P11 = ket2dm(basis(2, 1))
op=sigmax()




using LinearAlgebra
using Plots

# Parameter
T = 0
λ = 0.1       # Badkopplungsstärke
γ = 0.02      # Charakteristische Frequenz des Bades
dt = 0.01     # Zeitschritt
time = 0:dt:50-dt  # Zeitachse
P_e = zeros(length(time))  # Wahrscheinlichkeit des angeregten Zustands
P_e[1] = 1.0  # Anfangswahrscheinlichkeit

# Korrelationsfunktion C(t)
function C(t)
    if t >= 0
        return -0.5 * γ * λ * (
            2 * cosh(t * γ) * (1 + cosh(T * γ)) +
            cosh((t - T) * γ) * (
                floor((π - 2 * angle(t - T)) / (4 * π)) -
                floor(0.75 + angle(t - T) / (2 * π))
            )
        )
    else
        return 0.0
    end
end

# Diskretisiere die Korrelationsfunktion
C_vals = [C(t) for t in time]

# Numerische Lösung der DGL
for k in 2:length(time)
    convolution = sum(C_vals[1:k] .* P_e[k:-1:1]) * dt
    P_e[k] = P_e[k-1] - dt * P_e[k-1] - dt * convolution
end

# Plot der Lösung
Plots.plot(time, abs2.(P_e), label="P_e(t)", xlabel="Zeit t", ylabel="Wahrscheinlichkeit P_e(t)", title="Lösung der DGL mit Faltung (Julia)", grid=true)



function drude_lorentz_rwa_params(λ, γ, β, N,W)
    # λ: Reorganisation (Kopplungsstärke)
    # γ: Zerfallsrate
    # β: Inverse Temperatur
    # N: Anzahl der Matsubara-Terme
    # W: Wheight 
    function matsubara(n)
        return 2 * π * n / β
    end
    
    # Liste der Matsubara-Frequenzen erstellen
    ϵ = [matsubara(n) for n in 1:N]

    η_absorb = ComplexF64[λ*cosh(W*Td/2)^2 * W * (cot(W * β / 2.0)-0.5 - 1.0im)]
    η_emit = ComplexF64[λ*cosh(W*Td/2)^2 * W * (cot(W * β / 2.0)+0.5 - 1.0im)]
    γ_absorb = ComplexF64[W]
    γ_emit = ComplexF64[W]

    if N > 0 #absorb hat 1/2*coth-0.5 und emit ist nb+1 und somit 1/2*cot +0.5
        for l in 1:N
            append!(η_absorb, 4 * λ*cosh(ϵ[l]*Td/2)^2 * W * ϵ[l] * (kT^2) / 2*(((ϵ[l] * kT)^2) - W^2)-0.5) 
            append!(η_emit, 4 * λ*cosh(ϵ[l]*Td/2)^2 * W * ϵ[l] * (kT^2) / 2*(((ϵ[l] * kT)^2) - W^2)+0.5)  
            append!(γ_absorb, ϵ[l] * kT)
            append!(γ_emit, ϵ[l] * kT)
        end
    end
    #δ = _boson_drude_lorentz_approx_discrepancy(λ, W, kT, N, η, γ)
    return η_absorb, γ_absorb, η_emit, γ_emit
    

end
η_absorb, γ_absorb, η_emit, γ_emit=drude_lorentz_rwa_params(λ, W, kT, N, W)
Bad_rwa=BosonBathRWA(σm, η_absorb, γ_absorb, η_emit, γ_emit)



tier = 2
#M = M_Boson(H0, tier, bath)
M_even_cpu = M_Boson(H0, tier,Bad_rwa )
M_even_gpu = cu(M_even_cpu)
tlist = 0:0.1:60
noPulseSol = HEOMsolve(M_even_gpu, ψ0, tlist; e_ops = [ρ01]);

Plots.plot(
    tlist,
    [ abs(noPulseSol.expect[1, :])],
    label = ["bla"],
    #linestyle = [:solid :dot :dash],
    linewidth = 3,
    xlabel = L"t",
    ylabel = L"rho_{01}",
    grid = false,
)
Plots.plot(
    time, 
    real.(P_e), 
    label="P_e(t)", 
    xlabel="Zeit t", 
    ylabel="Wahrscheinlichkeit P_e(t)", 
    title="Lösung der DGL mit Faltung (Julia)",
     grid=true)

display(Plots.plot(
    tlist,
    [real(noPulseSol.expect[1, :])],
    label=["e[t]"],
    #linestyle=[:solid],
    linewidth=3,
    xlabel=L"t",
    ylabel=L"\rho_{01}",
    grid=false
))





