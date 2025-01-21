import QuantumToolbox
using HierarchicalEOM
import Plots
using LaTeXStrings



W=2

ϵ = -5
U = 10
σm = sigmam() ## σ-
σz = sigmaz() ## σz
II = qeye(2)  ## identity matrix
λ = 1.0          # Reorganisationsenergie
γ = 0.1          # Zerfallsrate
β = 1.0          # Inverse Temperatur (1/kT)
N = 5     
kT=1/β
# construct the annihilation operator for both spin-up and spin-down
# (utilize Jordan–Wigner transformation)
ω0 = 1
σz = sigmaz()
σx = sigmax()
H0 = 0.5 * ω0 * σz

# Define the operator that measures the 0, 1 element of density matrix
ρ01 = Qobj([0 1; 0 0])

ψ0 = (basis(2, 0) + basis(2, 1)) / √2;

P00 = ket2dm(basis(2, 0))
P11 = ket2dm(basis(2, 1))
op=sigmax()

function Boson_DrudeLorentz_Matsubara_self(op, λ, W, kT, N)
    β = 1.0 / kT

    function matsubara(n)
        return 2 * π * n / β
    end
    
    # Liste der Matsubara-Frequenzen erstellen
    ϵ = [matsubara(n) for n in 1:N]

    η = ComplexF64[λ * W * (cot(W * β / 2.0) - 1.0im)]
    γ = ComplexF64[W]

    if N > 0
        for l in 1:N
            append!(η, 4 * λ * W * ϵ[l] * (kT^2) / (((ϵ[l] * kT)^2) - W^2))
            append!(γ, ϵ[l] * kT)
        end
    end

    #δ = _boson_drude_lorentz_approx_discrepancy(λ, W, kT, N, η, γ)

    return BosonBath(op, η, γ)
end


bath=Boson_DrudeLorentz_Matsubara_self(op, λ, ω0, kT, N)
tier = 6
M = M_Boson(H0, tier, bath)
tlist = 0:0.5:5
noPulseSol = HEOMsolve(M, ψ0, tlist; e_ops = [ρ01]);

Plots.plot(
    tlist,
    [ real(noPulseSol.expect[1, :])],
    label = ["bla"],
    #linestyle = [:solid :dot :dash],
    linewidth = 3,
    xlabel = L"t",
    ylabel = L"rho_{01}",
    grid = false,
)
display(Plots.plot(
    tlist,
    [real(noPulseSol.expect[1, :])],
    label=["bla"],
    #linestyle=[:solid],
    linewidth=3,
    xlabel=L"t",
    ylabel=L"\rho_{01}",
    grid=false
))

Plots.plot(
    tlist,
    [real(noPulseSol.expect[1, :]), sin.(tlist)],  # Beispiel für zwei Kurven
    label=["Kurve 1" "Kurve 2"],
    linestyle=[:solid :dash],  # Jeweils ein Stil pro Kurve
    linewidth=3,
    xlabel=L"t",
    ylabel=L"\rho_{01}",
    grid=false
)

function drude_lorentz_rwa_params(λ, γ, β, N)
    # λ: Reorganisation (Kopplungsstärke)
    # γ: Zerfallsrate
    # β: Inverse Temperatur
    # N: Anzahl der Matsubara-Terme

    # Hauptterm
    η_absorb = [λ * γ / β * coth(β * γ / 2)]
    γ_absorb = [γ]
    η_emit = [λ * γ / β]
    γ_emit = [γ]

    # Matsubara-Terme
    for n in 1:N
        ω_n = 2 * π * n / β  # Matsubara-Frequenz
        η_n = 4 * λ * γ * ω_n / (β * (ω_n^2 + γ^2))
        push!(η_absorb, η_n)
        push!(γ_absorb, ω_n)
        push!(η_emit, η_n)
        push!(γ_emit, ω_n)
    end

    return η_absorb, γ_absorb, η_emit, γ_emit
end

# Beispielparameter
λ = 1.0          # Reorganisationsenergie
γ = 0.1          # Zerfallsrate
β = 1.0          # Inverse Temperatur (1/kT)
N = 5            # Anzahl der Matsubara-Terme

η_absorb, γ_absorb, η_emit, γ_emit = drude_lorentz_rwa_params(λ, γ, β, N)

println("η_absorb: $η_absorb")
println("γ_absorb: $γ_absorb")
println("η_emit: $η_emit")
println("γ_emit: $γ_emit")




#bath = BosonBathRWA(a_s, η_absorb, γ_absorb, η_emit, γ_emit)




