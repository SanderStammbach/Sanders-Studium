import QuantumToolbox
using HierarchicalEOM
import Plots
using LaTeXStrings
using LinearSolve # to change the solver for better GPU performance
using CUDA
CUDA.allowscalar(false) # Avoid unexpected scalar indexing


W=0.1

#ϵ = -5
#U = 10
σm = sigmam() ## σ-
σz = sigmaz() ## σz
II = qeye(2)  ## identity matrix
λ = 0.1          # Reorganisationsenergie
γ = 0.1          # Zerfallsrate  
N = 14     
kT=0.1
β = 1/kT ## Inverse Temperatur (1/kT)
# construct the annihilation operator for both spin-up and spin-down
# (utilize Jordan–Wigner transformation)
ω0 = 1
σz = sigmaz()
σx = sigmax() 
H0 = 0.5 * ω0 * σz
Td=1.5

# Define the operator that measures the 0, 1 element of density matrix
ρ01 = Qobj([0 1; 0 0])

ψ0 = (basis(2, 0) + basis(2, 1)) / √1;

P00 = ket2dm(basis(2, 0))
P11 = ket2dm(basis(2, 1))
op=sigmax()

function Boson_DrudeLorentz_Matsubara_self(op, λ, W, kT, N)
    

    function matsubara(n)
        return 2 * π * n / β
    end
    
    # Liste der Matsubara-Frequenzen erstellen
    ϵ = [matsubara(n) for n in 1:N]

    η = ComplexF64[λ*cosh(W*Td/2)^2 * W * (cot(W * β / 2.0) - 1.0im)]
    γ = ComplexF64[W]

    if N > 0
        for l in 1:N
            append!(η, 4 * λ*cosh(ϵ[l]*Td/2)^2 * W * ϵ[l] * (kT^2) / (((ϵ[l] * kT)^2) - W^2)) 
            append!(γ, ϵ[l] * kT)
        end
    end

    #δ = _boson_drude_lorentz_approx_discrepancy(λ, W, kT, N, η, γ)

    return BosonBath(op, η, γ)
end


bath=Boson_DrudeLorentz_Matsubara_self(op, λ, ω0, kT, N)
tier = 10
#M = M_Boson(H0, tier, bath)
M_even_cpu = M_Boson(H0, tier, bath)
M_even_cpu =addTerminator(M_even_cpu,bath) #with same coupling as λ?
M_even_gpu = cu(M_even_cpu)
tlist = 0:0.0002:2
noPulseSol = HEOMsolve(M_even_gpu, ψ0, tlist; e_ops = [ ρ01]);

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

#display(Plots.plot(
#    tlist,
#    [abs(noPulseSol.expect[1, :]*noPulseSol.expect[1, :])],
 #   label=["bla"],
#   #linestyle=[:solid],
#    linewidth=3,
 #   xlabel=L"t",
 #   ylabel=L"\rho_{01}",
#    grid=false
#))

Plots.plot(
    tlist,
    [real(noPulseSol.expect[1, :])],  # Beispiel für zwei Kurven
    label=["Kurve 1" "Kurve 2"],
    linestyle=[:solid :dash],  # Jeweils ein Stil pro Kurve
    linewidth=3,
    xlabel=L"t",
    ylabel=L"\rho_{01}",
    grid=false
)










