#!/bin/python3
import numpy as np
from qutip import *
import matplotlib.pyplot as plt

class Redfield:
    def __init__(self, H=None, rho0=None, Q=None, a=None, ad=None, w0=None, gamma=None, lam=None, tlist=None,spectrum=None):
        # Initialisierung der Attribute mit Standardwerten, wenn keine Parameter angegeben wurden
        self.a = a if a is not None else destroy(2)  # Vernichtungsoperator (Standard-Qubit)
        self.ad = ad if ad is not None else self.a.dag()  # Steigerungsoperator
        self._lam = lam if lam is not None else 0.01
        self._H = H if H is not None else w0 * (self.ad * self.a )+ np.sqrt(self._lam) * (self.a + self.a.dag())
        self._rho0 = rho0 if rho0 is not None else qutip.coherent_dm(2, 2.0)
        self.Q = Q if Q is not None else np.sqrt(self._lam)*(self.ad+self.a)
        self._w0 = w0 if w0 is not None else 1
        self._gamma = gamma if gamma is not None else 0.5

        self._tlist = tlist if tlist is not None else np.linspace(0, 10, 100)
        # self.ohmic_spectrum=spectrum if spectrum is not None else self.ohmic_spectrum

        if spectrum is None:
            self.spectrum = self.ohmic_spectrum
        else:
            self.spectrum = spectrum

        # self.a_ops = [[self.a,self.spectrum], [ self.ad, self.spectrum]]  # Liste von Kollapsoperatoren und spektralen Dichten
        self.a_ops=[[Q,self.spectrum]]
    def ohmic_spectrum(self, wk):
        """Ohmsche spektrale Dichte."""
        if wk == 0.0:  # Dephasierungs-induzierendes Rauschen
            return self._lam
        else:  # Relaxations-induzierendes Rauschen
            return (self._lam / 2) * (wk / (2 * np.pi)) * (wk > 0.0)

    def Rsolver(self):
        """Löst das System mit brmesolve."""
        e_ops = [self.a, self.ad]  # Beobachtungsoperatoren
        # Berechnung mit brmesolve
        result_brme = brmesolve(self._H, self._rho0, self._tlist, self.a_ops, e_ops=e_ops)
        return result_brme

    def plot_results(self, result):
        """Zeigt die Ergebnisse der Lösung an."""
        plt.figure()
        for idx, e_op in enumerate(result.expect):
            plt.plot(self._tlist, e_op, label=f"Expectation value {idx + 1}")
        plt.xlabel("Time")
        plt.ylabel("Expectation value")
        plt.legend()
        plt.show()


if __name__ == "__main__":
    # Definition des Hamiltonians
    #H = 0.5 * 2 * np.pi * sigmax()  # Beispiel: ein Qubit mit sigmax Hamiltonian

    # Anfangszustand des Systems
    psi0 = basis(2, 0)  # Qubit im Zustand |0>

    # Relaxationsoperator (z.B. sigmax)
    #O = sigmax()

    # Parameter des Systems
    w0 = 1.0  # Eigenfrequenz
    gamma = 0.1  # Dämpfungsrate
    lam = 1  # Kopplungsstärke
    kappa=0.4**2
    nB=1
    # Zeiten für die Lösung
    tlist = np.linspace(0.0, 30, 100)
    
    def spectrum_c(w):
        if w == 0.0:
            return kappa * nB
        else:
            return kappa * w * (nB + 1) * (w > 0) + kappa * nB * w * (w <= 0)

    # Instanziiere die Klasse
    redfield_solver = Redfield(H=None, rho0=psi0, Q=None, w0=w0, gamma=None, lam=lam, tlist=tlist)

    # Löse das System
    rho_t_redfield = redfield_solver.Rsolver()

    # Plotte die Ergebnisse
    redfield_solver.plot_results(rho_t_redfield)
