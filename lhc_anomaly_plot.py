import numpy as np
import matplotlib.pyplot as plt

def vacuum_shockwave_model():
    # Constants
    M_W_SM = 80357.0  # Standard Model W mass in MeV
    M_W_CDF = 80433.5 # CDF II measured mass in MeV
    Anomaly = M_W_CDF - M_W_SM # ~76.5 MeV shift
    
    # Critical Scale: Higgs VEV
    Lambda_crit = 246.0 * 1000 # Convert GeV to MeV
    
    # Energy Range (sqrt(s))
    # From low energy to HL-LHC (14 TeV) and FCC (100 TeV)
    energies_gev = np.logspace(1, 5, 100) # 10 GeV to 100 TeV
    energies_mev = energies_gev * 1000
    
    # Model: Density compression ratio
    # rho/rho_0 = 1 + kappa * (E / Lambda)^n
    # Mass shift delta M / M ~ delta rho / rho
    # We fit kappa to the CDF II result at Tevatron Energy (~1.96 TeV)
    
    E_Tevatron = 1.96 * 1e6 # MeV
    
    # Ansatz: Quadratic dependence on energy density? 
    # Or linear dependence on Energy (Shockwave pressure)?
    # Let's try a "Soft Turn-on" Sigmoid or simple power law for the pre-shock phase.
    # Assuming the "Shock" starts near Lambda_crit.
    
    # Let's use a fluid dynamic analog: Compressibility factor Z
    # delta_M = M_SM * A * (E / Lambda_crit)^2 / (1 + (E/Lambda_crit)^2)
    # This saturates at high E (Vacuum becomes incompressible 'solid'?)
    # Or maybe it diverges?
    
    # Let's try a simple quadratic onset that saturates.
    # delta(E) = A * (E^2) / (E^2 + Lambda^2)
    
    # Fit A to CDF II data
    # 76.5 = 80357 * A * (1.96e6^2) / (1.96e6^2 + 246000^2)
    # Since 1.96 TeV >> 246 GeV, the term (E^2)/(...) is approx 1.
    # This implies A ~ 76.5/80357 ~ 0.00095 (~ 0.1%)
    # This suggests the "Saturation" has already happened at Tevatron?
    # If so, why didn't LEP see it? LEP was at 209 GeV.
    
    # Check LEP point: E = 209 GeV < Lambda (246 GeV).
    # If we use the saturation model:
    # delta(209 GeV) = A * (209^2) / (209^2 + 246^2) ~ A * 0.4
    # Shift ~ 0.4 * 76 MeV ~ 30 MeV.
    # LEP precision was ~30 MeV. So this is consistent!
    
    # Let's define the model
    A = (Anomaly / M_W_SM) * ((E_Tevatron**2 + Lambda_crit**2) / E_Tevatron**2)
    
    def mass_shift(E):
        return M_W_SM * (1 + A * (E**2) / (E**2 + Lambda_crit**2))
    
    # Calculate values
    masses = mass_shift(energies_mev)
    
    # Plotting
    plt.figure(figsize=(10, 6))
    plt.semilogx(energies_gev, masses, label='IGD Vacuum Shockwave Model', color='blue', linewidth=2)
    
    # Reference Lines
    plt.axhline(y=M_W_SM, color='green', linestyle='--', label='Standard Model Prediction')
    plt.axhline(y=M_W_CDF, color='red', linestyle='--', label='CDF II Measurement (Anomaly)')
    
    # Key Points
    plt.scatter([1.96e3], [M_W_CDF], color='red', marker='*', s=150, zorder=5, label='Tevatron (CDF II)')
    plt.scatter([0.209e3], [mass_shift(209000)], color='orange', marker='o', s=100, zorder=5, label='LEP II Limit')
    plt.scatter([13.6e3], [mass_shift(13.6e6)], color='purple', marker='D', s=100, zorder=5, label='LHC Run 3')
    
    # Vertical Line for Lambda_crit
    plt.axvline(x=246.0, color='gray', linestyle=':', label=r'$\Lambda_{crit}$ (Higgs VEV)')
    
    # Annotations
    plt.title(r'W-Boson Mass Shift vs. Collision Energy ($\sqrt{s}$)', fontsize=14)
    plt.xlabel(r'Collision Energy $\sqrt{s}$ (GeV)', fontsize=12)
    plt.ylabel(r'W-Boson Mass (MeV)', fontsize=12)
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend()
    
    # Save plot
    plt.savefig('lhc_anomaly_prediction.png', dpi=300)
    print(f"Plot saved to lhc_anomaly_prediction.png")
    print(f"Model Parameter A (Saturation Amplitude): {A:.6e}")
    print(f"Predicted Mass at LEP (209 GeV): {mass_shift(209000):.2f} MeV")
    print(f"Predicted Mass at LHC (13.6 TeV): {mass_shift(13.6e6):.2f} MeV")

if __name__ == "__main__":
    vacuum_shockwave_model()
