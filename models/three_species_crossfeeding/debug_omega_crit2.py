#!/usr/bin/env python3
"""Debug script to understand omega_crit2 behavior"""

import numpy as np
import matplotlib.pyplot as plt

class ThreeSpeciesModel:
    """Simplified model for debugging"""

    def __init__(self, r_S=1.0, r_M=0.8, sigma_MS=1.5, sigma_SM=0.5,
                 sigma_GS=0.4, sigma_GM=0.4, sigma_SG=0.4, sigma_MG=0.4,
                 alpha_GS=0.3, alpha_GM=0.3, alpha_SG=0.3, alpha_MG=0.3):
        self.r_S = r_S
        self.r_M = r_M
        self.sigma_MS = sigma_MS
        self.sigma_SM = sigma_SM
        self.sigma_GS = sigma_GS
        self.sigma_GM = sigma_GM
        self.sigma_SG = sigma_SG
        self.sigma_MG = sigma_MG
        self.alpha_GS = alpha_GS
        self.alpha_GM = alpha_GM
        self.alpha_SG = alpha_SG
        self.alpha_MG = alpha_MG

    def net_parameters(self, omega):
        a = (1 - omega) * self.sigma_SG - omega * self.alpha_SG
        c = (1 - omega) * self.sigma_GS - omega * self.alpha_GS
        d = 2 * omega - 1
        return a, c, d

    def SG_equilibrium(self, omega):
        """S-G equilibrium when M is absent"""
        a, c, d = self.net_parameters(omega)

        denominator = 1 - a * c

        if abs(denominator) < 1e-10:
            return None, None

        s_star_SG = ((1 - omega) + a * d) / denominator
        g_star_SG = d + c * s_star_SG

        # Only return if both positive
        if s_star_SG > 0 and g_star_SG > 0:
            return s_star_SG, g_star_SG
        else:
            return None, None

    def invasion_fitness_M(self, omega):
        """M invasion fitness into S-G equilibrium"""
        s_star_SG, g_star_SG = self.SG_equilibrium(omega)

        if s_star_SG is None:
            return np.nan

        lambda_M = -self.r_M + self.sigma_MS * s_star_SG + self.sigma_MG * g_star_SG

        return lambda_M


# Create model and scan omega
model = ThreeSpeciesModel()

omega_range = np.linspace(0.01, 0.99, 500)
s_SG_values = []
g_SG_values = []
lambda_M_values = []

for om in omega_range:
    s, g = model.SG_equilibrium(om)
    lam = model.invasion_fitness_M(om)

    s_SG_values.append(s if s is not None else np.nan)
    g_SG_values.append(g if g is not None else np.nan)
    lambda_M_values.append(lam)

# Plot
fig, axes = plt.subplots(2, 1, figsize=(10, 8))

# Panel 1: S-G equilibrium densities
axes[0].plot(omega_range, s_SG_values, '-', label='$s^*_{SG}$', color='#1f77b4', linewidth=2)
axes[0].plot(omega_range, g_SG_values, '-', label='$g^*_{SG}$', color='#2ca02c', linewidth=2)
axes[0].set_ylabel('Equilibrium density', fontsize=10)
axes[0].set_title('S-G Equilibrium Densities vs ω', fontweight='bold')
axes[0].legend()
axes[0].grid(alpha=0.3)
axes[0].set_xlabel('ω', fontsize=10)

# Panel 2: M invasion fitness
axes[1].plot(omega_range, lambda_M_values, '-', color='#d62728', linewidth=2)
axes[1].axhline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
axes[1].set_ylabel('$\\lambda_M$ (M invasion fitness)', fontsize=10)
axes[1].set_xlabel('ω', fontsize=10)
axes[1].set_title('M Invasion Fitness into S-G Equilibrium', fontweight='bold')
axes[1].grid(alpha=0.3)

# Find zero crossings
zero_crossings = []
for i in range(len(lambda_M_values) - 1):
    if not np.isnan(lambda_M_values[i]) and not np.isnan(lambda_M_values[i+1]):
        if lambda_M_values[i] > 0 and lambda_M_values[i+1] < 0:
            zero_crossings.append((omega_range[i] + omega_range[i+1]) / 2)
            axes[1].axvline((omega_range[i] + omega_range[i+1]) / 2,
                           color='red', linestyle=':', linewidth=2)
            axes[1].text((omega_range[i] + omega_range[i+1]) / 2,
                        axes[1].get_ylim()[1]*0.9,
                        f'ω_crit2 ≈ {(omega_range[i] + omega_range[i+1]) / 2:.3f}',
                        rotation=90, va='top', fontsize=9)

plt.tight_layout()
plt.savefig('debug_omega_crit2.png', dpi=200, bbox_inches='tight')
print("✓ Debug figure saved: debug_omega_crit2.png")

print("\nDiagnostics:")
print(f"λ_M min: {np.nanmin(lambda_M_values):.4f}")
print(f"λ_M max: {np.nanmax(lambda_M_values):.4f}")
print(f"Zero crossings found: {len(zero_crossings)}")
if zero_crossings:
    for i, omega_c in enumerate(zero_crossings):
        print(f"  ω_crit2^({i+1}) = {omega_c:.6f}")
else:
    print("  No zero crossings found - λ_M may not change sign")
    print(f"\n  This means M always {'can' if np.nanmin(lambda_M_values) > 0 else 'cannot'} invade S-G equilibrium")
    print(f"  across the entire ω range, so ω_crit2 does not exist with these parameters.")

plt.show()
