# plot_eigen.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Inputs from your FDM solver
WAVE_CSV = "fdm_output.csv"   # columns: x,psi0,psi1,psi2,...
E_CSV    = "fdm_energies.csv" # columns: index,energy_J

df = pd.read_csv(WAVE_CSV)
E = pd.read_csv(E_CSV).sort_values("index")

x = df["x"].to_numpy()
psi_cols = [c for c in df.columns if c.startswith("psi")]

# Plot first few eigenstates
plt.figure()
for i, col in enumerate(psi_cols[:4]):  # first 4
    psi = df[col].to_numpy()
    # Optional: scale for visualization
    scale = 1.0 / (i + 1)
    plt.plot(x, psi * scale, label=f"{col}  (scaled ×{scale:.2f})")
plt.xlabel("x (m)")
plt.ylabel("ψ_n(x) [scaled]")
plt.title("FDM Eigenstates (shape only)")
plt.legend()
plt.grid(True)
plt.show()

# Print energies
print("Eigen-energies (J):")
for _, row in E.iterrows():
    print(f" n={int(row['index'])}  E={row['energy_J']:.6e} J")
