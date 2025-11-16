# animate_time.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Input CSV produced by timeEvolveCrankNicolson
CSV = "cn_time.csv"     # columns: x,t,RePsi,ImPsi,AbsPsi,Phase

df = pd.read_csv(CSV)
# Round t slightly to avoid float noise grouping
df["t_round"] = df["t"].round(12)

# Get sorted unique times
times = np.sort(df["t_round"].unique())

# Build frames: for each time, get x and |psi|^2
frames = []
for tt in times:
    sub = df[df["t_round"] == tt].sort_values("x")
    x = sub["x"].to_numpy()
    prob = (sub["AbsPsi"].to_numpy())**2
    frames.append((tt, x, prob))

# Set up the figure
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlabel("x (m)")
ax.set_ylabel("|ψ(x,t)|²")
ax.set_title("Time evolution")

# Determine fixed axes from all frames for stable view
xmin = min(f[1][0] for f in frames)
xmax = max(f[1][-1] for f in frames)
ymax = max(np.max(f[2]) for f in frames) * 1.05
ax.set_xlim(xmin, xmax)
ax.set_ylim(0, ymax)

time_text = ax.text(0.02, 0.95, "", transform=ax.transAxes)

def init():
    line.set_data([], [])
    time_text.set_text("")
    return line, time_text

def update(i):
    t, x, prob = frames[i]
    line.set_data(x, prob)
    time_text.set_text(f"t = {t:.3e} s")
    return line, time_text

ani = FuncAnimation(fig, update, frames=len(frames), init_func=init, blit=True, interval=60)

plt.show()

# To save as mp4 or gif:
# ani.save("evolution.mp4", fps=30)   # requires ffmpeg installed
# ani.save("evolution.gif", fps=20)   # requires pillow installed
