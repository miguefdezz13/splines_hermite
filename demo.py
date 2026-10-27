"""Demo: dibuja curvas cerradas con splines paramétricos y guarda demo.png.

Uso:  python demo.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from splines import plotspline, Resolución_Sistema_LU


def puntos_flor(cx, cy, r, petalos, n):
    t = np.linspace(0, 2 * np.pi, n)
    rr = r * (0.65 + 0.35 * np.cos(petalos * t))
    return np.column_stack((cx + rr * np.cos(t), cy + rr * np.sin(t)))


def puntos_espiral(cx, cy, n):
    t = np.linspace(0, 4 * np.pi, n)
    rr = 0.15 + 0.12 * t
    return np.column_stack((cx + rr * np.cos(t), cy + rr * np.sin(t)))


if __name__ == "__main__":
    # Comprobación rápida del solver LU frente a NumPy
    rng = np.random.default_rng(0)
    A = rng.random((6, 6)) + 6 * np.eye(6)
    b = rng.random(6)
    error = np.max(np.abs(Resolución_Sistema_LU(A, b) - np.linalg.solve(A, b)))
    print(f"Error máximo LU vs numpy.linalg.solve: {error:.2e}")

    curvas = [
        (puntos_flor(2.0, 2.0, 1.5, 5, 21), "gold"),
        (puntos_flor(2.0, 2.0, 0.6, 3, 13), "orangered"),
        (puntos_espiral(6.0, 2.0, 25), "deepskyblue"),
    ]

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.set_facecolor("k")
    for puntos, color in curvas:
        plotspline(puntos, color)
        ax.plot(puntos[:, 0], puntos[:, 1], "o", color="white", ms=3)
    ax.set_aspect("equal")
    ax.set_title("Splines paramétricos (puntos de control en blanco)")
    fig.tight_layout()
    fig.savefig("demo.png", dpi=150)
    print("Guardado demo.png")
