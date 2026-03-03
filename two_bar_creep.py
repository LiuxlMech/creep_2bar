#!/usr/bin/env python3
"""Two-bar thermo-creep coupling calculator.

Model:
- Bar A: room temperature, elastic only.
- Bar B: high temperature, elastic + thermal + creep.
- Equal axial strain due to rigid connection, zero resultant force at connection.

Usage examples:
    python3 two_bar_creep.py --hours 500
    python3 two_bar_creep.py --hours 500 --dt 0.1 --model norton --K 3.7e-11 --n 3
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, field
from typing import Callable, List, Tuple


# (Temperature °C, alpha 1/°C)
ALPHA_DATA: List[Tuple[float, float]] = [
    (204.0, 1.31e-5),
    (316.0, 1.35e-5),
    (427.0, 1.39e-5),
    (538.0, 1.44e-5),
    (649.0, 1.51e-5),
    (760.0, 1.57e-5),
    (871.0, 1.66e-5),
    (982.0, 1.73e-5),
]

# (Temperature °C, E MPa)
E_DATA: List[Tuple[float, float]] = [
    (20.0, 208000.0),
    (93.0, 201000.0),
    (204.0, 199000.0),
    (316.0, 191000.0),
    (427.0, 184000.0),
    (538.0, 176000.0),
    (649.0, 168000.0),
    (760.0, 157000.0),
    (871.0, 146000.0),
    (982.0, 129000.0),
]


def linear_interp(x: float, table: List[Tuple[float, float]]) -> float:
    """Piecewise linear interpolation with endpoint clamping."""
    if x <= table[0][0]:
        return table[0][1]
    if x >= table[-1][0]:
        return table[-1][1]
    for i in range(len(table) - 1):
        x0, y0 = table[i]
        x1, y1 = table[i + 1]
        if x0 <= x <= x1:
            r = (x - x0) / (x1 - x0)
            return y0 + r * (y1 - y0)
    raise RuntimeError("Interpolation failed unexpectedly.")


def alpha_at(temp_c: float) -> float:
    return linear_interp(temp_c, ALPHA_DATA)


def E_at(temp_c: float) -> float:
    return linear_interp(temp_c, E_DATA)


def integrate_alpha(T0_c: float, T1_c: float, n_steps: int = 2000) -> float:
    """Integrate alpha(T)dT using trapezoidal rule."""
    if T1_c == T0_c:
        return 0.0
    sign = 1.0
    a, b = T0_c, T1_c
    if b < a:
        sign = -1.0
        a, b = b, a

    h = (b - a) / n_steps
    total = 0.0
    for i in range(n_steps + 1):
        T = a + i * h
        w = 0.5 if i == 0 or i == n_steps else 1.0
        total += w * alpha_at(T)
    return sign * total * h


@dataclass
class NortonParams:
    K: float = 3.7e-11  # 1/(h*MPa^n)
    n: float = 3.0


@dataclass
class SimConfig:
    T_a_c: float = 20.0
    T_b_c: float = 800.0
    hours: float = 500.0
    dt_h: float = 0.01
    epsilon_creep0: float = 0.0
    norton: NortonParams = field(default_factory=NortonParams)


@dataclass
class SimResult:
    time_h: List[float]
    eps_creep: List[float]
    sigma_b_mpa: List[float]
    sigma_a_mpa: List[float]
    eps_th_b: float
    E_a_mpa: float
    E_b_mpa: float


class TwoBarCreepSolver:
    def __init__(self, cfg: SimConfig):
        self.cfg = cfg
        self.E_a = E_at(cfg.T_a_c)
        self.E_b = E_at(cfg.T_b_c)
        self.eps_th_b = integrate_alpha(cfg.T_a_c, cfg.T_b_c)
        self.C = 1.0 / (1.0 / self.E_a + 1.0 / self.E_b)

    def sigma_b(self, eps_creep_b: float) -> float:
        # Compression in B is negative.
        return -self.C * (eps_creep_b + self.eps_th_b)

    def creep_rate(self, eps_creep_b: float) -> float:
        sigma = self.sigma_b(eps_creep_b)
        p = self.cfg.norton
        # Generic odd/even n compatible form: eps_dot = K * sigma * |sigma|^(n-1).
        # Compression (sigma<0) gives negative creep strain; tension gives positive.
        return p.K * sigma * (abs(sigma) ** (p.n - 1.0))

    def solve(self) -> SimResult:
        n = int(round(self.cfg.hours / self.cfg.dt_h))
        t = 0.0
        eps_c = self.cfg.epsilon_creep0

        time_h = [t]
        eps_creep = [eps_c]
        sigma_b = [self.sigma_b(eps_c)]
        sigma_a = [-sigma_b[0]]

        # RK4 integration for better stability than forward Euler.
        for _ in range(n):
            dt = self.cfg.dt_h

            k1 = self.creep_rate(eps_c)
            k2 = self.creep_rate(eps_c + 0.5 * dt * k1)
            k3 = self.creep_rate(eps_c + 0.5 * dt * k2)
            k4 = self.creep_rate(eps_c + dt * k3)
            eps_c = eps_c + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

            t += dt
            sb = self.sigma_b(eps_c)

            time_h.append(t)
            eps_creep.append(eps_c)
            sigma_b.append(sb)
            sigma_a.append(-sb)

        return SimResult(
            time_h=time_h,
            eps_creep=eps_creep,
            sigma_b_mpa=sigma_b,
            sigma_a_mpa=sigma_a,
            eps_th_b=self.eps_th_b,
            E_a_mpa=self.E_a,
            E_b_mpa=self.E_b,
        )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Two-bar thermo-creep coupling calculator")
    p.add_argument("--T-a", type=float, default=20.0, help="Bar A temperature (°C)")
    p.add_argument("--T-b", type=float, default=800.0, help="Bar B temperature (°C)")
    p.add_argument("--hours", type=float, default=500.0, help="Simulation duration (h)")
    p.add_argument("--dt", type=float, default=0.01, help="Time step (h), suggest <=0.1 for stability")
    p.add_argument("--K", type=float, default=3.7e-11, help="Norton K, unit 1/(h*MPa^n)")
    p.add_argument("--n", type=float, default=3.0, help="Norton stress exponent n")
    return p


def main() -> None:
    args = build_parser().parse_args()
    cfg = SimConfig(
        T_a_c=args.T_a,
        T_b_c=args.T_b,
        hours=args.hours,
        dt_h=args.dt,
        norton=NortonParams(K=args.K, n=args.n),
    )
    result = TwoBarCreepSolver(cfg).solve()

    print("=== Two-Bar Creep Result ===")
    print(f"E_A = {result.E_a_mpa:.2f} MPa")
    print(f"E_B = {result.E_b_mpa:.2f} MPa")
    print(f"eps_th_B = {result.eps_th_b:.6e}")
    print(f"sigma_B(0) = {result.sigma_b_mpa[0]:.3f} MPa")
    print(f"sigma_B({args.hours:g} h) = {result.sigma_b_mpa[-1]:.3f} MPa")
    print(f"eps_creep_B({args.hours:g} h) = {result.eps_creep[-1]:.6e}")
    if not math.isfinite(result.sigma_b_mpa[-1]) or not math.isfinite(result.eps_creep[-1]):
        print("[WARN] Numerical overflow detected. Please reduce --dt (e.g., 0.01 or 0.001).")


if __name__ == "__main__":
    main()
