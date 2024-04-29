import math
import numpy as np

def evaluate_nusselt(Re, Pr, mode="heating"):
    if Re < 2300:
        Nu = 3.66
    elif 2300 <= Re < 10000:
        f = (1.58 * math.log(Re) - 3.28) ** (-2)
        Nu = (0.5 * f * (Re - 1000.) * Pr) / (1 + 12.7 * math.sqrt(0.5 * f) * (Pr**(2/3) - 1))
    else:
        n = 0.3 if mode == "heating" else 0.4
        Nu = 0.023 * Re**0.8 * Pr**n
    return Nu


def cross_section_Rb(θ1, θ2, θ3, σ, β, kgrout):
    M0 = β + math.log(θ2 / (2 * θ1 * (1 - θ1**4)**σ))
    M1 = θ3**2 * (1 - 4 * σ * θ1**4 / (1 - θ1**4))**2 / ((1 + β) / (1 - β) + θ3**2 * (1 + 16 * σ * θ1**4 / (1 - θ1**4)**2))
    return 1 / (4 * math.pi * kgrout) * (M0 - M1)#, M0, M1

def cross_section_Ra(θ1, θ2, θ3, σ, β, kgrout):
    M0 = β + math.log((1 + θ1**2)**σ / (θ3 * (1 - θ1**2)**σ))
    M1 = θ3**2 * (1 - θ1**4 + 4 * σ * θ1**2)**2 / ((1 + β) / (1 - β) * (1 - θ1**4)**2 + θ3**2 * (1 - θ1**4)**2 + 8 * σ * θ1**2 * θ3**2 * (1 + θ1**4))
    return 1 / (math.pi * kgrout) * (M0 - M1)#, M0, M1

