"""Physical constants and default model parameters (SI units unless noted)."""

EPS0 = 8.8541878128e-12   # vacuum permittivity, F/m
F = 96485.33212           # Faraday constant, C/mol
R = 8.314462618           # gas constant, J/(mol K)
T = 298.15                # measurement temperature, K (25 degC)

# Limiting molar ionic conductivities at 25 degC, S m^2 / mol
LAMBDA = {
    "H+": 349.65e-4,
    "Na+": 50.08e-4,
    "NO3-": 71.42e-4,
    "NO2-": 71.8e-4,
    "Cl-": 76.31e-4,
}
