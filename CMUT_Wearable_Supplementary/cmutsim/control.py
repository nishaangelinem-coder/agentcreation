"""Closed-loop focal-pressure regulation.

Plant: first-order-plus-dead-time (FOPDT) model of the drive-amplitude to
estimated-focal-pressure path. The effective dead time is half a pulse-
repetition interval (zero-order hold of a once-per-burst estimate at 100 Hz) and
the lag is the burst-averaging filter of the pressure estimator.
Controller: discrete PI with conditional-integration anti-windup and
actuator saturation. The pressure magnitude cannot go negative, so the
plant output is clipped at zero.
"""
import numpy as np


def simulate(Kp=4.0, Ki=55.0, K=1.0, tau=0.060, theta=0.005, dt=1e-3,
             T=2.0, setpoint=12.0, t_step=0.1, dist_t=1.0, dist_gain=0.6,
             u_max=40.0, noise_kpa=0.0, drift=None, rng=None, quant_kpa=0.0):
    """Return t, y (kPa), u, r. Disturbance: plant gain multiplied by
    dist_gain at dist_t (e.g. 0.6 = 40 % coupling loss from gel drying)."""
    n = int(round(T / dt))
    t = np.arange(n) * dt
    nd = int(round(theta / dt))
    y = np.zeros(n)
    ym = np.zeros(n)
    u = np.zeros(n)
    r = np.where(t >= t_step, setpoint, 0.0)
    gain = np.where(t >= dist_t, dist_gain, 1.0) * K
    if drift is not None:
        gain = gain * drift
    integ = 0.0
    a = np.exp(-dt / tau)
    rng = rng or np.random.default_rng(0)
    for i in range(1, n):
        ud = u[i - 1 - nd] if i - 1 - nd >= 0 else 0.0
        y[i] = a * y[i - 1] + (1 - a) * gain[i] * ud
        y[i] = max(y[i], 0.0)
        meas = y[i] + (noise_kpa * rng.standard_normal() if noise_kpa else 0)
        if quant_kpa:
            meas = np.round(meas / quant_kpa) * quant_kpa
        ym[i] = meas
        e = r[i] - meas
        u_unsat = Kp * e + integ + Ki * e * dt
        u_sat = np.clip(u_unsat, 0.0, u_max)
        if u_sat == u_unsat or np.sign(e) != np.sign(u_unsat - u_sat):
            integ += Ki * e * dt
        u[i] = np.clip(Kp * e + integ, 0.0, u_max)
    return t, y, u, r, ym


def metrics(t, y, r, t_step=0.1, dist_t=1.0, band=0.02):
    sp = r[-1]
    m1 = (t >= t_step) & (t < dist_t)
    seg_t, seg_y = t[m1], y[m1]
    over = max(0.0, (seg_y.max() - sp) / sp * 100)
    outside = np.where(np.abs(seg_y - sp) > band * sp)[0]
    ts = seg_t[outside[-1]] - t_step if outside.size else 0.0
    m2 = t >= dist_t
    out2 = np.where(np.abs(y[m2] - sp) > band * sp)[0]
    trec = t[m2][out2[-1]] - dist_t if out2.size else 0.0
    ess = abs(np.mean(y[(t > dist_t - 0.1) & (t < dist_t)]) - sp)
    dip = sp - y[m2].min()
    return dict(settling_s=ts, overshoot_pct=over, ss_error_kpa=ess,
                recovery_s=trec, disturbance_dip_kpa=dip)
