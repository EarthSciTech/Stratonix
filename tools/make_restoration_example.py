#!/usr/bin/env python3
"""Porosity log of examples/restoration_and_backstripping.yaml (stage 4B, STX-TN-15 section 4.5).

The shale trend of Sclater and Christie (1980) (phi0 0.63, c 0.51 1/km) below the H1 unconformity at
W1 (1000 m below the datum), shifted by 700 m of removed section: phi(z) = phi0 exp(-c (z + 700)) with
z the depth below the unconformity. No noise: the eroded-section estimate should return 700 m.
"""
import math
import os

PHI0, C, REMOVED, H1 = 0.63, 0.51e-3, 700.0, 1000.0
out = os.path.join(os.path.dirname(__file__), "..", "examples", "data", "restoration_backstripping", "w1_porosity.csv")
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w") as f:
    f.write("depth_m,porosity\n")
    for k in range(1, 39):
        z = H1 + 10 * k  # below the datum, within Shale_upper (H1 1000 m to H2 1400 m)
        f.write(f"{z:.1f},{PHI0 * math.exp(-C * ((z - H1) + REMOVED)):.6f}\n")
print(out)
