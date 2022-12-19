import math

class Parachute:
    __init__(self, area, deploy)

# Drogue Parameters
td_d = 0        # deploy at apogee (time zero)
#xd_d = 8748.0  # or deploy at altitude
vt_d = -40.0    # target velocity [m/s]
cd_d = 1.0      # drag coefficient
ot_d = 0.0      # opening time [s] (empirical value from testing)
# calculated area [m^2]
# size chute for target velocity at specified altitude rho(alt)
a_d = (2*g*m)/(rho(1000)*cd_d*vt_d**2)
print('drogue diameter: ', round(2*math.sqrt(a_d/math.pi), 2), 'm', sep='')

# Full Main (full) Parameters
xd_f = 350.0    # deploy altitude [m]
vt_f = -7.0     # target velocity [m/s]
cd_f = 0.92     # drag coefficient
ot_f = 0.2      # opening time (reefed to full) [s]
# calculated area [m^2]
a_f = (2*g*m)/(rho(0)*cd_f*vt_f**2)
print('full-main diameter: ', round(2*math.sqrt(a_f/math.pi), 2), 'm', sep='')

# parachute opening is modeled as constant acceleration process based on a known opening time
# drogue openness
o_d = 0.0
ov_d = 0.0
oa_d = 1.0/(ot_d**2) if ot_d > 0 else -1
# full openness
o_f =  0.0
ov_f = 0.0
oa_f = 1.0/(ot_f**2) if ot_f > 0 else -1