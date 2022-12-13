import math
import time
import matplotlib.pyplot as plt
import simulation as sim

# Environment Paramaters
rho_init = 1.2  # initial air density [kg/m^3]
g = 9.81        # accel due to gravity [m/s^2]
def rho(alt):   # air density approximation (valid for troposphere)
    return rho_init*math.e**(-alt/1040)
vw_10 = 5       # wind velocity at 10 meters [m/s]
vw_alpha = 0.27 # hellmann exponent
def wind_velocity(alt):    # wind velocity at specified altitude
    return vw_10*(max(alt,0)/10)**vw_alpha

# Rocket Parameters
m = 28.0        # rocket dry mass [kg]
v_init = 0.0    # initial velocity [m/s]
y_init = 7620.0 # initial altitude [m]

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

# Simulation
dt = 0.001      # time step [s]
t = 0.0       # simulation time [s] (time from apogee)
y = y_init    # rocket altitude [m]
v = v_init    # rocket velocity [m/s]
a = g         # rocket acceleration [m/s^2]
fd = 0.0      # parachute drag force [N]
x = 0         # rocket horizontal position [m]
# rocket horizontal velocity [m/s]
vx = wind_velocity(y_init)

# parachute opening is modeled as constant acceleration process based on a known opening time
# drogue openness
o_d = 0.0
ov_d = 0.0
oa_d = 1.0/(ot_d**2) if ot_d > 0 else -1
# full openness
o_f =  0.0
ov_f = 0.0
oa_f = 1.0/(ot_f**2) if ot_f > 0 else -1

while y > 0:
    if t > td_d: # drogue openness (opens near apogee and remains open)
        if o_d < 1.0 and oa_d > 0:
            ov_d = ov_d + sim.integrate(oa_d, 'oa_d', dt)
            o_d = o_d + sim.integrate(ov_d, 'ov_d', dt)
        else:
            o_d = 1.0

    if y < xd_f: # full main openness
        if o_f < 1.0 and oa_f > 0:
            ov_f = ov_f + sim.integrate(oa_f, 'oa_f', dt)
            o_f = o_f + sim.integrate(ov_f, 'ov_f', dt)
        else:
            o_f = 1.0

    drag_const = cd_d*a_d*o_d + cd_f*a_f*o_f
    fd = 0.5*rho(y)*drag_const*v**2 # parachute force
    a = (fd - m*g)/m # acceleration

    # time step
    vx = wind_velocity(y)
    v += sim.integrate(a, 'a', dt)
    y += sim.integrate(v, 'v', dt)
    x += sim.integrate(vx, 'vx', dt)
    t += dt

    sim.plot('Altitude vs. Time', t, y, 'Time', 's', 'Altitude', 'm')
    sim.plot('Velocity vs. Time', t, -v, 'Time', 's', 'Velocity', 'm/s', annotate_max = 0)
    sim.plot('Acceleration vs. Time', t, a, 'Time', 's', 'Acceleration', 'm/s^2', annotate_max = 1)
    sim.plot('Drag Force vs. Time', t, fd, 'Time', 's', 'Drag Force', 'N', annotate_max = 0)
    sim.plot('Altitude vs. Horizontal Position ({}m/s Base Wind Velocity)'.format(vw_10), x, y, 'Horizontal Position', 'm', 'Altitude', 'm')

    e_kinetic = 0.5*m*(-v)**2 * 10E-6
    e_potential = m*g*y * 10E-6
    e_total = e_kinetic + e_potential
    sim.plot('Energy vs. Time', t, e_kinetic, 'Time', 's', 'Energy', 'MJ', series='Kinetic', color='-r')
    sim.plot('Energy vs. Time', t, e_potential, series='Potential', color='-g')
    sim.plot('Energy vs. Time', t, e_total, series='Total', color='-b')

for y in range(int(y_init)):
    sim.plot('Wind Profile ({}m/s Base Wind Velocity)'.format(vw_10), wind_velocity(y), y, 'Wind Velocity', 'm/s', 'Altitude', 'm')

sim.draw_plots()
