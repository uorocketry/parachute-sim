import math
import matplotlib.pyplot as plt
import simulation as sim

# Environment Paramaters
rho_init = 1.2  # initial air density [kg/m^3]
g = 9.81        # accel due to gravity [m/s^2]
def rho(alt):   # air density approximation (valid for troposphere)
    return rho_init*math.e**(-alt/10400)
vw_10 = 0       # wind velocity at 10 meters [m/s]
vw_alpha = 0.27 # hellmann exponent
def wind_velocity(alt):    # wind velocity at specified altitude
    return vw_10*(max(alt,0)/10)**vw_alpha

# Rocket Parameters
m = 11.2        # dry mass [kg]
y_init = 304.80 # expected apogee [m]
vx_init = 50.5 # horizontal velocity at apogee [m/s]

# Drogue Parameters
td_d = 3.0      # deployment delay [s] (0 = deploy at apogee)
vt_d = -35.0    # target velocity [m/s]
cd_d = 1.8      # drag coefficient
ot_d = 0.0      # opening time [s] (empirical value from testing)
# calculated area [m^2]
a_d = (2*g*m)/(rho(1000)*cd_d*vt_d**2) # size parachute for target velocity at specified altitude rho(alt)
#a_d = 0.62**2*math.pi/4  # get parachute area based on diameter
print('drogue diameter: ', round(2*math.sqrt(a_d/math.pi), 2), 'm', sep='')

# Full Main (full) Parameters
xd_f = 450.0    # deploy altitude [m]
vt_f = -7.0     # target velocity [m/s]
cd_f = 2.2      # drag coefficient
ot_f = 0.0      # opening time (reefed to full) [s]
# calculated area [m^2]
a_f = (2*g*m)/(rho(0)*cd_f*vt_f**2)  # size parachute for target velocity at specified altitude rho(alt)
#a_f = 5.0**2*math.pi/4 # get parachute area based on diameter
print('full-main diameter: ', round(2*math.sqrt(a_f/math.pi), 2), 'm', sep='')

# Simulation
x = 0.0         # horizontal position [m]
vx = vx_init    # horizontal velocity [m/s]
y = y_init      # altitude [m]
vy = 0.0        # vertical velocity [m/s]
# parachute opening is modeled as constant acceleration process based on a known opening time
# drogue openness
o_d = 0.0
ov_d = 0.0
oa_d = 1.0/(ot_d**2) if ot_d > 0 else -1
# main openness
o_f =  0.0
ov_f = 0.0
oa_f = 1.0/(ot_f**2) if ot_f > 0 else -1

while y > 0:
    if sim.t > td_d: # drogue openness (opens near apogee and remains open)
        if o_d < 1.0 and oa_d > 0:
            ov_d += sim.integrate(oa_d, 'oa_d')
            o_d += sim.integrate(ov_d, 'ov_d')
        else:
            o_d = 1.0

    if y < xd_f: # main openness
        if o_f < 1.0 and oa_f > 0:
            ov_f += sim.integrate(oa_f, 'oa_f')
            o_f += sim.integrate(ov_f, 'ov_f')
        else:
            o_f = 1.0

    # wind velocity
    ux = -vx + wind_velocity(y)
    uy = -vy
    u = math.hypot(ux, vy) # airspeed

    # drag force
    cd = cd_d*a_d*o_d + cd_f*a_f*o_f
    fd = 0.5*rho(y)*cd*u**2
    fdx = 0 if u == 0 else fd*(ux/u)
    fdy = 0 if u == 0 else fd*(uy/u)

    # acceleration
    ax = fdx/m
    ay = fdy/m - g

    # time step
    vx += sim.integrate(ax, 'ax')
    vy += sim.integrate(ay, 'ay')
    x += sim.integrate(vx, 'vx')
    y += sim.integrate(vy, 'vy')

    # plot
    sim.plot(y, 'Altitude', 'm', annotate_max=0)
    sim.plot(-vy, 'Vertical Velocity', 'm/s', annotate_max=0)
    sim.plot(math.hypot(ax, ay), 'Total Acceleration', 'm/s$^2$', annotate_max = 0)
    sim.plot(ay, 'Vertical Acceleration', 'm/s$^2$', annotate_max = 0)
    sim.plot(fd, 'Total Drag Force', 'N', annotate_max = 0)
    sim.plot(math.hypot(vx, vy), 'Total Velocity', 'm/s', annotate_max = 0)

    #sim.plot(-vy, 'Velocity', 'm/s', series='y', color='-r')
    #sim.plot(vx, 'Velocity', 'm/s', series='x', color='-g')
    #sim.plot(math.hypot(vx, vy), 'Velocity', 'm/s', series='total', color='-b')

    #e_kinetic = 0.5*m*(math.hypot(vx, vy))**2/1000
    #e_potential = m*g*y/1000
    #sim.plot(e_kinetic, 'Energy', 'kJ', series='Kinetic', color='-r')
    #sim.plot(e_potential, 'Energy', 'kJ', series='Potential', color='-g')
    #sim.plot(e_kinetic + e_potential, 'Energy', 'kJ', series='Total', color='-b')

    sim.plotxy('Altitude vs. Horizontal Position ({} m/s Base Wind Velocity)'.format(vw_10), x, y, 'Horizontal Position', 'm', 'Altitude', 'm', annotate_max=0)

    # csv output
    if int(sim.t*(1/sim.TIMESTEP)) % 10 == 0: # 0.01s sample rate
        sim.csv_line('out.csv', round(sim.t, 2), x, y, vx, vy, ax, ay, fd)

    sim.step()

for y in range(int(y_init)):
    sim.plotxy('Wind Profile ({} m/s Base Wind Velocity)'.format(vw_10), wind_velocity(y), y, 'Wind Velocity', 'm/s', 'Altitude', 'm')
    sim.csv_line('wind_profile.csv', y, wind_velocity(y))

print('Writing CSVs')
sim.write_csvs()
print('Plotting')
sim.draw_plots()
