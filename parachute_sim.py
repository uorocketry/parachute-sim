import math
import matplotlib.pyplot as plt

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
t = [0.0]       # simulation time [s] (time from apogee)
y = [y_init]    # rocket altitude [m]
v = [v_init]    # rocket velocity [m/s]
a = [g]         # rocket acceleration [m/s^2]
fd = [0.0]      # parachute drag force [N]
x = [0]         # rocket horizontal position [m]
# rocket horizontal velocity [m/s]
vx = [wind_velocity(y_init)]

# parachute opening is modeled as constant acceleration process based on a known opening time
# drogue openness
o_d = [0.0]
ov_d = [0.0]
oa_d = [1.0/(ot_d**2) if ot_d > 0 else -1]
# full openness
o_f = [0.0]
ov_f = [0.0]
oa_f = [1.0/(ot_f**2) if ot_f > 0 else -1]

def integrate(f, dt):
    if (len(f) >= 3):   # quadratic
        return (dt/12)*(5*f[-1] +8*f[-2] -f[-3])
    elif (len(f) >= 2): # linear
        return (dt/2)*(f[-1] + f[-2])
    else:               # rectangular
        return dt*f[-1]

def derive(f, dt):
    if (len(f) >= 3):   # quadratic
        return (1/(2*dt))*(3*f[-1] -4*f[-2] +f[-3])
    elif (len(f) >= 2): # linear
        return (1/dt)*(f[-1] -f[-2])
    else:               # not defined
        return 0.0

while y[-1] > 0:
    if t[-1] > td_d: # drogue openness (opens near apogee and remains open)
        if o_d[-1] < 1.0 and oa_d[-1] > 0:
            ov_d.append(ov_d[-1] + integrate(oa_d, dt))
            o_d.append(o_d[-1] + integrate(ov_d, dt))
        else:
            o_d[-1] = 1.0

    if y[-1] < xd_f: # full main openness
        if o_f[-1] < 1.0 and oa_f[-1] > 0:
            ov_f.append(ov_f[-1] + integrate(oa_f, dt))
            o_f.append(o_f[-1] + integrate(ov_f, dt))
        else:
            o_f[-1] = 1.0

    drag_const = cd_d*a_d*o_d[-1] + cd_f*a_f*o_f[-1]
    fd.append(0.5*rho(y[-1])*drag_const*v[-1]**2) # parachute force
    a.append((fd[-1] - m*g)/m) # acceleration

    # time step
    vx.append(wind_velocity(y[-1]))
    v.append(v[-1] + integrate(a, dt))
    y.append(y[-1] + integrate(v, dt))
    x.append(x[-1] + integrate(vx, dt))
    t.append(t[-1] + dt)

def annotate_max(var, unit, dec):
    for i in range(1, len(t) - 1): # find and annotate local maxima
        if var[i-1] < var[i] and var[i] > var[i+1]:
            plt.annotate(str(round(var[i], dec)) + unit, (t[i] + 2, var[i]))

def plot_alt():
    plt.plot(t, y)
    plt.title('Altitude vs. Time')
    plt.ylabel('Altitude (m)')
    plt.xlabel('Time (s)')
    plt.figure()

def plot_vel():
    v_neg = [-vel for vel in v]
    plt.plot(t, v_neg)
    plt.title('Velocity vs. Time')
    annotate_max(v_neg, ' m/s', 0)
    plt.ylabel('Velocity (m/s)')
    plt.xlabel('Time (s)')
    plt.figure()

def plot_accel():
    plt.plot(t[1:], a[1:])
    plt.title('Acceleration vs. Time')
    annotate_max(a, ' m/s^2', 1)
    plt.ylabel('Acceleration (m/s^2)')
    plt.xlabel('Time (s)')
    plt.figure()

def plot_drag_force():
    plt.plot(t, fd)
    plt.title('Drag Force vs. Time')
    annotate_max(fd, ' N', 0)
    plt.ylabel('Drag Force (N)')
    plt.xlabel('Time (s)')
    plt.figure()

def plot_energy():
    e_kinetic = [0.5*m*vel**2 * 10E-6 for vel in v]
    e_potential = [m*g*h * 10E-6 for h in x]
    e_total = [k + p for k, p in zip(e_kinetic, e_potential)]
    plt.ylabel('Energy (MJ)')
    plt.xlabel('Time (s)')
    plt.plot(t, e_kinetic, "-r", label="Kinetic")
    plt.plot(t, e_potential, "-g", label="Potential")
    plt.plot(t, e_total, "-b", label="Total")
    plt.legend(loc="upper right")
    plt.figure()

def plot_drift():
    plt.plot(x, y)
    plt.title('Altitude vs. Horizontal Position ({}m/s Base Wind Velocity)'.format(vw_10))
    plt.xlabel('Horizontal Position (m)')
    plt.ylabel('Altitude (m)')
    plt.figure()

def plot_wind_profile():
    plt.plot([wind_velocity(x) for x in range(int(y_init))], range(int(y_init)))
    plt.title('Wind Profile ({}m/s Base Wind Velocity)'.format(vw_10))
    plt.ylabel('Altitude (m)')
    plt.xlabel('Wind Velocity (m/s)')
    plt.figure()

plot_alt()
plot_vel()
plot_accel()
plot_drag_force()
#plot_energy()
#plot_wind_profile()
plot_drift()

plt.close()
plt.show()
