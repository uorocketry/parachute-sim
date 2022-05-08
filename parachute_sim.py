import math
import matplotlib.pyplot as plt

# Environment Paramaters
rho_init = 1.2  # initial air density [kg/m^3]
g = 9.81        # accel due to gravity [m/s^2]
def rho(alt):   # air density approximation (valid for troposphere)
    return rho_init*math.e**(-alt/1040)

# Rocket Parameters
m = 45.0        # rocket dry mass [kg]
v_init = 0.0    # initial velocity [m/s]
x_init = 8748.0 # initial altitude [m]

# Drogue Parameters
td_d = 0        # deploy at apogee (time zero)
#xd_d = 8748.0  # deploy at apogee (initial altitude)
vt_d = -40.0    # target velocity [m/s]
cd_d = 1.1      # drag coefficient
ot_d = 0.2      # opening time [s] (empirical value from testing)
# calculated area [m^2]
# size chute for target velocity at specified altitude rho(alt)
a_d = (2*g*m)/(rho(1000)*cd_d*vt_d**2)
print('drogue diameter: ', round(2*math.sqrt(a_d/math.pi), 2), 'm', sep='')

# Full Main (full) Parameters
xd_f = 400.0    # deploy altitude [m]
vt_f = -7.0     # target velocity [m/s]
cd_f = 0.8      # drag coefficient
ot_f = 0.2      # opening time (reefed to full) [s]
# calculated area [m^2]
a_f = (2*g*m)/(rho(0)*cd_f*vt_f**2)
print('full-main diameter: ', round(2*math.sqrt(a_f/math.pi), 2), 'm', sep='')

# Partial Main (reefed) Parameters
a_r = a_f       # same a full main [m^2]
xd_r = 800.0    # deploy altitude [m]
vt_r = -18.0    # target velocity [m/s]
#cd_r = 0.1     # drag coefficient
ot_r = 0.2      # opening time [s]
# calculated area [m^2]
#a_r = (2*g*m)/(rho(450)*cd_r*vt_r**2)
#print('partial-main diameter: ', round(2*math.sqrt(a_r/math.pi), 2), 'm', sep='')
cd_r = (2*g*m)/(rho(450)*a_r*vt_r**2)
print('partial-main Cd: ', round(cd_r, 2), 'm', sep='')

# Simulation
dt = 0.001      # time step [s]
t = [0.0]       # simulation time [s] (time from apogee)
y = [x_init]    # rocket altitude [m]
v = [v_init]    # rocket velocity [m/s]
a = [g]         # rocket acceleration [m/s^2]
fd = [0.0]      # parachute drag force [N]

# parachute opening is modeled as constant acceleration process based on a known opening time
# drogue openness
o_d = 0.0
ov_d = 0.0
oa_d = 1.0/(ot_d**2) if ot_d > 0 else -1
# reefed openness
o_r = 0.0
ov_r = 0.0
oa_r = 1.0/(ot_r**2) if ot_r > 0 else -1
# full openness
o_f = 0.0
ov_f = 0.0
oa_f = 1.0/(ot_f**2) if ot_f > 0 else -1

while y[-1] > 0:
    if t[-1] > td_d: # drogue openness (opens near apogee and remains open)
        if o_d < 1.0 and oa_d > 0:
            ov_d += oa_d * dt
            o_d += ov_d * dt
        else:
            o_d = 1.0

    if y[-1] < xd_f: # full main openness
        if o_f < 1.0 and oa_f > 0:
            ov_f += oa_f * dt
            o_f += ov_f * dt
            # transition from reefed to full
            o_r = 1 - o_f
        else:
            o_f = 1.0
            o_r = 0.0

    elif y[-1] < xd_r: # partial main openness
        if o_r < 1.0 and oa_r > 0: # transition from reefed to full
            ov_r += oa_r * dt
            o_r += ov_r * dt
        else:
            o_r = 1.0

    drag_const = cd_d*a_d*o_d + cd_r*a_r*o_r + cd_f*a_f*o_f
    fd.append(0.5*rho(y[-1])*drag_const*v[-1]**2) # parachute force
    a.append((fd[-1] - m*g)/m) # acceleration

    # time step
    v.append(v[-1] + a[-1]*dt)
    y.append(y[-1] + v[-1]*dt)
    t.append(t[-1] + dt)

def annotate_max(var, unit, dec):
    for i in range(1, len(t) - 1): # find and annotate local maxima
        if var[i-1] < var[i] and var[i] > var[i+1]:
            plt.annotate(str(round(var[i], dec)) + unit, (t[i] + 2, var[i]))

plt.plot(t, y)
plt.title('Altitude vs. Time')
plt.ylabel('Altitude (m)')
plt.xlabel('Time (s)')
plt.figure()

v = [-vel for vel in v]
plt.plot(t, v)
plt.title('Velocity vs. Time')
annotate_max(v, ' m/s', 0)
plt.ylabel('Velocity (m/s)')
plt.xlabel('Time (s)')
plt.figure()

plt.plot(t[1:], a[1:])
plt.title('Acceleration vs. Time')
annotate_max(a, ' m/s^2', 1)
plt.ylabel('Altitude (m)')
plt.xlabel('Time (s)')
plt.figure()

#plt.plot(t, a)
#plt.title('Acceleration vs. Time')
#plt.ylabel('Acceleration (m/s^2)')
#plt.xlabel('Time (s)')
#plt.figure()

#e_kinetic = [0.5*m*vel**2 * 10E-6 for vel in v]
#e_potential = [m*g*h * 10E-6 for h in x]
#e_total = [k + p for k, p in zip(e_kinetic, e_potential)]
#plt.ylabel('Energy (MJ)')
#plt.xlabel('Time (s)')
#plt.plot(t, e_kinetic, "-r", label="Kinetic")
#plt.plot(t, e_potential, "-g", label="Potential")
#plt.plot(t, e_total, "-b", label="Total")
#plt.legend(loc="upper right")
#plt.figure()

plt.plot(t, fd)
plt.title('Drag Force vs. Time')
annotate_max(fd, ' N', 0)
plt.ylabel('Drag Force (N)')
plt.xlabel('Time (s)')
plt.show()
