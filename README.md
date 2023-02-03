# parachute-sim

A two degree-of-freedom parachute descent simulator used to determine approximate descent time, drift distance, and shock loading based on measured parachute and atmospheric parameters.

### Air Density Approximation

$\rho = \rho_0 e^{-y/10400}$

$\rho_0$: air density measured at ground level
$y$: the altitude in m

### Wind Velocity Profile

$v_w = v_{w10}(y/10)^{\alpha}$

$\alpha$: [Hellmann exponent](https://en.wikipedia.org/wiki/Wind_gradient#Wind_turbines)
$v_{w10}$: wind velocity measured at a standard 10m above ground level
$y$: altitude in m

### Parachute Opening

Parachute opening is modeled as constant acceleration process based on a known opening time. Opening time is determined empirically from deployment testing at predicted airspeed.

The opening process begins at a specified time or altitude. Parachute 'openness' starts at zero and increases quadratically until an 'openness' of one is reached at the specified opening time. Parachute drag is determined by multiplying base parachute drag with the 'openness' factor.

### Integration

quadratic (default)
$x_t = (\Delta t / 12) (5 \dot{x}_{t} + 8 \dot{x}_{t-1} - \dot{x}_{t-2})$

trapezoidal
$x_t = (\Delta t / 2) (\dot{x}_{t} + \dot{x}_{t-1})$

rectangular
$x_t = (\Delta t) (\dot{x}_{t})$
