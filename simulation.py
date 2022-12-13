import matplotlib.pyplot as plt

#INT_ORDER = 3 # 1:Rectangular 2:Linear 3:Quadratic

state_vec = {}
plot_vec = {}

class plot_data:
    def __init__(self, name, xlabel, ylabel):
        self.xs = []
        self.ys = []
        self.name = name
        self.xlabel = xlabel
        self.ylabel = ylabel

    def append(self, x, y):
        self.xs.append(x)
        self.ys.append(y)

def plot(name, x, y, xlabel = None, ylabel = None):
    if name not in plot_vec:
        plot_vec[name] = plot_data(name, xlabel, ylabel)
    plot_vec[name].append(x, y)

def update_state(f: float, idf: str) -> list[float]:
    if idf not in state_vec:
        state_vec[idf] = [None]*3
    state = state_vec[idf]
    state[0] = state[1]
    state[1] = state[2]
    state[2] = f
    #print('updb', state)
    return state

def integrate(f: float, idf: str, dt: float) -> float:
    state = update_state(f, idf)
    if (state[0] != None):
        return (dt/12)*(5*state[-1] +8*state[-2] -state[-3])
    elif (state[1] != None):
        return (dt/2)*(state[-1] + state[-2])
    else:
        return dt*state[-1]

def derive(f: float, dt: float) -> float:
    state = update_state(f)
    if (state[0] != None):
        return (1/(2*dt))*(3*state[-1] -4*state[-2] +state[-3])
    elif (state[1] != None):
        return (1/dt)*(state[-1] -state[-2])
    else:
        return 0.0

def draw_plots():
    for plot in plot_vec.values():
        plt.plot(plot.xs, plot.ys)
        plt.title(plot.name)
        if plot.ylabel != None:
            plt.ylabel(plot.ylabel)
        if plot.xlabel != None:
            plt.xlabel(plot.xlabel)
        plt.figure()
    plt.close()
    plt.show()

