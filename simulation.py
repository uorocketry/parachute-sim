import matplotlib.pyplot as plt

state_vec = {}
plot_vec = {}

class plot_data:
    def __init__(self, title, xlabel, xunits, ylabel, yunits, annotate_max):
        self.xs = []
        self.ys = []
        self.title = title
        self.xlabel = xlabel
        self.xunits = xunits
        self.ylabel = ylabel
        self.yunits = yunits
        self.annotate_max = annotate_max

    def append(self, x, y):
        self.xs.append(x)
        self.ys.append(y)

def plot(title: str, x: float, y: float, xlabel: str = '', xunits: str = '', ylabel: str = '', yunits: str = '', annotate_max: int = -1):
    if title not in plot_vec:
        plot_vec[title] = plot_data(title, xlabel, xunits, ylabel, yunits, annotate_max)
    plot_vec[title].append(x, y)

def update_state(f: float, idf: str) -> list[float]:
    if idf not in state_vec:
        state_vec[idf] = [None]*3
    state = state_vec[idf]
    state[0] = state[1]
    state[1] = state[2]
    state[2] = f
    return state

def integrate(f: float, idf: str, dt: float) -> float:
    state = update_state(f, idf)
    if (state[0] != None):
        return (dt/12)*(5*state[-1] +8*state[-2] -state[-3])
    elif (state[1] != None):
        return (dt/2)*(state[-1] + state[-2])
    else:
        return dt*state[-1]

def derive(f: float, idf: str, dt: float) -> float:
    state = update_state(f, idf)
    if (state[0] != None):
        return (1/(2*dt))*(3*state[-1] -4*state[-2] +state[-3])
    elif (state[1] != None):
        return (1/dt)*(state[-1] -state[-2])
    else:
        return 0.0

def an_max(data: plot_data):
    for i in range(1, len(data.xs) - 1): # find and annotate local maxima
        if data.ys[i-1] < data.ys[i] and data.ys[i] > data.ys[i+1]:
            plt.annotate('{} {}'.format(round((data.ys[i]), data.annotate_max), data.yunits), (data.xs[i] + 2, data.ys[i]))

def draw_plots():
    for plot in plot_vec.values():
        plt.plot(plot.xs, plot.ys)
        plt.title(plot.title)
        if plot.xlabel != '' or plot.xunits != '':
            plt.xlabel('{} ({})'.format(plot.xlabel, plot.xunits))
        if plot.ylabel != '' or plot.yunits != '':
            plt.ylabel('{} ({})'.format(plot.ylabel, plot.yunits))
        if plot.annotate_max >= 0:
            an_max(plot)
        plt.figure()
    plt.close()
    plt.show()

