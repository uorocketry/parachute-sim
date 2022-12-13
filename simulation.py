import matplotlib.pyplot as plt

state_vec = {}
plot_vec = {}

class PlotData:
    class Series():
        def __init__(self, name, color):
            self.name = name
            self.color = color
            self.xs = []
            self.ys = []
        def append(self, x, y):
            self.xs.append(x)
            self.ys.append(y)

    def __init__(self, title, xlabel, xunits, ylabel, yunits, annotate_max):
        self.data = {}
        self.title = title
        self.xlabel = xlabel
        self.xunits = xunits
        self.ylabel = ylabel
        self.yunits = yunits
        self.annotate_max = annotate_max

    def append(self, x, y, series, color):
        if series not in self.data:
            self.data[series] = PlotData.Series(series, color)
        self.data[series].append(x, y)

def plot(title: str, x: float, y: float, xlabel: str = '', xunits: str = '', ylabel: str = '', yunits: str = '', 
        annotate_max: int = -1, series: str = '', color = ''):
    if title not in plot_vec:
        plot_vec[title] = PlotData(title, xlabel, xunits, ylabel, yunits, annotate_max)
    plot_vec[title].append(x, y, series, color)

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

def an_max(series: PlotData.Series, pd: PlotData):
    for i in range(1, len(series.xs) - 1): # find and annotate local maxima
        if series.ys[i-1] < series.ys[i] and series.ys[i] > series.ys[i+1]:
            plt.annotate('{} {}'.format(round((series.ys[i]), pd.annotate_max), pd.yunits), (series.xs[i] + 2, series.ys[i]))

def plot_series(series: PlotData.Series):
    if series.color != '':
        plt.plot(series.xs, series.ys, series.color, label=series.name, linewidth=0.75)
    else:
        plt.plot(series.xs, series.ys, label=series.name, linewidth=0.75)

def draw_plots():
    for plot in plot_vec.values():
        plt.title(plot.title)
        if plot.xlabel != '' or plot.xunits != '':
            plt.xlabel('{} ({})'.format(plot.xlabel, plot.xunits))
        if plot.ylabel != '' or plot.yunits != '':
            plt.ylabel('{} ({})'.format(plot.ylabel, plot.yunits))

        for series in plot.data.values():
            plot_series(series)
            if plot.annotate_max >= 0:
                an_max(series, plot)

        if len(plot.data) > 1:
            plt.legend(loc="upper right")
        plt.figure()

    plt.close()
    plt.show()

