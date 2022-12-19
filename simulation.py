import matplotlib.pyplot as plt
from math import log10, floor, sqrt

t = 0.0
TIMESTEP = 0.001

_state_vec = {}
_plot_vec = {}

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

def step():
    global t
    t += TIMESTEP

def plot(y: float, name: str, yunits: str = '', annotate_max: int = -1, series: str = '', color = 'k'):
    plotxy('{} vs. Time'.format(name), t, y, 'Time', 's', name, yunits, annotate_max, series, color)

def plotxy(title: str, x: float, y: float, xlabel: str = '', xunits: str = '', 
        ylabel: str = '', yunits: str = '', annotate_max: int = -1, series: str = '', color = 'k'):
    if title not in _plot_vec:
        _plot_vec[title] = PlotData(title, xlabel, xunits, ylabel, yunits, annotate_max)
    _plot_vec[title].append(x, y, series, color)

def _update_state(f: float, idf: str) -> list[float]:
    if idf not in _state_vec:
        _state_vec[idf] = [None]*3
    state = _state_vec[idf]
    state[0] = state[1]
    state[1] = state[2]
    state[2] = f
    return state

def integrate(f: float, idf: str) -> float:
    state = _update_state(f, idf)
    if (state[0] != None):
        return (TIMESTEP/12)*(5*state[-1] +8*state[-2] -state[-3])
    elif (state[1] != None):
        return (TIMESTEP/2)*(state[-1] + state[-2])
    else:
        return TIMESTEP*state[-1]

def derive(f: float, idf: str) -> float:
    state = _update_state(f, idf)
    if (state[0] != None):
        return (1/(2*TIMESTEP))*(3*state[-1] -4*state[-2] +state[-3])
    elif (state[1] != None):
        return (1/TIMESTEP)*(state[-1] -state[-2])
    else:
        return 0.0

def _an_max(series: PlotData.Series, pd: PlotData):
    for i in range(1, len(series.xs) - 1): # find and annotate local maxima
        if series.ys[i-1] < series.ys[i] and series.ys[i] > series.ys[i+1]:
            plt.annotate('{} {}'.format(round((series.ys[i]), pd.annotate_max), pd.yunits), (series.xs[i] + 2, series.ys[i]))

def _plot_series(series: PlotData.Series):
    if series.color != '':
        plt.plot(series.xs, series.ys, series.color, label=series.name, linewidth=0.75)
    else:
        plt.plot(series.xs, series.ys, label=series.name, linewidth=0.75)

def _scale_axis(max: float):
    magnitude = 10**floor(log10(max))
    dig = max//magnitude
    err = max%(dig*magnitude)
    if dig > 3 or err/magnitude > 0.5:
        # ex. 7432 -> 8000
        return magnitude*(dig+1)
    else:
        # ex. 2432 -> 2500
        return magnitude*(dig+0.5)

def _scale_plots(plot: PlotData):
    x_lims = [0.0, 0.0]
    y_lims = [0.0, 0.0]
    for series in plot.data.values():
        x_lims[0] = min(x_lims[0], min(series.xs))
        x_lims[1] = max(x_lims[1], max(series.xs))
        y_lims[0] = min(y_lims[0], min(series.ys))
        y_lims[1] = max(y_lims[1], max(series.ys))

    y_lims[1] = _scale_axis(y_lims[1])

    plt.xlim(x_lims)
    plt.ylim(y_lims)

def _douglas_peucker(xs: list[float], ys: list[float], epsilon: float):
    dmax = 0
    imax = 0
    for i in range(2, len(xs)):
        x1 = xs[0]
        x2 = xs[-1]
        y1 = ys[0]
        y2 = ys[-1]
        dx = x2-x1
        dy = y2-y1
        d = abs((y1-ys[i])*dx - (x1-xs[i])*dy) / sqrt(dx**2 + dy**2)
        if (d > dmax):
            imax = i
            dmax = d

    if (dmax > epsilon):
        xs_1, ys_1 = _douglas_peucker(xs[:imax], ys[:imax], epsilon)
        xs_2, ys_2 = _douglas_peucker(xs[imax:], ys[imax:], epsilon)
        return ((xs_1+xs_2), (ys_1+ys_2))
    else:
        return ([xs[0], xs[-1]], [ys[0], ys[-1]])

def _decimate(data: PlotData.Series):
    delta = max(max(data.ys)-min(data.xs), max(data.xs)-min(data.xs))
    epsilon = delta/100000
    old_len = len(data.xs)
    data.xs, data.ys = _douglas_peucker(data.xs, data.ys, epsilon)
    print('old:{}, new:{}, ({}%)'.format(old_len, len(data.xs), round(100*len(data.xs)/old_len, 4)))


def draw_plots():
    for plot in _plot_vec.values():
        plt.title(plot.title)
        plt.gcf().canvas.manager.set_window_title(plot.title)

        if plot.annotate_max >= 0:
            for series in plot.data.values():
                _an_max(series, plot)

        for series in plot.data.values():
            _decimate(series)

        for series in plot.data.values():
            _plot_series(series)

        if plot.xlabel != '' or plot.xunits != '':
            plt.xlabel('{} ({})'.format(plot.xlabel, plot.xunits))
        if plot.ylabel != '' or plot.yunits != '':
            plt.ylabel('{} ({})'.format(plot.ylabel, plot.yunits))
        if len(plot.data) > 1:
            plt.legend(loc="upper right")

        _scale_plots(plot)
        plt.figure()

    plt.close()
    plt.show()
