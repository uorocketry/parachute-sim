import matplotlib.pyplot as plt
from math import log10, floor, sqrt
import csv

t = 0.0
TIMESTEP = 0.001

_state_vec = {}
_plot_vec = {}
_csv_vec = {}

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

def csv_line(file: str, *vars):
    if file not in _csv_vec:
        _csv_vec[file] = []
        _csv_vec[file].append(vars)
    else:
        _csv_vec[file].append(vars)

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

def differentiate(f: float, idf: str) -> float:
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
            plt.annotate('{} {}'.format(round((series.ys[i])), pd.yunits), (series.xs[i] + 2, series.ys[i]))

def _plot_series(series: PlotData.Series):
    if series.color != '':
        plt.plot(series.xs, series.ys, series.color, label=series.name, linewidth=0.75)
    else:
        plt.plot(series.xs, series.ys, label=series.name, linewidth=0.75)

def _scale_axis(max: float):
    magnitude = 10**floor(log10(max)) # max rounded down to nearest power of 10
    dig = max//magnitude # first digit of max
    err = max%(dig*magnitude)
    if dig > 3 or err/magnitude > 0.5:
        # ex. 7432 -> 8000
        if err/magnitude < 0.8: # ensure clearance for annotations
            return magnitude*(dig+1)
        else:
            return magnitude*(dig+2)
    else:
        # ex. 2432 -> 2500
        if err/magnitude < 0.4: # ensure clearance for annotations
            return magnitude*(dig+0.5)
        else:
            return magnitude*(dig+1)

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

def _perp_dist(xs: list[float], ys: list[float]):
    dmax = 0
    imax = 0
    x1 = xs[0]
    x2 = xs[-1]
    y1 = ys[0]
    y2 = ys[-1]
    dx = x2-x1
    dy = y2-y1
    for i in range(2, len(xs)):
        d = abs((y1-ys[i])*dx - (x1-xs[i])*dy) / sqrt(dx**2 + dy**2)
        if (d > dmax):
            imax = i
            dmax = d
    return imax, dmax

def _decimate_douglas_peucker(data: PlotData.Series, epsilon: float):
    xs_stack = [data.xs]
    ys_stack = [data.ys]
    xs_out = []
    ys_out = []

    while len(xs_stack) > 0:
        xs = xs_stack.pop()
        ys = ys_stack.pop()
        imax, dmax = _perp_dist(xs, ys)
        if (dmax > epsilon):
            xs_stack.append(xs[imax:])
            ys_stack.append(ys[imax:])
            xs_stack.append(xs[:imax])
            ys_stack.append(ys[:imax])
        else:
            xs_out.append(xs[0])
            ys_out.append(ys[0])

    xs_out.append(data.xs[-1])
    ys_out.append(data.ys[-1])
    data.xs = xs_out
    data.ys = ys_out

def _decimate_preserve_maxima(data: PlotData.Series, resolution: float):
    xs = data.xs
    ys = data.ys

    e_x = (max(xs) - min(xs)) / resolution
    e_y = (max(ys) - min(ys)) / resolution

    xs_out = [xs[0]]
    ys_out = [ys[0]]
    # sign_x = xs[0] < xs[1]
    # sign_y = ys[0] < ys[1]

    # slope direction change tracking
    is_increasing = ys[1] > ys[0]

    for i in range(1, len(xs) - 1):
        dx = xs[i] - xs_out[-1] 
        dy = ys[i] - ys_out[-1] 

        # local max or min idenification
        is_local_max = ys[i] > ys[i - 1] and ys[i] > ys[i + 1]
        is_local_min = ys[i] < ys[i - 1] and ys[i] < ys[i + 1]

        # determination of if the point is being kept: significant dx/dy, local max/min/ or change in slope direction
        if (abs(dx) > e_x or abs(dy) > e_y or is_local_max or is_local_min or
                (is_increasing != (ys[i] < ys[i + 1]))):
            is_increasing = ys[i] < ys[i + 1]  # directional update 
            xs_out.append(xs[i])
            ys_out.append(ys[i])

    xs_out.append(xs[-1])
    ys_out.append(ys[-1])
    data.xs = xs_out
    data.ys = ys_out

def draw_plots():
    for plot in _plot_vec.values():
        for series in plot.data.values():
            old_len = len(series.xs)
            _decimate_preserve_maxima(series, 1000)
            #delta = min(max(series.ys)-min(series.ys), max(series.xs)-min(series.xs))
            #_decimate_douglas_peucker(series, delta/100000)
            print('old:{}, new:{}, ({}%)'.format(old_len, len(series.xs), round(100*len(series.xs)/old_len, 4)))

    for plot in _plot_vec.values():
        plt.title(plot.title)
        plt.gcf().canvas.manager.set_window_title(plot.title)

        for series in plot.data.values():
            _plot_series(series)
            if plot.annotate_max >= 0:
                _an_max(series, plot)

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

def write_csvs():
    for file in _csv_vec.keys():
        if not file.lower().endswith('.csv'):
            file += '.csv'
        with open(file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            for line in _csv_vec[file]:
                writer.writerow(line)