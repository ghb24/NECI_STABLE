#!/usr/bin/python
'''Plot data from and FCIMCStats file, using OUTPUT file as needed

See "plot.py -h" for more details
'''

g_ver_str = '0.2'

# Required modules
from pylab import *
from numpy import *
from matplotlib import rc, rcParams
from os import path
import os
import sys
import matplotlib.axes
try:
    from argparse import ArgumentParser
except:
    from locargparse import ArgumentParser
import re
import time
import subprocess
import tempfile
from scipy import optimize

def _general_function(params, xdata, ydata, function):
    return function(xdata, *params) - ydata

def _weighted_general_function(params, xdata, ydata, function, weights):
    return weights * (function(xdata, *params) - ydata)

# Would be nice if this was here...#
def curve_fit(f, xdata, ydata, p0=None, sigma=None, **kw):
    """
    Use non-linear least squares to fit a function, f, to data.

    Assumes ``ydata = f(xdata, *params) + eps``

    Parameters
    ----------
    f : callable
        The model function, f(x, ...).  It must take the independent
        variable as the first argument and the parameters to fit as
        separate remaining arguments.
    xdata : An N-length sequence or an (k,N)-shaped array
        for functions with k predictors.
        The independent variable where the data is measured.
    ydata : N-length sequence
        The dependent data --- nominally f(xdata, ...)
    p0 : None, scalar, or M-length sequence
        Initial guess for the parameters.  If None, then the initial
        values will all be 1 (if the number of parameters for the function
        can be determined using introspection, otherwise a ValueError
        is raised).
    sigma : None or N-length sequence
        If not None, it represents the standard-deviation of ydata.
        This vector, if given, will be used as weights in the
        least-squares problem.


    Returns
    -------
    popt : array
        Optimal values for the parameters so that the sum of the squared error
        of ``f(xdata, *popt) - ydata`` is minimized
    pcov : 2d array
        The estimated covariance of popt.  The diagonals provide the variance
        of the parameter estimate.

    Notes
    -----
    The algorithm uses the Levenburg-Marquardt algorithm:
    scipy.optimize.leastsq. Additional keyword arguments are passed directly
    to that algorithm.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.optimize import curve_fit
    >>> def func(x, a, b, c):
    ...     return a*np.exp(-b*x) + c

    >>> x = np.linspace(0,4,50)
    >>> y = func(x, 2.5, 1.3, 0.5)
    >>> yn = y + 0.2*np.random.normal(size=len(x))

    >>> popt, pcov = curve_fit(func, x, yn)

    """
    if p0 is None or isscalar(p0):
        # determine number of parameters by inspecting the function
        import inspect
        args, varargs, varkw, defaults = inspect.getargspec(f)
        if len(args) < 2:
            msg = "Unable to determine number of fit parameters."
            raise ValueError(msg)
        if p0 is None:
            p0 = 1.0
        p0 = [p0]*(len(args)-1)

    args = (xdata, ydata, f)
    if sigma is None:
        func = _general_function
    else:
        func = _weighted_general_function
        args += (1.0/asarray(sigma),)

    res = optimize.leastsq(func, p0, args=args, full_output=1, **kw)
    (popt, pcov, infodict, errmsg, ier) = res

    if ier not in [1,2,3,4]:
        msg = "Optimal parameters not found: " + errmsg
        raise RuntimeError(msg)

    if (len(ydata) > len(p0)) and pcov is not None:
        s_sq = (func(popt, *args)**2).sum()/(len(ydata)-len(p0))
        pcov = pcov * s_sq
    else:
        pcov = inf

    return popt, pcov


####################
# A colour management class for plots
####################
class colour_manager:
    def __init__ (self):
        #self.colours = matplotlib.axes._process_plot_var_args.defaultColors
        self.colours = rcParams['axes.color_cycle']
        self.ind = 0
        self.ref_point = 0
        self.frozen = False
        self.ind_max = 0

    def __call__ (self):
        ind = self.ind
        if not self.frozen:
            self.ind = self.ind + 1
        return self.colours[(ind-1) % len(self.colours)]

    def reset_point(self):
        self.ref_point = self.ind

    def freeze (self, state=True):
        self.frozen = state

    def reset (self, value=None):
        self.ind_max = max(self.ind, self.ind_max)
        if value is None:
            self.ind = self.ref_point
        else:
            self.ref_point = value
            self.ind = value

    def jump_max (self):
        self.ind = max(self.ind, self.ind_max)


# Load a file, and read the columns (up to a specified limit)
def read_cols (f, last_iter=None, last_im_time=None):

    # The list of available column titles. Some are duplicates to
    # deal with where the labelling has adjusted, or to cope with
    # both FCIQMCStats and fciqmc_stats conventions.
    max_runs = 20
    re_label = re.compile('^#\s*((\d+\.\s*.*?)\s*)+$')
    re_label_split = re.compile('#?\s*\d+\.')
    col_labels = {
        'iter': ('Step', 'Iter.'),
        'shift': (('Shift. (cyc)', 'Shift') +
                    tuple(('Shift ({0})'.format(run+1) for run in range(max_runs)))),
        'growth': ('GrowRate', 'Growth fac.'),
        'parts': (('TotWalkers', 'Tot. parts') +
                    tuple(('Parts ({0})'.format(run+1) for run in range(max_runs)))),
        'ref_parts': (('NoatHF', 'Tot. ref') +
                        tuple(('Ref ({0})'.format(run+1) for run in range(max_runs)))),
        'im_time': ('TotImagTime', 'Im. time'),
        'it_time': ('Iter. time', "IterTime"),
        'proje_corr': ('Proj.E.ThisCyc', 'Proj. E (cyc)'),
        'proje_tot': (('Tot-Proj.E.ThisCyc', 'Tot. Proj. E') +
                    tuple(('Tot ProjE ({0})'.format(run+1) for run in range(max_runs)))),
        'exp_av_proje': (),
        'av_proje': ('Proj.E',),
        'av_shift': ('Av.Shift',),
        'born': ('NoBorn', 'No. born'),
        'died': ('NoDied', 'No. died'),
        'annihil': ('Annihil', 'No. annihil'),
        'S2': ('Inst S^2', ),
        'S2_init': (),
        'accept': ('AccRat', 'Acc. rate'),
    }
    
    # Which data columns can be repeated in an output files, for multiple run calculations
    multi_run_cols = [
        'shift',
        'parts',
        'ref_parts',
        'proje_tot'
    ]

    # Invert the label dictionary for easy lookup
    label_lookup = dict([(hdr, title) for title in col_labels
                                for hdr in col_labels[title]])

    rows = []
    iter_col = None
    im_time_col = None
    found_hdr = False
    for line in f:

        # Ignore comment lines
        if (line[0] == '#'):

            # The first line contains the column header. Split this up
            # into sections
            if not found_hdr and re_label.match(line):
                found_hdr = True
                col_hdrs =  [s.strip() for s in re_label_split.split(line)
                                       if s != '']
                col_titles = [label_lookup.get(hdr, None) for hdr in col_hdrs]
                try:
                    iter_col = col_titles.index('iter')
                except:
                    pass
                try:
                    im_time_col = col_titles.index('im_time')
                except:
                    pass

                try:
                    s2_col = col_titles.index('S2')
                    if col_titles[s2_col+1] == 'S2':
                        col_titles[s2_col+1] = 'S2_init'
                except:
                    pass

            continue

        bits = [float(s) for s in line.split()]

        # Have we specified an iteration/time not to read after?
        if last_iter and iter_col and bits[iter_col] > last_iter:
            break
        if last_im_time and im_time_col and bits[im_time_col] > last_im_time:
            break

        # Add the row to the store
        rows.append(bits)

    # Flip this around so we have a list of columns to return
    cols = zip(*rows)

    # And push the data into our custom object
    data = column_data()
    for i, title in enumerate(col_titles):
        if title is not None:
            if title in multi_run_cols:
                if title not in data:
                    data[title] = []
                data[title].append(cols[i])
            else:
                data[title] = cols[i]

    return data



def output_file (fin):
    '''Given an FCIMCStats file, determine the output file'''

    # Process the filename
    stem = path.dirname(fin)
    stem = stem + '/' if stem != '' else stem
    if 'FCIMCStats' in fin:
        parts = [x.strip('.') for x in path.basename(fin).partition('FCIMCStats')]
    elif 'FCIQMCStats_' in fin:
        parts = [x.strip('.') for x in path.basename(fin).partition('FCIQMCStats_')]
    elif 'fciqmc_stats' in fin:
        parts = [x.strip('.') for x in path.basename(fin).partition('fciqmc_stats')]
    else:
        parts = [x.strip('.') for x in path.basename(fin).partition('FCIQMCStats')]

    # File bases for output files
    out_base_in = ['OUT', 'OUTPUT', 'OUT.FCIMC', 'OUTPUT.FCIMC', 'out', 'output.FCIMC', 'output', 'out_', 'neci.out', 'neci.out_']
    out_base = []
    map(out_base.extend, [(a, '.'+a, a+'.', '.'+a+'.') for a in out_base_in])

    # Generate list of files to try
    out_names = []
    for f in out_base:
        out_names.append(stem + parts[0] + f + parts[2])
        out_names.append(stem + parts[2] + f + parts[0])
    out_names.extend(out_base_in)

    # Does any of these possible filenames exist?
    for f in out_names:
        if path.isfile(f):
            return f
    return None



#
# A data store
# --> Something dict-like that receives read-in data. This can then be
#     accessed during the run.
class column_data(dict):
    '''
    This basically acts as a dict, but if some of the data is not
    available, then it returns an iterator of the correct length so
    that we don't fuck up the plots.
    '''

    def __getitem__(self, key):
        '''
        If we have stored the column values, then return them. Otherwise
        we need to return a null column of the correct length, which will
        be plotted as 'no-data'
        '''

        try:
            val = super(column_data, self).__getitem__(key)
        except KeyError:
            # If we haven't defined 'iter' then this will throw again.
            iters = super(column_data, self).__getitem__('iter')
            val = zeros_like(iters)

        return val



####################
# Our global plotter object. Does all of the work
####################
class plotter:

    class in_file:
        def __init__ (self, fn=None):

            # Are we using the default files?
            if fn is None:
                if path.exists('fciqmc_stats'):
                    fn = 'fciqmc_stats'
                elif path.exists('FCIMCStats'):
                    fn = 'FCIMCStats'
                elif path.exists('FCIQMCStats'):
                    fn = 'FCIQMCStats'
                else:
                    fn = 'FCIQMCStats_'

            # Is this a file on a cluster which needs mounting?
            self.tmp_dir = None
            if ':' in fn:

                # Mount the relevant files.
                (host, jobid) = fn.split(':')
                self.tmp_dir = tempfile.mkdtemp()
                mount_script = "%s-job-mount" % host
                subprocess.call([mount_script, jobid, self.tmp_dir])

                # Find the output file to use!
                if path.exists("%s/FCIMCStats" % self.tmp_dir):
                    fn = "%s/FCIMCStats" % self.tmp_dir
                elif path.exists("%s/FCIQMCStats" % self.tmp_dir):
                    fn = "%s/FCIQMCStats" % self.tmp_dir
                else:
                    fn = "%s/FCIQMCStats_" % self.tmp_dir


            self.fn = fn
            self.ofile = output_file(fn)


        def __del__ (self):
            '''Clean up any temporary filesystem objects created'''

            # If this is being called as part of cleanup code, then global
            # variables may have been destroyed --> modules may not exist
            # any more. Just reload them in case.
            import subprocess, os

            if self.tmp_dir:
                subprocess.call(["fusermount", "-u", self.tmp_dir])
                os.rmdir(self.tmp_dir)



    def __init__ (self):
        '''Default arguments'''

        class ax_data:
            ax = None
            plot = False
            legend_pos = 1
            shared = False
            label_left = True
            log_scale = False
            label = ''
            limits = None

        # Storage for the axes objects
        self.axes = []
        self.fig = None

        # Energy plot options
        self.E = ax_data()
        self.E.plot = True
        self.E.tgt_E = []
        self.E.label = 'Energy / $E_h$'
        self.E.plot_averages = False
        self.E.total_energies = True
        self.E.plot_lines = False
        self.E.plot_exp_average = False
        self.axes.append(self.E)

        # Walker plot options
        self.W = ax_data()
        self.W.plot = True
        self.W.log_scale = True
        self.W.plot_growth_cpts = False
        self.W.label = 'No. of Walkers'
        self.share_walkers_E = False
        self.axes.append(self.W)

        # Spin plot options
        self.S = ax_data()
        self.S.plot = False
        self.S.log_scale = False
        self.S.label = '$S^2$'
        self.share_walkers_spin = True
        self.axes.append(self.S)

        # Time plot options
        self.T = ax_data()
        self.T.plot = False
        self.T.log_scale = False
        self.T.plot_tau = False
        self.T.plot_iter_time = False
        self.T.label = 'time / s (or a.u.)'
        #self.share_walkers_time = True
        self.axes.append(self.T)

        # Survival plots
        self.U = ax_data()
        self.U.plot = False
        self.U.log_scale = False
        self.U.label = 'Spawning survival ratio'
        self.axes.append(self.U)

        # Growth rate comparisons
        self.G = ax_data()
        self.G.plot = False
        self.G.label = "Growth rate ratio"
        self.axes.append(self.G)

        # Plot data
        self.in_files = [plotter.in_file()]
        self.last_iter = None
        self.x_time = False
        self.x_itime = False
        self.x_walkers = False
        self.fit_walkers_E = False
        self.fix_colours = False
        self.plot_legend = True
        self.legend_strs = None
        self.verbose = False

        # Storage to keep track of temporary directories which
        # need removal.
        self.temp_mount_dirs = []





    def proc_args (self):
        '''Read in the command line options. For default arguments see
           __init__'''

        parser = ArgumentParser(description="Plot data from and "\
                            "FCIMCStats file using OUTPUT file as needed")

        parser.add_argument("in_files", metavar='file',
                            default=[f.fn for f in self.in_files],
                            help="FCIMCStats files to plot.", nargs='*')

        parser.add_argument("-S", "--plot-S2", default=self.S.plot,
                            help="Plot the instantaneous total spin (S^2)",
                            action='store_true')

        parser.add_argument("-l", "--linear", default=(not self.W.log_scale),
                            help="Plot walkers on a linear scale",
                            action='store_true')

        parser.add_argument("-t", "--time", default=self.x_time,
                            help="Plot using cumulative time on the x-axis",
                            action='store_true')

        parser.add_argument("-T", "--imag-time", default=self.x_itime,
                            help="Plot using elapsed imaginary time on the "\
                            "x-axis", action='store_true')

        parser.add_argument("--spin-log", default=self.S.log_scale,
                            help="Plot the instantaneous spin on a log scale",
                            action='store_true')

        parser.add_argument("--tgt-E", default=self.E.tgt_E,
                            help="A list of target energies, or a file "\
                                "containing such a list (as its first "\
                                "column)", action='append')

        parser.add_argument("-L", "--no-legends", default=not self.plot_legend,
                            help="Disable printing of legends on the plots",
                            action='store_true')

        parser.add_argument("-I", "--last-iter", default=self.last_iter,
                            help="Last iteration (or imag time) to include "\
                                 "in the plot", type=int)

        parser.add_argument("-W", "--walker-x", default=self.x_walkers, 
                            help="Use the walker number as the x-axis",
                            action='store_true')

        parser.add_argument("--no-walkers", default=not self.W.plot,
                            help="Don't plot any walker-number details",
                            action='store_true')

        parser.add_argument("-f", "--fit-walkers-E", action='store_true', 
                            default=self.fit_walkers_E,
                            help="Fit a stretched exponential to the energy "\
                                 "and walker data.")

        parser.add_argument("-c", "--correlation-E", action='store_true',
                            default=not self.E.total_energies,
                            help="Use the correlation energy directly. Don't"\
                                 " convert into total energies.")

        parser.add_argument("-a", "--average-E", action='store_true',
                            default=self.E.plot_averages,
                            help="Plot the cumulative average projected "\
                                 "energy and shift.")

        parser.add_argument("-C", "--fix-colours", action='store_true',
                            default=self.fix_colours,
                            help="Use the same line colours for all the "\
                                 "lines from the same file")

        parser.add_argument("-E", "--E-range", default=self.E.limits,
                            help="Specify maximum and minimum energies "\
                                 "for plot.", type=float, nargs=2)

        parser.add_argument("-V", "--verbose", action='store_true',
                            default=self.verbose,
                            help="Use verbose output")

        parser.add_argument("-v", "--version", action='version',
                            version='%%(prog)s %s' % g_ver_str)

        parser.add_argument("--separate-spin", action='store_true',
                            default=not self.share_walkers_spin,
                            help="Plot S^2 on a different set of axes to "\
                                 "the number of walkers")

        parser.add_argument("--plot-tau", action='store_true',
                            default=self.T.plot_tau,
                            help="Plot the live value of tau")

        parser.add_argument("--plot-iter-time", action='store_true',
                            default=self.T.plot_iter_time,
                            help="Plot the live time-per-iteration")

        parser.add_argument("--plot-ref-rel-growrate", action='store_true',
                            default=self.G.plot,
                            help="Plot the relative growth ratio of all "\
                                 "walkers compared to the reference")

        parser.add_argument("--time-log", action='store_true',
                            default=self.T.log_scale,
                            help="Plot tau or iteration times on a log scale")

        parser.add_argument("-G", "--plot-growth-cpts", action='store_true',
                            default=self.W.plot_growth_cpts, 
                            help="Plot the number of walkers born, died and "\
                                 "annihilated")
        
        parser.add_argument("--plot-acc-ratio", action='store_true',
                            default=self.U.plot,
                            help="Plot the acceptance/survival ratio")

        parser.add_argument("--plot-energy-lines", action='store_true',
                            default=self.E.plot_lines,
                            help="Use lines rather than points to plot the "\
                                 "energy values")

        parser.add_argument("--plot-exp-average", action='store_true',
                            default=self.E.plot_exp_average,
                            help="Plot the exponentially weighted average of "\
                                 "the projected energy")

        parser.add_argument("--legend-strs", default=self.legend_strs, 
                            help="What labels should we use in the legend?",
                            type=str, nargs=1)

        args = parser.parse_args()

        # Store obtained data
        self.S.plot = args.plot_S2
        self.W.log_scale = not args.linear
        self.x_time = args.time
        self.x_walkers = args.walker_x
        self.S.log_scale = args.spin_log
        self.E.tgt_E = map(to_float_if_float, args.tgt_E)
        self.plot_legend = not args.no_legends
        self.last_iter = args.last_iter
        self.fit_walkers_E = args.fit_walkers_E
        self.E.total_energies = not args.correlation_E
        self.share_walkers_spin = not args.separate_spin
        self.fix_colours = args.fix_colours
        self.E.limits = args.E_range
        self.E.plot_averages = args.average_E
        self.x_itime = args.imag_time
        self.T.plot_tau = args.plot_tau
        self.T.plot_iter_time = args.plot_iter_time
        self.T.plot = self.T.plot_tau or self.T.plot_iter_time
        self.T.log_scale = args.time_log
        self.G.plot = args.plot_ref_rel_growrate
        self.W.plot = not args.no_walkers
        self.W.plot_growth_cpts = args.plot_growth_cpts
        self.U.plot = args.plot_acc_ratio
        self.E.plot_lines = args.plot_energy_lines
        self.E.plot_exp_average = args.plot_exp_average
        self.verbose = args.verbose

        if args.legend_strs:
            self.legend_strs = args.legend_strs[0].split(',')

        if self.E.limits:
            if self.E.limits[1] < self.E.limits[0]:
                self.E.limits = self.E.limits[1], self.E.limits[0]

        self.in_files = [plotter.in_file(f) \
                               for f in args.in_files]



    def init_axes (self):
        '''Initialise the axes object for later usage'''

        # If we are plotting by walker number, then we don't want to plot
        # the total number of walkers...
        assert not(self.x_walkers and self.x_time)
        if self.x_walkers:
            self.W.plot = False
            self.share_walkers_E = False
            self.share_walkers_spin = False

        # Update the axis data structures to be consistent with applied
        # options.
        # TODO: Update legend position, unless customised.
        if self.share_walkers_E and self.E.plot and self.W.plot:
            self.E.shared = True
            self.E.label_left = True
            self.E.legend_pos = 4
            self.W.shared = True
            self.W.label_left = False
            self.W.legend_pos = 1

        if self.share_walkers_spin and self.W.plot and self.S.plot:
            assert not (self.share_walkers_E and self.E.plot and self.W.plot)
            self.W.shared = True
            self.W.label_left = True
            self.W.legend_pos = 4
            self.S.shared = True
            self.S.label_left = False
            self.S.legend_pos = 1

        # Some details for time plots
        self.T.shared = False
        self.T.label_left = True
        self.T.legend_pos = 4

        # Some details for success/acceptance ratio plots
        self.U.shared = False
        self.U.label_left = True
        self.U.legend_pos = 4

        # Some details for growth ratio plots
        self.G.shared = False
        self.G.label_left = True
        self.G.legend_pos = 4
        
        # How many sections do we need? Count active axes, and shared axes.
        naxes = sum([1 if x.plot else 0 for x in self.axes])
        nshared = sum([1 if x.shared else 0 for x in self.axes]) / 2
        nsplits = naxes - nshared

        # New figure
        self.fig = figure()
        self.fig.subplots_adjust(hspace=0)

        # Keep track of what we have achieved
        pos_curr = nsplits
        ax_x = None

        # Create the axes objects
        for x in reversed(self.axes):
            if x.plot:
                framep = not (x.label_left and x.shared)
                x.ax = self.fig.add_subplot (nsplits, 1, pos_curr, 
                                             sharex = ax_x, frameon=framep)
                if not ax_x:
                    ax_x = x.ax
                if x.label_left or not x.shared:
                    pos_curr -= 1



    def clear_axes (self):

            for x in reversed(self.axes):
                if x.plot:
                    x.ax.cla()

            
    def axis_labels (self):
        '''Set the axis labels. This needs to be done after plotting
           otherwise we can't correctly eliminate labels which depend on the
           y-values'''

        # We only need to plot the x-coord once, so work from the bottom
        # and stop trying once we have printed it.
        xlabel_ax = None

        # S^2 is not allowed to be -ve.
        if self.S.plot:
            ax = self.S.ax
            lim = ax.get_ylim()
            ax.set_ylim (max(1.0e-6, lim[0]), lim[1])

        # Work through the axes from bottom to top, so we can find the 
        # axis to apply x-labels to first.
        for x in reversed(self.axes):
            if x.plot:
                # Set axis limits if specified
                if x.limits:
                    x.ax.set_ylim (x.limits)

                # Remove unecessary duplicated x-labels
                if xlabel_ax:
                    for label in x.ax.get_xticklabels():
                        label.set_visible(False)
                else:
                    xlabel_ax = x.ax

                # Apply log-scales if desired
                if x.log_scale:
                    x.ax.set_yscale('log')

                # Where do the labels appear? Apply them
                x.ax.set_ylabel(x.label)
                if x.label_left:
                    if x.shared:
                        x.ax.yaxis.tick_left()
                    x.ax.yaxis.set_label_position('left')
                else:
                    x.ax.yaxis.tick_right()
                    x.ax.yaxis.set_label_position('right')

                # Remove top y-axis value on all plots except the top one.
                x.ax.get_yticklabels()[-1].set_visible(False)

                # Apply legend as required
                if self.plot_legend and x.legend_pos:
                    x.ax.legend(loc=x.legend_pos)
                    #x.ax.legend(loc=x.legend_pos, ncol=2)


        # Apply the correct x-coordinate label to the correct axes.
        assert xlabel_ax is not None
        if self.x_time:
            xlabel = "Cumulative Time / s"
        elif self.x_walkers:
            xlabel = "No. of Walkers"
        elif self.x_itime:
            xlabel = "Imaginary Time / a.u."
        else:
            xlabel = "Iteration"
        xlabel_ax.set_xlabel(xlabel)


    def do_plot (self):
        '''Overall plotting wrapper'''

        self.init_axes ()
        self.plot_data ()
        self.axis_labels ()

    def do_replot (self):

        self.clear_axes ()
        self.plot_data ()
        self.axis_labels ()



    def plot_data (self):
        '''Actually plot the data'''

        # Colour management object
        col = colour_manager()

        # Plot data from the plot files.
        for i, fl in enumerate(self.in_files):

            # Give ourselves some output
            if self.verbose:
                print 'Plotting "%s" with output file "%s"' % (fl.fn, fl.ofile)

            # Do we want to use the same colour for all lines in the file?
            if self.fix_colours:
                col.freeze(False)
                col()
                col.freeze(True)

            # Store the colour used for this file, so we can repeat it.
            col.reset_point()

            with open(fl.fn, 'r') as f:

                # Read data from file, stopping at the correct point
                if self.x_itime:
                    cols = read_cols(f, last_im_time=self.last_iter)
                else:
                    cols = read_cols(f, last_iter=self.last_iter)

                # Is there an output file to process?
                tau_changes = None
                if fl.ofile is not None and self.E.total_energies:
                    ref_E, E_final, tau_changes = process_output (fl.ofile)

                    # Adjust certain columns if they are there.
                    if 'shift' in cols:
                        for i in range(len(cols['shift'])):
                            # This is a bit of a hack...
                            if abs(cols['shift'][i][-1]) < 10:
                                cols['shift'][i] = apply_ref_Es(cols['shift'][i],
                                                                 cols['iter'], ref_E)

                    if self.E.plot_exp_average and 'exp_av_projE' in cols:
                        cols['exp_av_projE'] = \
                                apply_ref_Es(cols['exp_av_projE'],
                                             cols['iter'], ref_E)

                    if self.E.plot_averages and 'av_projE' in cols:
                        cols['av_projE'] = apply_ref_Es(cols['av_projE'],
                                                        cols['iter'], ref_E)

                    if self.E.plot_averages and 'av_shift' in cols:
                        cols['av_projE'] = apply_ref_Es(cols['av_shift'],
                                                        cols['iter'], ref_E)

                # Are we using iteration, time or walker number as x-coord?
                if self.x_time:
                    x = cumulative_time(cols['iter'], cols['it_time'])
                elif self.x_walkers:
                    x = cols['parts'][0]
                elif self.x_itime:
                    x = cols['im_time']
                else:
                    x = cols['iter']

                # Do we have a legend prefix?
                if self.legend_strs:
                    leg_pre = self.legend_strs[i]
                else:
                    leg_pre = ''

                # Plot items on the energy plot.
                if self.E.plot:
                    ax = self.E.ax
                    fmt = '' if self.E.plot_lines else ','
                    for s in cols['shift']:
                        ax.plot (x, s, col()+fmt, label=' '.join([leg_pre,'Inst. Shift']))
                    if self.E.total_energies:
                        for p in cols['proje_tot']:
                            ax.plot (x, p, col()+fmt, label=' '.join([leg_pre,'Proj. E']))
                    else:
                        ax.plot (x, cols['proje_corr'], col()+fmt, label=' '.join([leg_pre,'Proj. E']))
                    if self.fit_walkers_E:
                        fit, efit, err = stretched_exp_fit(cols['parts'], cols['proje_corr'])
                        ax.plot (x, fit, col(), label=' '.join([leg_pre,'FIT']))
                        ax.axhline(efit, color=col(), label=' '.join([leg_pre,"E-fit"]))
                    if self.E.plot_averages:
                        ax.plot (x, cols['av_proje'], col()+fmt, label=' '.join([leg_pre,'Av. Proj. E']))
                        ax.plot (x, cols['av_shift'], col()+fmt, label=' '.join([leg_pre,'Av. Shift']))
                    if self.E.plot_exp_average:
                        ax.plot (x, cols['exp_av_proje'], col()+fmt, label=' '.join([leg_pre,'Exp. Average']))

                # Plot items on the walker plot
                if self.W.plot:
                    if not self.share_walkers_E:
                        col.reset()
                    ax = self.W.ax
                    for p in cols['parts']:
                        ax.plot (x, p, col(), label=' '.join([leg_pre,'No. Walkers']))
                    for r in cols['ref_parts']:
                        ax.plot (x, r, col(), label=' '.join([leg_pre,'No. at Ref']))

                    if self.W.plot_growth_cpts:
                        ax.plot (x, cols['born'], col(), label=' '.join([leg_pre,'No. Born']))
                        ax.plot (x, cols['died'], col(), label=' '.join([leg_pre,'No. Died']))
                        ax.plot (x, cols['annihil'], col(), label=' '.join([leg_pre,'No. Annihil']))

                # Plot items on the spin plot
                if self.S.plot:
                    if not self.share_walkers_spin:
                        col.reset()
                    ax = self.S.ax
                    if self.S.log_scale:
                        X, S = zip(*[(i, s) for i, s in zip(x, cols['S2']) if s > 0])
                    else:
                        X, S = x, S2
                    ax.plot (X, S, col(), label=' '.join([leg_pre,'Inst. $S^2$']))
                    
                    if self.S.log_scale:
                        X, S_init = zip(*[(i, s) for i, s in zip(x, cols['S2_init']) if s > 0])
                    else:
                        X, S_init = x, S2_init
                    ax.plot (X, S_init, col(), label=' '.join([leg_pre,'Inst. $S^2$ init']))

                # Plot time related details
                if self.T.plot:
                    ax = self.T.ax
                    col.reset()
                    if self.T.plot_tau:
                        timestep = get_timestep(cols['iter'], cols['im_time'], tau_changes)
                        ax.plot (x, timestep, col(), label=' '.join([leg_pre,'Timestep (a.u.)']))
                    if self.T.plot_iter_time:
                        ax.plot (x, cols['it_time'], col(), label=' '.join([leg_pre,'Iteration Time']))

                # Plot survival/acceptance related details
                if self.U.plot:
                    ax = self.U.ax
                    col.reset()
                    ax.plot (x, cols['accept'], col(), label=' '.join([leg_pre,'Acceptance rate']))

                if self.G.plot:
                    # Get the growth rate for the reference
                    grow_ref = [1.0] + [r2 / r1 for (r1, r2) in zip(cols['ref_parts'],
                                                                    cols['ref_parts'][1:])]
                    growth_ratio = array(cols['growth']) / array(grow_ref)

                    alpha = 0.08
                    growth_ratio_smoothed = [growth_ratio[0]]
                    elem = growth_ratio[0]
                    for g in growth_ratio[1:]:
                        elem = (1 - alpha) * elem + alpha * g
                        growth_ratio_smoothed.append(elem)

                    ax = self.G.ax
                    col.reset()
                    #ax.plot (x, growth_ratio, col(), label=' '.join([leg_pre,'Growth ratio']))
                    ax.plot (x, growth_ratio_smoothed, col(), label=' '.join([leg_pre,'Growth ratio']))

            # Jump to the highest colour used.
            col.jump_max()

        # Plot the supplied energy values
        if self.E.plot:
            for f in self.E.tgt_E:
                vals = file_E_list(f)
                for v in vals:
                    self.E.ax.axhline(v, color=col())

        # Plot a horizontal line for the growth ratios
        if self.G.plot:
            self.G.ax.axhline(1.0, color=col())






####################################
# Make the plotter global --> accessible to callback functions
####################################
plot = plotter()



def usage ():
    '''Print the usage statement'''
    print __doc__



def cumulative_time (it, itime):
    '''Calculate a cumulative time field from the iteration number and
       iteration time fields'''

    # Initialise output list
    cum_time = list()
    cum_time.append(0.0)

    # Loop over all iterations/times and sum.
    for i in range(1, len(it)):
        elem = (it[i] - it[i-1]) * itime[i]
        elem += cum_time[i-1]
        cum_time.append(elem)

    return cum_time



def get_timestep (it, im_time, tau_changes):
    '''Calculate the timestep for each iteration from the value of the
       imaginary time at each output step'''

    # Initialise output list
    tsteps = list()

    if tau_changes is not None:

        # TODO: initial tau?
        tau = 0.1
        pos = 0
        for i in it:
            if pos < len(tau_changes):
                if i >= tau_changes[pos][0]:
                    tau = tau_changes[pos][1]
                    pos += 1

            tsteps.append(tau)
    
    else:
        tsteps.append(0.0)

        # Loop over all iterations
        for i in range(1, len(it)):
            elem = (im_time[i] - im_time[i-1]) / float(it[i] - it[i-1])
            tsteps.append(elem)

    return tsteps



def process_output (fout):
    '''Process a specified output files
       --> A list of reference energies (and the iterations they apply from)
       --> The final (averaged) energy if it has been reached and printed'''

    # Return a zero reference energy, and no final energy if no output file
    if not fout:
        return [(0, 0)], None

    # Useful regular expressions
    # re_startiter = re.compile('^\s*St(ep)*\s+Shift')
    re_startiter = re.compile('^\s*Initial memory allocation sucessful...')
    re_HF_D0 = re.compile('^\s*\<D0\|H\|D0\>\s*=\s*(.+)\s*$')
    re_iter = re.compile('^\s*(\d+) ')
    re_ref_E = re.compile('^\s*Reference [Ee]nergy (now )*set to:\s*(.+)$')
    re_final_E = re.compile('^\s*Summed approx E\(Beta\)=\s*(.+)$')
    re_new_tau = re.compile('^\s*(New (tau|timestep):|Updating time-step\. New time-step =)\s+(.+)$')
    re_init_tau = re.compile('^\s*(From analysis of reference determinant '\
                             'and connections, an upper bound for the '\
                             'timestep is|Using initial time-step):\s*(.+)$')
    re_start_fciqmc = re.compile('^\s*Performing Parallel FCIQMC....')

    # Initial values and 
    HF_energies = []
    tau_changes = []
    E_final = None
    it = None
    E_HF = None
    new_tau = None
    with open(fout, 'r') as f:
        for line in f:
            if it is not None:
                # Have we restarted everything?
                m = re_start_fciqmc.match(line)
                if m:
                    HF_energies = []
                    tau_changes = []
                    E_final = None
                    it = None
                    E_HF = None
                    new_tau = None
                    continue

                # If we have hit an iteration line, update the counter. If the
                # reference energy has changed since the last occasion, we need
                # to append it to the list
                m = re_iter.match(line)
                if m:
                    # Extract iteration. Catch for output bug.
                    it = int(m.group(1))
                    it = 0 if it == 1 else it

                    if E_HF is not None:
                        HF_energies.append((it, E_HF))
                        E_HF = None

                    if new_tau is not None:
                        tau_changes.append((it, new_tau))
                        new_tau = None

                # Have we changed reference det?
                m = re_ref_E.match(line)
                if m:
                    E_HF = float(m.group(2))
                    continue

                # The final energy as calculated in the simulation
                m = re_final_E.match(line)
                if m:
                    E_final = float(m.group(1))
                    continue

                # Change of tau?
                m = re_new_tau.match(line)
                if m:
                    new_tau = float(m.group(3))
                    continue

            else:
                # Do we need to start counting iterations?
                if re_startiter.match(line):
                    it = 0
                    continue

                # Have we hit the initial setting of the HF energy?
                m = re_HF_D0.match(line)
                if m:
                    E_HF = float(m.group(1))
                    continue

                # Change of reference?
                m = re_ref_E.match(line)
                if m:
                    E_HF = float(m.group(2))
                    continue

                # Change of tau?
                m = re_init_tau.match(line)
                if m:
                    new_tau = float(m.group(2))
                    continue


    return HF_energies, E_final, tau_changes

def apply_ref_Es (E, it, ref_Es):
    '''Calculate the plottable (total) energies from the correlation energies
       in the FCIMCStats files, the reference energies from the output files
       and the list eof iteration numbers.'''

    # Main loop
    E_out = []
    ref_it, E_next = -1, ref_Es[0][1]
    ref_curr = -1
    for i, ecurr in zip(it, E):

        # Update reference energy in the list
        if (ref_it and i >= ref_it) or ref_curr == -1:
            ref_curr += 1
            ref_E = E_next
            if len(ref_Es) > ref_curr + 1:
                ref_it, E_next = ref_Es[ref_curr+1]
            else:
                ref_it = None

        # Add to the list.
        E_out.append(ecurr + ref_E)
    return E_out


def keypress_callback (event):
    '''Respond to a keypress event'''

    if event.key == 'e':
        # Turn on interactive mode, temporarily
        ion()

        plot.do_replot()
        draw()

        # Return to normal, blocking mode (allows normal exiting)
        ioff()


def resize_callback (event):
    '''Respond to a resize event'''
    
    print 'Resize event'
    print 'name', event.name
    print 'canvas', event.canvas
    print 'guiEvent', event.guiEvent
    
        
def scroll_callback (event):
    '''Respond to a scroll event'''
    
    print 'Scroll event'
    print 'name', event.name
    print 'canvas', event.canvas
    print 'guiEvent', event.guiEvent


def to_float_if_float(s):
    '''If a provided string value is a float, then return it as a float.
       Otherwise, return the string unchanged.'''
    try:
        f = float(s)
        return f
    except:
        return s


def file_E_list(f):
    '''Return a list of the float values in the (first column of) the 
       specified file, or the float value provided.'''

    if type(f) == float:
        return [f]
    else:
        ret = []
        with open(f, 'r') as fl:
            for line in fl:
                if line:
                    ret.append(float(line.split()[0]))
        return ret

def stretched_exp_fit (wlk, proj):
    '''Fit the number of walkers/projected energy to a stretched exponential
    function of the form:

    E_p(Nw) = E_0 * (1 - exp(-(Nw/c)**beta))

    This corresponds to the error taking the form:

    d(delta E)/dN = -c(delta E) / N**alpha

    Possibly try with another constant before the exp?'''

    # Fitting function, and derived error function
    fitfunc = lambda N, E, c, al: E * (1 - exp(-((N/c)**al)))

    # Initial guesses for parameters
    # E_0, c, alpha
    p0 = [-20.8, 10000.0, 0.141]

    # Fitting data. We don't want to include early iterations with far too
    # few walkers (need some stability)
    zipped_data = zip(wlk, proj)
    W, P = zip(*[z for z in zipped_data if z[0] > 10000])

    # Fit the data!
    popt, pcov = curve_fit (fitfunc, W, P, p0=p0, maxfev=500000)
#    print 'popt', popt
#    print 'pcov', pcov

    # Obtain the fitted curve to plot.
    #popt = [-20.88292143, 10000, 0.141398]
    fit_data = map(lambda N: fitfunc(N, popt[0], popt[1], popt[2]), wlk)
    
     # Output the fitting in a human readable form
    # n.b. diagonals of the covariance matrix provide the variance of the
    #      parameter estimates.
    err = sqrt(pcov[0,0]) if pcov is not inf else -1
    print 'Fit data; E_0 = %f (+-%f), c = %f, alpha = %f' % \
                                       (popt[0], err, popt[1], popt[2])
    if pcov is inf:
        print 'Poor fit obtained (covariance matrix == inf)'

    return fit_data, popt[0], err
                


# If we are running this directly, execute main.
if __name__ == "__main__":

    plot.proc_args()

    plot.do_plot()

    # Enter plotting main loop
    plot.fig.canvas.mpl_connect('key_release_event', keypress_callback)
#    plot.fig.canvas.mpl_connect('resize_event', resize_callback)
#    plot.fig.canvas.mpl_connect('scroll_event', scroll_callback)
    show()
