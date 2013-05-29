#!/usr/bin/python
'''Plot data from an FCIMCStats file (or similar, e.g. NodeFile)

Usage:
	plot.py [args] file1 [file2 ...]

Arguments:
	--average, -a       Plot average shift and projected energy values
	--reset-colours, -c When using multiple plots (split), restart the colours (within
	                    each data file) in each plot window.
	--total-energy, -e  Plot the total energy rather than the correlation energy. Uses
	                    data from the output file (if found) to do this.
	--energy-limits, -E Use the specified limits for the energy axis. Format as min,max
	                    or "min max"
	--atom-first, -F    Assume that the first data set corresponds to an atom, when we
	                    are examining dimers --> double the energy values to give a
	                    comparative dissociated value.
	--final, -f         Plot the final energy calculated on the plot
	--help, -h          Display this help
	--iterations=, -I   Which iterations to plot. This also works as times if -t is 
	                    specified. Specify as:
	                       - Max iter only
	                       - min,max
	                       - "min max"
	--linear, -l        Use a linear rather than a log scale for walker number
	--legend-pos, -L    Specify a list of positions for the legend(s)
	--output-file, -O   Use this list of output files ("A B C..." or A,B,C)
	--E-other, -o       Plot the discrete energies in file E-file as horizontal lines
	                    on the plot. Each line in the file gives a new line. The energy
						may be followed by a name to put in the legend.
	--no-projE, -P      Don't plot the projected energy
	--eps, -p           Produce an eps output file, rather than outputting to screen
	--print-ref-E, -r   Print the reference energies being used.
	--repeat, -R        Repeat the plot every (arg, default=10) seconds. 
	--split, -s         Plot walkers and energies in separate axes
	--no-shift, -S      Don't plot the shift
	--time, -t          Plot using cumulative time rather than iteration on x-axis
	--no-walkers, -W    Don't plot walker numbers
'''

from pylab import *
from numpy import *
from matplotlib import rc, rcParams
from os import path
import matplotlib.axes
import getopt
import re
import time

class plotter:
	def __init__ (self):
		'''Default arguments'''
		self.x_cum_time = False
		self.disp_avg = False
		self.linear_walker = False
		self.process_out = False
		self.plot_E_final = False
		self.plot_walkers = True
		self.atom_first = False
		self.split_plot = False
		self.reset_colours = False
		self.print_ref_E = False
		self.total_E = True
		self.files = ['FCIMCStats']
		self.fig = None
		self.ax1 = None
		self.ax2 = None
		self.bounds = None
		self.output_files = None
		self.E_limits = None
		self.num_lines = 0
		self.E_other = None
		self.repeat = None
		self.legend_pos = None
		self.plot_shift = True
		self.plot_projE = True
		self.eps_out = None

	def proc_args (self, args):
		'''Read in the command line options. For default arguments see __init__'''

		# Use getopt to process the arguments
		try:
			opts, self.files = getopt.getopt(args, "htalE:efFWsI:O:cro:R:L:SPp:", ["help", "time", "average", "linear", "total-energy", "energy-limits=", "atom-first", "final", "no-walkers", "split", "iterations=", "output-file=", "reset-colours", "reset-colors", "print-ref-E", "E-other=", "repeat", "legend-pos", "no-shift", "no-projE", "eps="])
		except getopt.GetoptError, err:
			print str(err)
			usage()
			sys.exit(2)

		if not self.files:
			self.files.append("FCIMCStats")

		for o, a in opts:
			if o in ("-h", "--help"):
				usage()
				sys.exit()
			elif o in ("-t", "--time"):
				self.x_cum_time = True
			elif o in ("-a", "--average"):
				self.disp_avg = True
			elif o in ("-l", "--linear"):
				self.linear_walker = True
			elif o in ("-e", "--total-energy"):
				self.total_E = True
				self.process_out = True
			elif o in ("-E", "--energy-limits"):
				limits = a.replace(',', ' ').split()
				if len(limits) == 2:
					self.E_limits = (float(limits[0]), float(limits[1]))
				else:
					print 'Invalid limits specified: %s' % a
					usage()
					sys.exit(2)
			elif o in ("-F", "--atom-first"):
				self.atom_first = True
			elif o in ("-W", "--no-walkers"):
				self.plot_walkers = False
			elif o in ("-s", "--split"):
				self.split_plot = True
			elif o in ("-I", "--iterations"):
				limits = a.replace(',', ' ').split()
				if len(limits) == 2:
					self.bounds = (int(limits[0]), int(limits[1]))
				elif len(limits) == 1:
					self.bounds = (0, int(limits[0]))
				else:
					print 'Invalid limits specified: %s' % a
					usage()
					sys.exit(2)
			elif o in ("-f", "--final"):
				self.plot_E_final = True
				self.total_E = True
				self.process_out = True
				print 'Printing final E on plot. Requires printing total energies, setting --total-energy (-E)'
			elif o in ("-O", "--output-file"):
				files = a.replace(',', ' ').split()
				if len(files) == len(self.files):
					self.output_files = {}
					for i in range(len(files)):
						self.output_files[self.files[i]] = files[i]
				else:
					print 'Insufficient output files specified. Need to specify %d files', len(self.files)
					print 'Input files: ', self.files
					usage()
					sys.exit(2)
			elif o in ("-o", "--E-other"):
				self.E_other = a
			elif o in ("-c", "--reset-colours", "--reset-colors"):
				self.reset_colours = True
			elif o in ("-r", "--print-ref-E"):
				self.print_ref_E = True
			elif o in ("-R", "--repeat"):
				if len(a) > 0:
					self.repeat = int(a)
				else:
					self.repeat = 10
			elif o in ("-L", "--legend-pos"):
				positions = a.replace(',', ' ').split()
				if len(positions) != 2:
					print 'The positions of 2 legends must be supplied.'
					usage()
					sys.exit(2)
				else:
					self.legend_pos = []
					for i in range(len(positions)):
						self.legend_pos.append(int(positions[i]))
			elif o in ("-S", "--no-shift"):
				self.plot_shift = False
			elif o in ("-P", "--no-projE"):
				self.plot_projE = False
			elif o in ("-p", "--eps"):
				self.eps_out = a
			else:
				assert False, "Unhandled Option"	

	def get_colour (self):
		'''Get the colour to use for a given line'''

		# Extract the default set of colours
		colors = matplotlib.axes._process_plot_var_args.defaultColors
		
		ind = self.num_lines % len(colors)
		self.num_lines += 1
		return colors[ind]

	@staticmethod
	def update_energies (E, it, HF_E):
		'''Calculate the total energies from the correlation energies provided using
		   the list of HF energies by iteration provided'''

		E_out = []
		if HF_E:
			E_curr = 0
			HF_pos = 0
			for i in range(len(it)):
				if HF_pos is not None and it[i] >= HF_E[HF_pos][0]:
					E_curr = HF_E[HF_pos][1]
					#print 'updated E_curr', it[i], E_curr
					HF_pos += 1
					if HF_pos >= len(HF_E):
						HF_pos = None

				E_out.append(E[i] + E_curr)

		return array(E_out)

	@staticmethod
	def process_output (filename, find_out, print_ref_E):
		'''Process the OUT, OUTPUT, OUTPUT.FCIMC etc. file to get additional details to plot.
		   If required, we need to work out what the filename is from the FCIMCStats filename...'''

		fout = None
		if find_out:
			if path.basename(filename)[0:10] == 'FCIMCStats':
				out_names = ["OUT", "OUTPUT", "OUT.FCIMC", "OUTPUT.FCIMC", "OUT.FCIMC.init", "OUT.FCIMC.csf", "OUTPUT.FCIMC.init", "OUTPUT.FCIMC.csf"]
				if len(path.basename(filename)) > 10:
					for i in range(len(out_names)):
						out_names[i] = out_names[i] + path.basename(filename)[10:]
				dirname = path.dirname(filename)
				if dirname != '':
					dirname = dirname + '/'

				for tail in out_names:
					if path.isfile(dirname + tail):
						fout = dirname + tail
						break
				
		else:
			fout = filename

		if not fout:
			print 'Cannot find associated output file for stats file: %s' % filename


		
		# Regular expressions for analysing output file
		re_startiter = re.compile('^\s*St(ep)*\s+Shift')
		re_HF_D0 = re.compile('^\s*\<D0\|H\|D0\>\s*=\s*(.+)\s*$')
		re_iter = re.compile('^\s*(\d+)')
		re_ref_E = re.compile('^\s*Reference [Ee]nergy (now )*set to:\s*(.+)$')
		re_final_E = re.compile('^\s*Summed approx E\(Beta\)=\s*(.+)$')

		# Output lists
		HF_energies = []
		E_final = 0.0

		if fout:
			print 'Using output filename: %s' % fout

			with open(fout, 'r') as f:
				it = None     # What iteration are we on
				E_HF = None   # What is the current HF energy
				for line in f:
					if it is not None:
						# If we have hit an iteration line, then update counter
						# and output as necessary.
						m = re_iter.match(line)
						if m:
							it = int(m.group(1))
							# For some reason, the first iteration is output as 1 not 0 in
							# the output (not the FCIMCStats) file if we have restarted a
							# calculation.
							if it == 1:
								it = 0
							if E_HF is not None:
								if print_ref_E:
									print 'Changed reference energy: %f' % E_HF
								HF_energies.append((it, E_HF))
								E_HF = None

						# Have we changed reference det. Need to deal with a restart later.
						m = re_ref_E.match(line)
						if m:
							E_HF = float(m.group(2))
							continue

						# The final energy as calculated in the simulation
						m = re_final_E.match(line)
						if m:
							E_final = float(m.group(1))
							continue
					else:
						# Do we need to start counting iterations
						if re_startiter.match(line):
							it = 0
							continue

						# Have we hit the initial setting of the HF energy?
						m = re_HF_D0.match(line)
						if m:
							E_HF = float(m.group(1))
							continue

						# Or with a change of reference?
						m = re_ref_E.match(line)
						if m:
							E_HF = float(m.group(2))
							continue

		return HF_energies, E_final

	def do_plot (self):
		'''Actually draw the plot!'''

		# Enable TeX output and set default line formatting
		rc('text', usetex=True)
		if self.eps_out:
			rc('lines', linewidth=2)
			rc('axes', labelsize=10)
			rc('text', fontsize=10)
			rc('legend', fontsize=10)
			rc('xtick', labelsize=10)
			rc('ytick', labelsize=10)
		else:
			rc('lines', linewidth=1)
			rc('legend', fontsize=rcParams['axes.labelsize'])

		# Generate the axes
		new_plot = False
		if not self.fig:
			self.fig = figure()
			if self.split_plot:
				self.ax1 = self.fig.add_subplot('212')
				self.ax2 = self.fig.add_subplot('211', sharex=self.ax1)
				self.fig.subplots_adjust(hspace=0)
			else:
				self.ax1 = self.fig.add_subplot('111')
				self.ax2 = self.fig.add_subplot('111', sharex=self.ax1, frameon=False)
			new_plot = True

		# Create aliases to the stored objects for ease.
		ax1 = self.ax1
		ax2 = self.ax2

		# Save and restore the coordinates as needed.
		# TODO: Need to shift coordinates from the window if adjusted...
		savecoords=()
		if not new_plot:
			savecoords = (ax1.get_xlim(), ax1.get_ylim(), ax2.get_xlim(), ax2.get_ylim())
			ax1.cla()
			ax2.cla()

			ax1.set_xlim(savecoords[0])
			ax1.set_ylim(savecoords[1])
			#ax2.set_xlim(savecoords[2]) # Not needed as axes share x-limits
			ax2.set_ylim(savecoords[3])

		# Reset plot count to get colours right
		self.num_lines = 0

		# Iterate through all the specified FCIMCStats files and display
		done_first = False
		for fl in self.files:
			file_start_col = self.num_lines
			file_top_col = self.num_lines
			with open(fl, 'r') as f:
				it, sft, wlk, avProj, avSft, proj, atRef, itime = \
						loadtxt (f, usecols=(0,1,4,8,9,10,11,15), unpack=True)

				# Do we want to append details from the output file?
				E_final = 0.0
				if self.process_out:
					lookup = (self.output_files[fl],False) if self.output_files else (fl,True)
					HF_energies, E_final = plotter.process_output (lookup[0], lookup[1], self.print_ref_E)

					if self.total_E and HF_energies:
						sft = plotter.update_energies(sft, it, HF_energies)
						proj = plotter.update_energies(proj, it, HF_energies)
						avProj = plotter.update_energies(avProj, it, HF_energies)
						avSft = plotter.update_energies(avSft, it, HF_energies)

				# Are we considering an atom first, amongst dimers?
				if self.atom_first and not done_first:
					done_first = True
					sft = 2 * sft
					proj = 2 * proj
					avProj = 2 * avProj
					avSft = 2 * avSft

				# Are we using cumulative time as one axis?
				if self.x_cum_time:
					x = cumulative_time(it, itime)
				else:
					x = it

				# Plot the walker counts
				if self.plot_walkers:
					ax1.plot(x, wlk, self.get_colour(), label='No. Walkers')
					ax1.plot(x, atRef, self.get_colour(), label='No. at Ref')

				# Do we want the colours of each of the files to match between the
				# walker and energy plots?
				if self.reset_colours:
					file_top_col = self.num_lines
					self.num_lines = file_start_col

				# Plot the energies
				if self.plot_projE:
					ax2.plot(x, proj, self.get_colour(), label='Proj. E')
				if self.plot_shift:
					ax2.plot(x, sft, self.get_colour(), label='Inst. Shift')

				# Plot the final energy
				if self.plot_E_final:
					ax2.plot([x[0],x[len(x)-1]], [E_final, E_final], self.get_colour(), label='Final E')

				# Plot the averages
				if self.disp_avg and self.plot_projE:
					ax2.plot(x, avProj, self.get_colour(), label='Av. Proj. E')
				if self.disp_avg and self.plot_shift:
					ax2.plot(x, avSft, self.get_colour(), label='Av. Shift')

				self.num_lines = max(self.num_lines, file_top_col)

		# Plot extra energy values as requested
		if self.E_other:
			xlim = ax1.get_xlim()
			with open(self.E_other, 'r') as f:
				for line in f:
					bits = line.split()
					ecurr = float(bits[0])
					bits = "".join(bits[1:])

					ax2.plot([xlim[0],xlim[1]], [ecurr,ecurr], self.get_colour(), label=bits)


		# Axis labels etc.
		ax1.set_ylabel("Number of Walkers")
		if not self.linear_walker:
			ax1.set_yscale('log')

		ax2.set_ylabel("Energy / $E_h$")

		if self.split_plot:
			# Don't display axis labels on the top plot
			for label in ax2.get_xticklabels():
				label.set_visible(False)

			# Remove the last y-axis label on the lower plot
			ax1.get_yticklabels()[len(ax1.get_yticklabels())-1].set_visible(False)
		else:
			ax1.yaxis.tick_left()
			ax2.yaxis.tick_right()
			ax2.yaxis.set_label_position("right")

		if self.x_cum_time:
			ax1.set_xlabel("Cumulative Time / s")
		else:
			ax1.set_xlabel("Iteration")

		# If bounds have been set, apply them here
		if self.bounds:
			ax1.set_xlim(self.bounds)
		if self.E_limits:
			ax2.set_ylim(self.E_limits)

		# Put the legends on, in the right places, and display
		if self.legend_pos:
			ax1.legend(loc=self.legend_pos[0], ncol=2)
			ax2.legend(loc=self.legend_pos[1], ncol=2)
		elif self.split_plot:
			ax1.legend(loc=4, ncol=2)
			ax2.legend(loc=1, ncol=2)
		else:
			ax1.legend(loc=2)
			ax2.legend(loc=1)

		# Connect callback(s)
		self.fig.canvas.mpl_connect('key_release_event', keypress_callback)

		if self.eps_out:
			savefig(self.eps_out)
		else:
			show()

###########################################
# Declare the plot object as global, so the callbacks can access it!
###########################################
plot = plotter()


def usage ():
	'''Print the usage statement'''
	print __doc__

def cumulative_time (it, itime):
	'''Calculate a cumulative time field from the iteration number and iteration time fields'''

	cum_time = list()
	cum_time.append(0.0)
	for i in range(1, len(it)):
		elem = (it[i] - it[i-1]) * itime[i]
		elem += cum_time[i-1]
		cum_time.append(elem)
	return cum_time

def keypress_callback (event):
	'''Respond to a keypress event'''

	if event.key == 'r':
		plot.do_plot()



# If we are running this directly, execute main.
if __name__ == "__main__":

	plot.proc_args(sys.argv[1:])

	plot.do_plot()

	if plot.repeat:
		while(True):
			time.sleep(5)
			plot.do_plot()

