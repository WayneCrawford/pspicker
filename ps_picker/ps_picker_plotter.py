import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


class PSPickerPlotter():
    """
    Class to plot figures representing the work in PSPicker

    There are three types of plots:
    Figure 1: Global plot showing all stations and the entire time
    Figure 2: (Global) Pick plot showing all stations but zoomed in on
              pick window
    Figures 3+: Station plots, each one divided in four:
        1) Picks
        2) Data
        3) Kurtosis
        4) Energy and signal-to-noise level
    """
    def __init__(self, show_plots=True):
        """
        :param show_plots: make the plots or not?
        :kind show_plots: bool
        """
        self.show_plots = show_plots
        self.fig1 = None  # Global plot figure
        self.ax1 = None   # Global plot axis
        self.fig2 = None  # Global pick plot figure
        self.ax2 = None   # Global pick plot axis
        self.fig3 = None  # Station plot figure
        self.axs3 = None  # Station plot axes (list)

    def global_window_setup(self):
        """
        Create Figure 1 (global picking window) (figure 1)
        """
        if not self.show_plots:
            return
        self.fig1, self.ax1 = plt.subplots(num='PSPicker: Overview')

    def global_window_onetrace(self, trace, iTrace, extrema_s):
        """
        :param trace: trace to plot
        :param iTrace: y-index of the trace (used to separate traces)
        :param extrema: offset time of all extrema for this trace
        """
        # Pick_Function.m: ~160?
        if not self.show_plots:
            return
        ax = self.ax1
        t_begin = trace.stats.starttime.matplotlib_date
        norm_trace = trace.copy().detrend().normalize()
        ax.plot_date(norm_trace.times(type="matplotlib"),
                     (0.5 * norm_trace.data) + iTrace,
                     fmt='-',
                     label=trace.stats.station)
        distri_times = [t_begin+s/86400 for s in extrema_s]
        ax.vlines(distri_times, iTrace + 0.5, iTrace - 0.5, color='k')
        plt.draw()
        plt.show(block=False)

    def global_window_timebounds(self, rect_start_time, rect_end_time,
                                 stations):
        """
        Plot a rectangle representing the global time bounds and set axis limits
        
        :param rect_start_time: global window start time
        :param rect_end_time: global window end time
        :param stations: list of stations in same order as iTrace above
        """
        if not self.show_plots:
            return
        # Pick_Function.m:231
        ax = self.ax1
        w1_s = rect_start_time.matplotlib_date
        w1_e = rect_end_time.matplotlib_date
        N = len(stations)
        patches = [Rectangle((w1_s, -0.5), width=w1_e-w1_s, height=N + 1,
                   color='lightblue', alpha=0.1, edgecolor=None)]
        ax.add_collection(PatchCollection(patches))
        ax.set_title('Phase 1: Calculate global window')
        ax.set_xlabel(rect_start_time.strftime('%Y-%m-%d'))
        ax.set_ylim(-0.5, N + 0.5)
        ax.set_yticks(np.arange(0, N))
        ax.set_yticklabels(stations)
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)

    def pick_window_start(self, starttime, endtime, stations):
        """
        :param starttime: plot start time
        :param endtime: plot end time
        :param stations: sorted list of stations
        """
        if not self.show_plots:
            return
        # Pick_Function.m:250
        fig, ax = plt.subplots(num='PSPicker: all picks')
        N = len(stations)
        fig.clf()
        ax.set_xlim(starttime.matplotlib_date, endtime.matplotlib_date)
        ax.set_ylim(-0.5, N - 0.5)
        ax.set_yticks(np.arange(0, N))
        ax.set_yticklabels(stations)
        ax.set_title('Phase 3: Select final picks based on clustering')
        self.fig2 = fig
        self.ax2 = ax
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)
        # datetick('keeplimits')
        # hold('on')

    def station_window_start(self, trace, iter):
        """
        Plot pick phases on Figure 3+

        :param trace: trace to plot
        :param iter: station iteration number
        """
        if not self.show_plots:
            return
        # Pick_Function.m:334
        #fig, axs = plt.subplots(4, 1, num=2+iter)
        fig, axs = plt.subplots(4, 1, num=trace.stats.station)
        fig.clf()
        axs[0].plot(trace.times(type="matplotlib"), trace.data, 'k')
        axs[0].set_title('Phase 2: Pick phases (station {})'
                         .format(trace.stats.station))
        self.fig3 = fig
        self.axs3 = axs
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)

    def station_window_add_snr_nrg(self, picker, snr, energy, all_mean_M,
                                   kurto_modif):
        """
        Figure 3+: individual station plots of:
            a) Pick windows
            b) Data
            c) kurtosis
            d) Energy and signal-to-noise
        """
        if not self.show_plots:
            return
        # Pick_Function.m:443
        axs = self.axs3
        fig = self.fig3
        sr = snr.stats.sampling_rate
        imin, imax = np.nonzero(np.isfinite(all_mean_M))[0][[0, -1]]
        xmin, xmax = snr.times(type="matplotlib")[[imin, imax]]
        # xmin = find(not np.isnan(all_mean_M), 1)[0] / sr
        # xmax = find(not np.isnan(all_mean_M), 1)[-1] / sr
        w1_s = picker.run.global_first_time.matplotlib_date
        w1_e = picker.run.global_last_time.matplotlib_dates
        # plt.figure(2 + iter)
        # plt.subplot(4, 1, 1)
        # Pick windows
        patches = [Rectangle((xmin, picker.loop.dmin),
                             width=xmax-xmin,
                             height=picker.loop.dmax - picker.loop.dmin,
                             color=[.9, .9, .9], edgecolor=None)]
        axs[0].add_collection(PatchCollection(patches))
        # plt.patch(concat([xmin, xmax, xmax, xmin, xmin]),
        #       concat([self.loop.dmin, self.loop.dmin, self.loop.dmax,
        #               self.loop.dmax, self.loop.dmin]),
        #       dot(-0.1, ones(1, 5)), concat([0.9, 0.9, 0.9]),
        #       'EdgeColor', 'none')
        axs[0].set_xlim(w1_s, w1_e)
        axs[0].set_ylim(picker.loop.dmin, picker.loop.dmax)

        # Data
        tr = picker.loop.datP[0].copy()
        tr.data[np.isnan(all_mean_M)] = np.nan
        axs[1].plot(tr.times(type='matplotlib'), tr.data, 'k')
        axs[1].set_xticklabels([])

        # Kurtosis
        axs[2].plot(all_mean_M.times(type='matplotlib'), all_mean_M.data, 'b')
        axs[2].set_ylabel('Kurtosis', color='b')
        axs[2].tick_params(axis='y', labelcolor='b')
        axs[2].set_xticklabels([])
        ax2b = axs[2].twinx()  # Create a twin axis with a different y-scale
        ax2b.plot(kurto_modif[0].times(type='matplotlib'), kurto_modif[0].data, 'r--')
        ax2b.plot(kurto_modif[-1].times(type='matplotlib'), kurto_modif[-1].data, 'r')
        ax2b.set_ylabel('Kurtosis Extrema', color='r')
        ax2b.tick_params(axis='y', labelcolor='r')
        fig.tight_layout()   # (to avoid clipping right label)
        self.kmin, self.kmax = np.min(all_mean_M), np.max(all_mean_M)
        # axa, h1, h2 = plotyy(xax, all_mean_M, xax, kurto_modif(arange(),
        #                      concat([1, end()])), nargout=3)
        # plt.hold(axa(1), 'on')
        # plt.set(h1(1), 'Color', 'b')
        # plt.set(axa(1), 'YColor', 'b')
        # plt.set(h2(1), 'Color', 'r', 'LineStyle', '--')
        # plt.set(h2(2), 'Color', 'r')
        # plt.set(axa(2), 'YColor', 'r')
        # kMinMax = concat([min(all_mean_M), max(all_mean_M)])
        # plt.ylabel(axa(1), 'Kurtosis')
        # plt.ylabel(axa(2), 'Kurtosis Extrem')
        # plt.set(gca, 'XTickLabel', [])

        # Energy and SNR
        nrgplot = energy.copy()
        nrgplot.data = 20 * np.log10(nrgplot.data)
        nrgplot.data[np.isnan(all_mean_M)] = np.nan
        snrplot = snr.copy()
        snrplot.data[np.isnan(all_mean_M)] = np.nan
        snrthresh = snrplot.copy()
        snrthresh[np.isfinite(snrthresh)] = picker.param.SNR_thres[0]
        axs[3].plot_date(snr.times(type='matplotlib'), snrplot.data, 'b')
        axs[3].plot_date(snr.times(type='matplotlib'), snrthresh.data, 'b--')
        axs[3].set_ylabel('SNR (dB)', color='b')
        axs[3].tick_params(axis='y', labelcolor='b')
        ax3b = axs[3].twinx()  # Create a twin axis with a different y-scale
        ax3b.plot_date(nrgplot.times(type='matplotlib'), nrgplot.data, 'r')
        ax3b.set_ylabel('Energy (dB)', color='r')
        ax3b.tick_params(axis='y', labelcolor='r')
        # axb, h1, h2 = plt.plotyy(
        #     xax, concat([ravel(snrplot), ravel(snrthresh)]), xax, nrgplot)
        # plt.set(h1(1), 'Color', 'b')
        # plt.set(h1(2), 'Color', 'b', 'LineStyle', '--')
        # plt.set(axb(1), 'YColor', 'b')
        # plt.set(h2(1), 'Color', 'r')
        # plt.set(axb(2), 'YColor', 'r')
        # plt.ylim(axb(1), [0, max(snrplot)])
        # plt.ylabel(axb(1), 'SNR (dB)')
        # plt.ylabel(axb(2), 'Energy (dB)')
        axs[3].xlabel('Time')
        plt.draw()
        plt.show(block=False)

    def station_window_add_extrema(self, picker, extrema):
        """
        Add extrema to figure 3+ axis[2]
        """
        if not self.show_plots:
            return
        # Pick_Function.m:502
        for extremum in extrema:
            tval = picker.loop.index_to_time(extremum['i']).matplotlib_date
            ax = self.axs[0]
            #  plt.subplot(4, 1, 1)
            ax.vlines(tval, picker.loop.dmin, picker.loop.dmax,
                      color=[0.8, 0.8, 0.8])
            ax = self.axs[1]
            # plt.subplot(4, 1, 2)
            ax.vlines(tval, picker.loop.dmin, picker.loop.dmax,
                      color=[0.8, 0.8, 0.8])
            # h = plt.plot(axa(1), [extremum['i'], extremum['i']])
            # set(h, 'Color', [0.8, 0.8, 0.8])

    def station_window_Ptrace(self, picker, iter):
        """
        Plot P traces for one station

        iter is used to separate the traces on the y axis
        """
        if not self.show_plots:
            return
        ax2 = self.ax2
        # Pick_Function.m:601
        # for j in length(datP(1,arange())).reshape(-1):
        norm_trace = picker.loop.datP[0].copy()
        norm_trace.normalize()
        ax2.plot_date(norm_trace.times(type="matplotlib"),
                      norm_trace.data + iter, 'k')
        ax2.set_xlim(picker.run.global_first_time.matplotlib_date,
                     picker.run.global_last_time.matplotlib_date)

    def plot_onsets(self, picker, iter, onset_P, onset_S):
        """
        Plot onsets for one station

        onsets are indices
        iter is used to separate the traces on the y axis
        """
        if not self.show_plots:
            return
        axs = self.axs3
        ax2 = self.ax2
        # Pick_Function.m:576

        # Pick Window: add picks to subplots 1, 2 and 3
        if onset_P is not None:
            p_time = picker.loop.index_to_time(onset_P).matplotlib_date
            axs[0].vlines(p_time, picker.loop.dmin, picker.loop.dmax, 'b')
            axs[1].vlines(p_time, picker.loop.dmin, picker.loop.dmax, 'b')
            axs[2].vlines(p_time, self.kmin, self.kmax, 'b')
            ax2.vlines(p_time, iter-0.5, iter+0.5, 'b:')
        if onset_S is not None:
            s_time = picker.loop.index_to_time(onset_S).matplotlib_date
            axs[0].plot(s_time, picker.loop.dmin, picker.loop.dmax, 'r')
            axs[1].plot(s_time, picker.loop.dmin, picker.loop.dmax, 'r')
            axs[2].plot(s_time, self.kmin, self.kmax, 'r')
            ax2.vlines(s_time, iter-0.5, iter+0.5, 'r:')

    def pick_window_add_picks(self, picker, picks):
        """
        Plot initial and final P and S picks for all stations

        station waveforms are plotted at positions corresponding to the order
        of stations in self.run.stations
        """
        if not self.show_plots:
            return
        # Pick_Function.m:811
        ax = self.ax2
        delY = 1
        for phase, color, rect_color in zip(['P', 'S'],
                                            ['b', 'r'],
                                            [[.8, .8, 1.], [1., .8, .8]]):
            picks = [p for p in picks if p.phase_hint[0] == phase]
            if len(picks) > 0:
                window_start = picks[0].time.matplotlib_date
                for pick in picks:
                    j = picker.run.stations.index(pick.station)
                    if pick.time < window_start:
                        window_start = pick.time
                    ax.vlines(pick.time.matplotlib_date,
                              j - 0.5*delY, j + 0.5*delY, colors=color)
                p = [Rectangle((window_start, 0.5),
                               width=picker.param.assoc_cluster_windows[phase]/86400,
                               height=picker.run.n_stations,
                               color=rect_color, edgecolor=None, alpha=0.5)]
                ax.add_collection(PatchCollection(p))
        ax.set_xlabel(picker.loop.t_begin.strftime('%Y-%m-%d'))
        plt.draw()
        plt.show(block=False)

# def move_figure(f, x, y):
#     """Move figure's upper left corner to pixel (x, y)"""
#     backend = matplotlib.get_backend()
#     if backend == 'TkAgg':
#         f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
#     elif backend == 'WXAgg':
#         f.canvas.manager.window.SetPosition((x, y))
#     else:
#         # This works for QT and GTK
#         # You can also use window.setGeometry
#         f.canvas.manager.window.move(x, y)
# 
# def move_figure(position="top-right"):
#     '''
#     Move and resize a window to a set of standard positions on the screen.
#     Possible positions are:
#     top, bottom, left, right, top-left, top-right, bottom-left, bottom-right
#     '''
# 
#     mgr = plt.get_current_fig_manager()
#     mgr.full_screen_toggle()  # primitive but works to get screen size
#     py = mgr.canvas.height()
#     px = mgr.canvas.width()
# 
#     d = 10  # width of the window border in pixels
#     if position == "top":
#         # x-top-left-corner, y-top-left-corner, x-width, y-width (in pixels)
#         mgr.window.setGeometry(d, 4*d, px - 2*d, py/2 - 4*d)
#     elif position == "bottom":
#         mgr.window.setGeometry(d, py/2 + 5*d, px - 2*d, py/2 - 4*d)
#     elif position == "left":
#         mgr.window.setGeometry(d, 4*d, px/2 - 2*d, py - 4*d)
#     elif position == "right":
#         mgr.window.setGeometry(px/2 + d, 4*d, px/2 - 2*d, py - 4*d)
#     elif position == "top-left":
#         mgr.window.setGeometry(d, 4*d, px/2 - 2*d, py/2 - 4*d)
#     elif position == "top-right":
#         mgr.window.setGeometry(px/2 + d, 4*d, px/2 - 2*d, py/2 - 4*d)
#     elif position == "bottom-left":
#         mgr.window.setGeometry(d, py/2 + 5*d, px/2 - 2*d, py/2 - 4*d)
#     elif position == "bottom-right":
#         mgr.window.setGeometry(px/2 + d, py/2 + 5*d, px/2 - 2*d, py/2 - 4*d)
# 
# import pyfig as fig
# for ix in range(6): f = p.figure(ix)
# fig.stack('all')
# fig.stack(1,2)
# fig.hide(1)
# fig.restore(1)
# fig.tile()
# fig.pile()
# 
# if __name__ == '__main__':
# 
#     # Usage example for move_figure()
# 
#     plt.figure(1)
#     plt.plot([0, 1])
#     move_figure("top-right")
# 
#     plt.figure(2)
#     plt.plot([0, 3])
#     move_figure("bottom-right")