import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


class Plotter():
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
    def __init__(self, plot_global=True, plot_stations=True):
        """
        :param plot_global: show global selection and picks plots
        :kind plot_global: bool
        :param plot_stations: show individual station plots
        :kind plot_stations: bool
        """
        self.plot_stations = plot_stations
        self.plot_global = plot_global
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
        if not self.plot_global:
            return
        self.fig1, self.ax1 = plt.subplots(num='PSPicker: Overview')

    def global_window_onetrace(self, trace, iTrace, extrema_s):
        """
        :param trace: trace to plot
        :param iTrace: y-index of the trace (used to separate traces)
        :param extrema_s: offset time of all extrema for this trace
        """
        # Pick_Function.m: ~160?
        if not self.plot_global:
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
        Plot the global time bounds and set axis limits

        :param rect_start_time: global window start time
        :param rect_end_time: global window end time
        :param stations: list of stations in same order as iTrace above
        """
        if not self.plot_global:
            return
        # Pick_Function.m:231
        ax = self.ax1
        N = len(stations)
        patches = [Rectangle((rect_start_time.matplotlib_date, -0.5),
                             width=rect_end_time.matplotlib_date
                             -rect_start_time.matplotlib_date,
                             height=N + 1, color='lightblue', alpha=0.1,
                             edgecolor=None)]
        ax.add_collection(PatchCollection(patches))
        ax.set_title('Phase 1: Calculate global window')
        ax.set_xlabel(rect_start_time.strftime('%Y-%m-%d'))
        ax.set_ylim(-0.5, N + 0.5)
        ax.set_yticks(np.arange(0, N))
        ax.set_yticklabels(stations)
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)

    def pick_window_setup(self, starttime, endtime, stations):
        """
        :param starttime: plot start time
        :param endtime: plot end time
        :param stations: sorted list of stations
        """
        if not self.plot_global:
            return
        # Pick_Function.m:250
        fig, ax = plt.subplots(num='PSPicker: all picks')
        N = len(stations)
        # fig.clf()
        ax.set_title('Phase 3: Select final picks based on clustering')
        # ax.set_xlim(starttime.matplotlib_date, endtime.matplotlib_date)
        ax.set_xlim((starttime-30).datetime, endtime.datetime)
        ax.set_ylim(-0.5, N - 0.5)
        ax.set_yticks(np.arange(0, N))
        ax.set_yticklabels(stations)
        self.fig2 = fig
        self.ax2 = ax
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)

    def pick_window_add_picks(self, picks, stations, t_begin,
                              P_cluster_width, S_cluster_width):
        """
        Plot initial and final P and S picks for all stations

        :picks: all picks
        :stations: ordered list of stations.  The y position depends on the
                   position of the pick's station in this list
        :param t_begin: begin time for all plots
        :param P_cluster_width: width of P cluster window
        :param S_cluster_width: with of S cluster window
        """
        if not self.plot_global:
            return
        # Pick_Function.m:811
        self._plot_picks_phase(picks, 'P', 'b', [.8, .8, 1.],
                               P_cluster_width, stations)
        self._plot_picks_phase(picks, 'S', 'r', [1., .8, .8],
                               S_cluster_width, stations)
        self.ax2.set_xlabel(t_begin.strftime('%Y-%m-%d'))
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)

    def _plot_picks_phase(self, picks, phase, color, rect_color,
                          width, stations):
        """
        :param phase: "P" or "S"
        :param color: color for picks
        :param rec_color: color for pick window
        :parm width: width of pick window in seconds
        :param stations: ordered list of stations
        """
        ax = self.ax2
        picks = [p for p in picks if p.phase_hint[0] == phase]
        if len(picks) > 0:
            window_start = picks[0].time
            for pick in picks:
                if pick.time < window_start:
                    window_start = pick.time
                j = stations.index(pick.waveform_id.station_code)
                # ax.vlines(pick.time.matplotlib_date,
                ax.vlines(pick.time.datetime, j - 0.5, j + 0.5, colors=color)
                print('plotting {} pick at {}, {} to {}'.format(
                      color, pick.time, j-0.5, j+0.5))
            p = [Rectangle((window_start.matplotlib_date, 0.5),
                           width=width / 86400, height=len(stations),
                           color=rect_color, edgecolor=None, alpha=0.5)]
            print('plotting rectangle from {},{}, width={}, height={}'
                  .format(window_start, 0.5, width / 86400, len(stations)))
            ax.add_collection(PatchCollection(p))

    def station_window_setup(self, trace, iter):
        """
        Plot pick phases on Figure 3+

        :param trace: trace to plot
        :param iter: station iteration number
        """
        if not self.plot_stations:
            return
        # Pick_Function.m:334
        # fig, axs = plt.subplots(4, 1, num=2+iter)
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

    def station_window_add_snr_nrg(self, starttime, endtime, snr_threshold,
                                   trace, snr, energy, all_mean_M,
                                   kurto_modif):
        """
        Figure 3+: individual station plots of:
            a) Pick windows
            b) Data
            c) kurtosis
            d) Energy and signal-to-noise

        :param starttime: starttime for all plots
        :param endtime: endtime for all plots
        :param snr_threshold: minimum acceptable SNR
        :param trace: data trace to plot
        :param snr: signal-to-noise level trace
        :param energy: energy trace
        :param all_mean_M: mean Kurtosis trace?
        :kurto_modif: modified kurtosis trace?
        """
        if not self.plot_stations:
            return
        # Pick_Function.m:443
        axs = self.axs3
        fig = self.fig3
        imin, imax = np.nonzero(np.isfinite(all_mean_M))[0][[0, -1]]
        tmin, tmax = all_mean_M.times(type="matplotlib")[[imin, imax]]
        dmin, dmax = np.nanmin(trace.data), np.nanmax(trace.data)
        # Pick windows
        patches = [Rectangle((tmin, dmin),  width=tmax - tmin,
                             height=dmax - dmin,
                             color=[.9, .9, .9], edgecolor=None)]
        axs[0].add_collection(PatchCollection(patches))
        axs[0].set_xlim(starttime.matplotlib_date, endtime.matplotlib_date)
        axs[0].set_ylim(dmin, dmax)

        # Data
        tr = trace.slice(tmin, tmax)
        axs[1].plot(tr.times(type='matplotlib'), tr.data, 'k')
        axs[1].set_xticklabels([])

        # Kurtosis
        axs[2].plot(all_mean_M.times(type='matplotlib'), all_mean_M.data, 'b')
        axs[2].set_ylabel('Kurtosis', color='b')
        axs[2].tick_params(axis='y', labelcolor='b')
        axs[2].set_xticklabels([])
        ax2b = axs[2].twinx()  # Create a twin axis with a different y-scale
        ax2b.plot(kurto_modif[0].times(type='matplotlib'),
                  kurto_modif[0].data, 'r--')
        ax2b.plot(kurto_modif[-1].times(type='matplotlib'),
                  kurto_modif[-1].data, 'r')
        ax2b.set_ylabel('Kurtosis Extrema', color='r')
        ax2b.tick_params(axis='y', labelcolor='r')
        fig.tight_layout()   # (to avoid clipping right label)
        self.kmin, self.kmax = np.min(all_mean_M.data), np.max(all_mean_M.data)

        # Energy and SNR
        nrgplot = energy.slice(tmin, tmax)
        nrgplot.data = 20 * np.log10(nrgplot.data)
        snrplot = snr.slice(tmin, tmax)
        axs[3].plot_date(snrplot.times(type='matplotlib'), snrplot.data, 'b')
        axs[3].axhline(snr_threshold, color='b', ls='--')
        axs[3].set_ylabel('SNR (dB)', color='b')
        axs[3].tick_params(axis='y', labelcolor='b')
        ax3b = axs[3].twinx()  # Create a twin axis with a different y-scale
        ax3b.plot_date(nrgplot.times(type='matplotlib'), nrgplot.data, 'r')
        ax3b.set_ylabel('Energy (dB)', color='r')
        ax3b.tick_params(axis='y', labelcolor='r')
        axs[3].set_xlabel('Time')
        plt.draw()
        plt.show(block=False)

    def station_window_add_extrema(self, picker, extrema):
        """
        Add extrema to figure 3+ axis[2]
        """
        if not self.plot_stations:
            return
        # Pick_Function.m:502
        for extremum in extrema:
            tval = picker.loop.index_to_time(extremum['i']).matplotlib_date
            ax = self.axs3[0]
            #  plt.subplot(4, 1, 1)
            ax.vlines(tval, picker.loop.dmin, picker.loop.dmax,
                      color=[0.8, 0.8, 0.8])
            ax = self.axs3[1]
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
        if not self.plot_stations:
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
        if not self.plot_stations and not self.plot_global:
            return
        if self.plot_global:
            ax2 = self.ax2
        if self.plot_stations:
            axs = self.axs3
        # Pick_Function.m:576

        # Pick Window: add picks to subplots 1, 2 and 3
        if onset_P is not None:
            t = picker.loop.index_to_time(onset_P).matplotlib_date
            if self.plot_stations:
                axs[0].vlines(t, picker.loop.dmin, picker.loop.dmax, 'b')
                axs[1].vlines(t, picker.loop.dmin, picker.loop.dmax, 'b')
                axs[2].vlines(t, self.kmin, self.kmax, 'b')
            if self.plot_global:
                ax2.vlines(t, iter-0.5, iter+0.5, color='b', ls=':')
        if onset_S is not None:
            t = picker.loop.index_to_time(onset_S).matplotlib_date
            if self.plot_stations:
                axs[0].plot(t, picker.loop.dmin, picker.loop.dmax, 'r')
                axs[1].plot(t, picker.loop.dmin, picker.loop.dmax, 'r')
                axs[2].plot(t, self.kmin, self.kmax, 'r')
            if self.plot_global:
                ax2.vlines(t, iter-0.5, iter+0.5, color='r', ls=':')

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
