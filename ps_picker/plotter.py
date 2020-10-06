import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


class Plotter():
    """
    Class to plot figures representing the work in PSPicker
    """
    def __init__(self, plot_global=True, plot_stations=True):
        """
        :param plot_global: show global selection and picks plots
        :kind plot_global: bool
        :param plot_stations: show individual station plots
        :kind plot_stations: bool
        """
        self.gw = Global_Window(plot_global)
        self.pw = Picks_Window(plot_global)
        self.sw = Station_Window(plot_stations)


class Global_Window():
    """
    Global window Plotter

    Shows the entire time span for all stations, the initial picks and the
    selected pick window
    """
    def __init__(self, plot=True):
        self.plot = plot
        if not plot:
            self.fig = None
            self.ax = None
        else:
            self.fig, self.ax = plt.subplots(num="PSPicker: Overview")

    def plot_trace(self, trace, i_station, extrema_t):
        """
        Plot one trace
        :param trace: trace to plot
        :param i_station: y-index of the trace (used to separate traces)
        :param extrema_t: offset time of all extrema for this trace
        """
        # Pick_Function.m: ~160?
        if not self.plot:
            return
        ax = self.ax
        norm_trace = trace.copy().detrend().normalize()
        ax.plot_date(norm_trace.times(type="matplotlib"),
                     (0.5 * norm_trace.data) + i_station,
                     fmt='-',
                     label=trace.stats.station)
        ax.vlines([t.matplotlib_date for t in extrema_t],
                  i_station + 0.45, i_station - 0.45, color='k')
        plt.draw()
        plt.show(block=False)

    def plot_timebounds(self, rect_start_time, rect_end_time, stations):
        """
        Plot the global time bounds and set axis limits

        :param rect_start_time: global window start time
        :param rect_end_time: global window end time
        :param stations: list of stations in same order as iTrace above
        """
        if not self.plot:
            return
        # Pick_Function.m:231
        ax = self.ax
        N = len(stations)
        width = rect_end_time - rect_start_time
        patch = Rectangle((rect_start_time.matplotlib_date, -0.5),
                          width=width/86400, height=N + 1, color='lightblue',
                          alpha=0.1, edgecolor=None)
        ax.add_patch(patch)
        ax.set_title('Phase 1: Calculate global window')
        ax.set_xlabel(rect_start_time.strftime('%Y-%m-%d'))
        ax.set_ylim(-0.5, N + 0.5)
        ax.set_yticks(np.arange(0, N))
        ax.set_yticklabels(stations)
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)


class Picks_Window():
    """
    Picks window Plotter

    Shows picks in the selected pick window, for all stations
    """
    def __init__(self, plot=True):
        self.plot = plot
        if not plot:
            self.fig = None
            self.ax = None
        else:
            self.fig, self.ax = plt.subplots(num="PSPicker: all picks")

    def setup(self, starttime, endtime, stations):
        """
        :param starttime: plot start time
        :param endtime: plot end time
        :param stations: sorted list of stations
        """
        if not self.plot:
            return
        # Pick_Function.m:250
        N = len(stations)
        # fig.clf()
        self.ax.set_title('Phase 3: Select final picks based on clustering')
        # ax.set_xlim(starttime.matplotlib_date, endtime.matplotlib_date)
        self.ax.set_xlim(starttime.datetime, endtime.datetime)
        self.ax.set_ylim(-0.5, N-0.5)
        self.ax.set_yticks(np.arange(0, N))
        self.ax.set_yticklabels(stations)
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)

    def Ptraces_onsets(self, traces, onset_P, onset_S, i_station):
        """
        Plot P trace(s) for one station

        :param traces: Stream of traces to plot
        :param onset_P: P onset time or None
        :param onset_S: S onset time or None
        :param i_station: sequence number of the station
        """
        if not self.plot:
            return
        ax = self.ax
        # Pick_Function.m:601
        for tr in traces:
            norm_trace = tr.copy().normalize()
            ax.plot_date(norm_trace.times(type="matplotlib"),
                          norm_trace.data + i_station, 'k', zorder=5.)
        if onset_P is not None:
            ax.vlines(onset_P.matplotlib_date, i_station-0.5, i_station+0.5,
                      color='b', lw=5, alpha=0.3, zorder=1.)
        if onset_S is not None:
            ax.vlines(onset_S.matplotlib_date, i_station-0.5, i_station+0.5,
                      color='r', lw=5, alpha=0.3, zorder=1.)
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)

    def picks(self, picks, stations, t_begin, P_clust_s, S_clust_s):
        """
        Plot initial and final P and S picks for all stations

        :picks: all picks
        :stations: ordered list of stations.  The y position depends on the
                   position of the pick's station in this list
        :param t_begin: begin time for all plots
        :param P_clust_s: width of P cluster window in seconds
        :param S_clust_s: width of S cluster window in seconds
        """
        if not self.plot:
            return
        # Pick_Function.m:811
        p_picks = [p for p in picks if p.phase_hint[0] == 'P']
        s_picks = [p for p in picks if p.phase_hint[0] == 'S']
        self._picks(p_picks, 'b', '#a0a0ff', P_clust_s, stations)
        self._picks(s_picks, 'r', '#ffa0a0', S_clust_s, stations)
        self.ax.set_xlabel(t_begin.strftime('%Y-%m-%d'))
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)

    def _picks(self, picks, color, face_color, width, stations):
        """
        Plot picks for one station as solid lines
        
        :parm picks: list of Picks
        :param color: color for picks
        :param face_color: face color for pick window
        :parm width: width of pick window in seconds
        :param stations: ordered list of stations
        """
        ax = self.ax
        if len(picks) > 0:
            window_start = picks[0].time
            for pick in picks:
                if pick.time < window_start:
                    window_start = pick.time
                j = stations.index(pick.waveform_id.station_code)
                # ax.vlines(pick.time.matplotlib_date,
                ax.vlines(pick.time.datetime, j - 0.5, j + 0.5, colors=color,
                          zorder=10.)
            p = Rectangle((window_start.matplotlib_date, -0.5),
                           width=width / 86400, height=len(stations),
                           fc=face_color, ec=color, alpha=0.5, zorder=0.)
            ax.add_patch(p)
            # ax.add_collection(PatchCollection([p]))
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)


class Station_Window():
    """
    Station window Plotter

    Detailed plot of picks, data, kurtosis, energy and SNR for one station
    Figures 3+: Station plots, each one divided in four:
        1) Picks
        2) Data
        3) Kurtosis
        4) Energy and signal-to-noise level
    """
    def __init__(self, plot=True):
        """
        Initialize station window object (but not the window)
        """
        self.plot = plot
        self.fig = None
        self.axs = None

    def setup(self, trace, iter):
        """
        Plot pick phases on Figure 3+

        :param trace: trace to plot
        :param iter: station iteration number
        """
        if not self.plot:
            return
        # Pick_Function.m:334
        # fig, axs = plt.subplots(4, 1, num=2+iter)
        self.fig, self.axs = plt.subplots(4, 1, sharex=True,
                                          num=trace.stats.station)
        # fig.clf()
        self.axs[0].plot(trace.times(type="matplotlib"), trace.data, 'k')
        self.axs[0].set_title('Phase 2: Pick phases (station {})'
                              .format(trace.stats.station))
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)

    def data(self, starttime, endtime, snr_threshold, trace, snr, energy,
                all_mean_M, kurto_grad):
        """
        Plot the data, snr, energy, kurtosis, etc
        
        :param starttime: starttime for all plots
        :param endtime: endtime for all plots
        :param snr_threshold: minimum acceptable SNR
        :param trace: data trace to plot
        :param snr: signal-to-noise level trace
        :param energy: energy trace
        :param all_mean_M: mean cumulative Kurtosis
        :kurto_grad: cumulative kurtosis gradients, for different smoothings
        """
        if not self.plot:
            return
        # Pick_Function.m:443
        axs = self.axs
        fig = self.fig
        imin, imax = np.nonzero(np.isfinite(all_mean_M))[0][[0, -1]]
        tmin, tmax = all_mean_M.times(type="utcdatetime")[[imin, imax]]
        dmin, dmax = np.nanmin(trace.data), np.nanmax(trace.data)
        # Pick windows
        patch = Rectangle((tmin.matplotlib_date, dmin),
                          width=(tmax - tmin)/86400, height=dmax - dmin,
                          color=[.9, .9, .9], edgecolor=None)
        axs[0].add_patch(patch)
        axs[0].set_xlim(starttime.datetime, endtime.datetime)
        axs[0].set_ylim(dmin, dmax)
        axs[0].set_ylabel('Picks')

        # Data subplot
        tr = trace.slice(tmin, tmax)
        axs[1].plot_date(tr.times(type='matplotlib'), tr.data, 'k')
        axs[1].set_xticklabels([])
        axs[1].set_ylabel('Data')

        # Kurtosis
        axs[2].plot(all_mean_M.times(type='matplotlib'), all_mean_M.data, 'b')
        axs[2].set_ylabel('Kurtosis', color='b')
        axs[2].tick_params(axis='y', labelcolor='b')
        axs[2].set_xticklabels([])
        ax2b = axs[2].twinx()  # Create a twin axis with a different y-scale
        ax2b.plot(kurto_grad[0].times(type='matplotlib'),
                  kurto_grad[0].data, 'r--')
        ax2b.plot(kurto_grad[-1].times(type='matplotlib'),
                  kurto_grad[-1].data, 'r')
        ax2b.set_ylabel('Extrema', color='r')
        ax2b.tick_params(axis='y', labelcolor='r')
        fig.tight_layout()   # (to avoid clipping right label)
        self.kmin, self.kmax = np.min(all_mean_M.data), np.max(all_mean_M.data)

        # Energy and SNR
        nrgplot = energy.slice(tmin, tmax)
        nrgplot.data = 20 * np.log10(nrgplot.data)
        snrplot = snr.slice(tmin, tmax)
        ax3b = axs[3].twinx()  # Create a twin axis with a different y-scale
        axs[3].plot_date(nrgplot.times(type='matplotlib'), nrgplot.data, 'r',
                       zorder=0., alpha=0.5)
        axs[3].set_ylabel('Energy (dB)', color='r')
        axs[3].tick_params(axis='y', labelcolor='r')
        ax3b.plot_date(snrplot.times(type='matplotlib'), snrplot.data, 'b',
                         zorder=1.)
        ax3b.axhline(snr_threshold, color='b', ls='--', zorder=0.5)
        ax3b.set_ylabel('SNR', color='b')
        ax3b.tick_params(axis='y', labelcolor='b')
        ax3b.set_xlabel('Time')
        plt.draw()
        plt.show(block=False)

    def extrema(self, extrema_t, data_limits):
        """
        Add extrema to figure 3+ axis[2]

        :param extrema: extrema times
        :param data_limits: minimum and maximum data values
        """
        if not self.plot:
            return
        # Pick_Function.m:502
        dmin, dmax = data_limits[0], data_limits[1]
        for t in extrema_t:
            ax = self.axs[0]    # plt.subplot(4, 1, 1)
            ax.vlines(t.matplotlib_date, dmin, dmax, color=[0.8, 0.8, 0.8])
            ax = self.axs[1]   # plt.subplot(4, 1, 2)
            ax.vlines(t.matplotlib_date, dmin, dmax, color=[0.8, 0.8, 0.8])

    def onsets(self, onset_P, onset_S, data_limits):
        """
        Plot onsets for one station, on pick and station windows

        :param onset_P: P onset time or None
        :param onset_S: S onset time or None
        :param data_limits: minimum and maximum data values
        """
        if not self.plot:
            return
        axs = self.axs
        # Pick_Function.m:576
        # Station Window: add picks to subplots 1, 2 and 3
        if onset_P is not None:
            t = onset_P.matplotlib_date
            axs[0].vlines(t, data_limits[0], data_limits[1], 'b')
            axs[1].vlines(t, data_limits[0], data_limits[1], 'b')
            axs[2].vlines(t, self.kmin, self.kmax, 'b')
        if onset_S is not None:
            t = onset_S.matplotlib_date
            axs[0].plot(t, data_limits[0], data_limits[1], 'r')
            axs[1].plot(t, data_limits[0], data_limits[1], 'r')
            axs[2].plot(t, self.kmin, self.kmax, 'r')
