import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from obspy.core.stream import Stream

from .logger import log


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

    def plot_trace(self, trace, i_station, candidates):
        """
        Plot one trace
        :param trace: trace to plot
        :param i_station: y-index of the trace (used to separate traces)
        :param candidates: candidate picks
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
        ax.vlines([x.time.matplotlib_date for x in candidates],
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
        # print(f'plot_timebounds: {rect_start_time}, {rect_end_time},
        #       {len(stations)}')
        width = rect_end_time - rect_start_time
        patch = Rectangle((rect_start_time.matplotlib_date, -0.5),
                          width=width/86400, height=N + 1, color='lightblue',
                          alpha=0.1, edgecolor=None)
        ax.axvline(rect_start_time.matplotlib_date)
        ax.axvline(rect_end_time.matplotlib_date)
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
    zorders: 0: windows: 1: trace, 2: candidates, 3: initial p-s candidates,
             4: cluster candidates, 5: final picks
    """
    def __init__(self, plot=True):
        self.plot = plot
        if not plot:
            self.fig = None
            self.ax = None
        else:
            self.fig, self.ax = plt.subplots(num="PSPicker: all picks")
        self.stations = []  # ordered list of stations

    def setup(self, starttime, endtime, stations):
        """
        :param starttime: plot start time
        :param endtime: plot end time
        :param stations: sorted list of stations
        """
        if not self.plot:
            return
        # Pick_Function.m:250
        self.stations = sorted(stations)
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

    def traces_candidates(self, traces, c_P, c_S, candidates, station):
        """
        Plot P trace(s) and candidates for one station

        :param traces: Stream of traces to plot
        :param c_P: P PickCandidate or None
        :param c_S: S PickCandidate  or None
        :param candidates: list of PickCandidates
        :param station: station name
        """
        if not self.plot:
            return
        ax = self.ax
        i_sta = self.stations.index(station)
        # Pick_Function.m:601
        for tr in traces:
            norm_trace = tr.copy().normalize()
            ax.plot_date(norm_trace.times(type="matplotlib"),
                         norm_trace.data + i_sta,
                         color='gray', ls='-', marker=None, zorder=1)
        for cand in candidates:
            ax.vlines(cand.time.matplotlib_date, i_sta-0.5, i_sta+0.5,
                      color='gray', lw=5, alpha=0.3, zorder=2.9)
            
        if c_P is not None:
            ax.vlines(c_P.time.matplotlib_date, i_sta-0.5, i_sta+0.5,
                      color='b', lw=5, alpha=0.3, zorder=3)
        if c_S is not None:
            ax.vlines(c_S.time.matplotlib_date, i_sta-0.5, i_sta+0.5,
                      color='r', lw=5, alpha=0.3, zorder=3)
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)

    def picks(self, picks, t_begin, p_clust, s_clust, o_clust):
        """
        Plot initial and final P and S picks for all stations

        :picks: all P and S PickCandidates
        :param t_begin: begin time for all plots
        :param p_clust: p-cluster times {station: time}
        :param s_clust: s-cluster times {station: time}
        :param o_clust: origin-cluster times {station: time}
        """
        if not self.plot:
            return
        # Pick_Function.m:811
        p_picks = [p for p in picks if p.phase_guess == 'P']
        s_picks = [p for p in picks if p.phase_guess == 'S']
        self._picks(p_picks, 'b')
        self._picks(s_picks, 'r')
        self._add_cluster_rectangle(p_clust, 'b', '#a0a0ff')
        self._add_cluster_rectangle(s_clust, 'r', '#ffa0a0')
        self._add_cluster_rectangle(o_clust, 'g', '#a0ffa0')
        self.ax.set_xlabel(t_begin.strftime('%Y-%m-%d'))
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)

    def _picks(self, picks, color):
        """
        Plot picks for one station as solid lines

        :parm picks: list of Picks
        :param color: color for picks
        """
        if len(picks) == 0:
            return
        for pick in picks:
            j = self.stations.index(pick.station)
            self.ax.vlines(pick.time.datetime, j - 0.5, j + 0.5, colors=color,
                           zorder=5)
        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)

    def _add_cluster_rectangle(self, cluster, color, face_color):
        if len(cluster) == 0:
            return
        for sta, t in cluster.items():
            self.ax.plot(t.matplotlib_date, self.stations.index(sta),
                         color=color, marker='x', ls=None, zorder=4)
        t = [x.matplotlib_date for x in cluster.values()]
        log(cluster, 'debug')
        log(min(t), 'debug')
        log(max(t), 'debug')
        p = Rectangle((min(t), -0.5), width=max(t) - min(t),
                      height=len(self.stations), fc=face_color, ec=color,
                      alpha=0.5, zorder=0)
        self.ax.add_patch(p)


class Station_Window():
    """
    Station window Plotter

    Detailed plot of picks, data, kurtosis, energy and SNR for one station
    Figures 3+: Station plots, each one divided in four:
        1) Picks
        2) Data
        3) Kurtosis
        4) Energy and signal-to-noise level
        5) Polarity
    """
    def __init__(self, plot=True):
        """
        Initialize station window object (but not the window)
        """
        self.plot = plot
        self.fig = None
        self.ax_picks = None
        self.ax_data = None
        self.ax_kurt = None
        self.ax_snr = None
        self.ax_pol = None

    def setup(self, trace, iter):
        """
        Plot pick phases on Figure 3+

        :param trace: trace to plot
        :param iter: station number (ignored)
        """
        if not self.plot:
            return
        # Pick_Function.m:334
        name = trace.stats.station
        self.fig, axs = plt.subplots(5, 1, sharex=True, num=name)
        self.fig.subplots_adjust(hspace=0)
        self.ax_picks = axs[0]
        self.ax_cand = axs[1]
        self.ax_kurt = axs[2]
        self.ax_snr = axs[3]
        self.ax_pol = axs[4]

        # fig.clf()
        self.ax_picks.set_title('Phase 2: Pick phases (station {})'
                                .format(trace.stats.station))

    def plot_data(self, starttime, endtime, snr_quality_thresholds,
                  trace, energy, kurto_mean, kurto_grad, dip_rect,
                  dip_thresh_P, dip_thresh_S):
        """
        Plot the waveforms, energy, kurtosis, etc

        :param starttime: starttime for all plots
        :param endtime: endtime for all plots
        :param snr_quality_thresholds: list of SNR quality thresholds
        :param trace: data trace to plot
        :param energy: WaveformEnergy object
        :param kurto_mean: mean Kurtosis
        :kurto_grad: list of cumulative kurtosis gradient Traces, for
                     different smoothings
        :dip_rect: dip-rectilinearity trace or None
        :dip_thresh_P: P-wave dip-rect threshold
        :dip_thresh_S: S-wave dip-rect threshold
        """
        if not self.plot:
            return
        # Pick_Function.m:443
        # imin, imax = np.nonzero(np.isfinite(all_mean_M))[0][[0, -1]]
        # tmin, tmax = all_mean_M.times(type="utcdatetime")[[imin, imax]]
        tmin, tmax = kurto_mean.stats.starttime, kurto_mean.stats.endtime
        self.kmin, self.kmax = np.min(kurto_mean.data), np.max(kurto_mean.data)
        kurto_grad_stream = Stream(kurto_grad).slice(tmin, tmax)

        # Plot the subplots
        self.ax_picks.set_xlim(starttime.datetime, endtime.datetime)
        self.fig.tight_layout()   # (to avoid clipping right labels)
        self._plot_picks(self.ax_picks, trace.slice(tmin, tmax))
        self._plot_cand(self.ax_cand, trace.slice(tmin, tmax))
        self._plot_kurt(self.ax_kurt, kurto_mean, kurto_grad_stream)
        self._plot_snr(self.ax_snr, energy.slice(tmin, tmax),
                       snr_quality_thresholds)
        self._plot_pol(self.ax_pol, dip_rect, dip_thresh_P, dip_thresh_S)

        self.fig.tight_layout()   # (to avoid clipping right labels)
        self.fig.subplots_adjust(hspace=0.1)

        plt.draw()
        plt.show(block=False)
        plt.pause(0.001)

    @staticmethod
    def _plot_picks(ax, trace):
        dmin, dmax = np.nanmin(trace.data), np.nanmax(trace.data)
        # patch = Rectangle((tmin.matplotlib_date, dmin),
        #                   width=(tmax - tmin)/86400, height=dmax - dmin,
        #                   color=[.9, .9, .9], edgecolor=None)
        ax.plot(trace.times(type="matplotlib"), trace.data, 'k')
        # ax.add_patch(patch)
        ax.set_ylim(dmin, dmax)
        ax.set_ylabel('Picks')

    @staticmethod
    def _plot_cand(ax, trace):
        pass

    @staticmethod
    def _plot_kurt(ax, kurto_mean, kurto_grad):
        ax.plot(kurto_mean.times(type='matplotlib'), kurto_mean.data, 'b')
        ax.set_ylabel('Kurtosis', color='b')
        ax.tick_params(axis='y', labelcolor='b')
        # ax.set_xticklabels([])
        axb = ax.twinx()  # Create a twin axis with a different y-scale
        axb.plot(kurto_grad[0].times(type='matplotlib'),
                 kurto_grad[0].data, 'r--')
        axb.plot(kurto_grad[-1].times(type='matplotlib'),
                 kurto_grad[-1].data, 'r')
        axb.set_ylabel('Extrema', color='r')
        axb.tick_params(axis='y', labelcolor='r')

    @staticmethod
    def _plot_snr(ax, energy, snr_quality_thresholds):
        nrgdB = energy.nrg.copy()
        nrgdB.data = 20 * np.log10(nrgdB.data)
        ax.plot_date(energy.nrg.times(type='matplotlib'), energy.nrg.data,
                     'r', zorder=0., alpha=0.5)
        ax.set_ylabel('Energy', color='r')
        ax.tick_params(axis='y', labelcolor='r')
        axb = ax.twinx()  # Create a twin axis with a different y-scale
        axb.plot_date(energy.snr.times(type='matplotlib'), energy.snr.data,
                      'b', zorder=1.)
        if energy.snr_threshold is not None:
            axb.axhline(energy.snr_threshold, color='b', ls='--', zorder=0.5)
        for thresh in snr_quality_thresholds:
            axb.axhline(thresh, color='b', ls='--', lw=0.1, zorder=0.5)
        axb.set_ylabel('SNR', color='b')
        axb.tick_params(axis='y', labelcolor='b')

    @staticmethod
    def _plot_pol(ax, trace, dip_thresh_P, dip_thresh_S):
        """
        :param dip_thresh_P: minimum P-wave dip-rect
        :param dip_thresh_S: maximum S-wave dip-rect
        """
        if trace is not None:
            tr = trace.copy()
            tr.data[tr.data == 0] = np.nan
            ax.plot_date(tr.times(type='matplotlib'), tr.data, 'k')
            ax.axhline(0, color='k', ls='-', lw=0.1, zorder=0.5)
            ax.axhline(dip_thresh_P, color='b', ls='--', zorder=0.5)
            ax.axhline(dip_thresh_S, color='r', ls='--', zorder=0.5)
            ax.set_ylim(-1, 1)
        else:
            ax.text(0.5, 0.5, 'Polarity not calculated',
                    ha='center', va='center', transform=ax.transAxes)
        ax.set_ylabel('Dip-Rect')

    def candidates(self, candidates, DR_thresh_P, DR_thresh_S):
        """
        Add pick candidates to all axes

        :param candidates: list of PickCandidate
        """
        if not self.plot:
            return
        # Extrema picks are grey
        color = [0.8, 0.8, 0.8]
        # Plot candidate times across  all axes
        for c in candidates:
            self.ax_picks.axvline(c.time.matplotlib_date, color=color)
            self.ax_kurt.axvline(c.time.matplotlib_date, color=color)
            self.ax_snr.axvline(c.time.matplotlib_date, color=color)
            self.ax_pol.axvline(c.time.matplotlib_date, color=color)

        # Detail candidates on dedicated axis
        ax = self.ax_cand
        axb = ax.twinx()  # Create a twin axis with a different y-scale
        ax.stem([x.time.matplotlib_date for x in candidates],
                [x.picker_value for x in candidates], markerfmt='ro',
                linefmt='r-', basefmt='r:', use_line_collection=True)
        ax.set_ylabel('picker_value', color='r')
        ax.tick_params(axis='y', labelcolor='r')
        print('Doh!')
        for x in candidates:
            print(x)
        axb.stem([x.time.matplotlib_date for x in candidates],
                 [x.DR for x in candidates], markerfmt='bx',
                 linefmt='b-', basefmt='b:', use_line_collection=True)
        axb.axhline(DR_thresh_P, color='b', ls='--')
        axb.axhline(DR_thresh_S, color='b', ls='--')
        axb.set_ylabel('dip-rect', color='b')
        axb.tick_params(axis='y', labelcolor='b')
        axb.set_ylim(-1, 1)

    def onsets(self, onset_P, onset_S, data_limits):
        """
        Plot onsets for one station, on pick and station windows

        :param onset_P: P PickCandidate or None
        :param onset_S: S PickCandidate or None
        :param data_limits: minimum and maximum data values
        """
        if not self.plot:
            return
        # Pick_Function.m:576
        # Station Window: add picks to subplots 1, 2 and 3
        if onset_P is not None:
            t = onset_P.time.matplotlib_date
            self.ax_picks.vlines(t, data_limits[0], data_limits[1], 'b')
            self.ax_kurt.vlines(t, self.kmin, self.kmax, 'b')
            # self.ax_extrem.vlines(t, data_limits[0], data_limits[1], 'b')
        if onset_S is not None:
            t = onset_S.time.matplotlib_date
            self.ax_picks.plot(t, data_limits[0], data_limits[1], 'r')
            # self.ax_extrem.plot(t, data_limits[0], data_limits[1], 'r')
            self.ax_kurt.plot(t, self.kmin, self.kmax, 'r')
