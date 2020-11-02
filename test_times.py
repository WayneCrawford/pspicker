from time import perf_counter
from obspy import read

tr = read('http://examples.obspy.org/RJOB_061005_072159.ehz.new')[0]
print(f'Running Trace.times() for {tr.stats.npts} samples')
for type in ['utcdatetime', 'matplotlib', 'timestamp', 'relative']:
    tic = perf_counter()
    times = tr.times(type)
    toc = perf_counter()
    print(f"{1.e3*(toc-tic):7.2f} ms for '{type}'")
