from obspy import read
from obspy.core.stream import Stream

tr = read('http://examples.obspy.org/RJOB_061005_072159.ehz.new')[0]
tr2 = tr.copy()
tr3 = tr.copy()
# tr3.data = tr3.data[:-1]
tr3.stats.starttime += 0.001
st1 = Stream(traces=[tr, tr2])
st2 = Stream(traces=[tr, tr3])
print(st1)
print(st2)
st1.stack()
st2.stack()
print(st1)
print(st2)
