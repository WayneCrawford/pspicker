# Generated with SMOP  0.41-beta
from smop.libsmop import *
# seisan2degree.m

    
@function
def seisan2degree(A_lat=None,A_lon=None,reflat=None,reflon=None,*args,**kwargs):
    varargin = seisan2degree.varargin
    nargin = seisan2degree.nargin

    # contstants
    
    coefx=88.8
# seisan2degree.m:5
    coefy=111.195
# seisan2degree.m:6
    A_lat_new=A_lat(arange(),1) + multiply(sign(A_lat(arange(),1)),A_lat(arange(),2)) / 60
# seisan2degree.m:8
    A_lon_new=A_lon(arange(),1) + multiply(sign(A_lon(arange(),1)),A_lon(arange(),2)) / 60
# seisan2degree.m:9
    y=A_lat_new - reflat
# seisan2degree.m:12
    x=A_lon_new - reflon
# seisan2degree.m:13
    y=multiply(y,coefy)
# seisan2degree.m:15
    x=multiply(x,coefx)
# seisan2degree.m:16
    return x,y
    
if __name__ == '__main__':
    pass
    