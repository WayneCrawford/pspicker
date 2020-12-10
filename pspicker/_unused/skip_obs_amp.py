# Generated with SMOP  0.41-beta
from smop.libsmop import *
# skip_obs_amp.m

    #### program to skip OBS amplitude remove AMP line
#### isn't seen by seisan when locate
    
    
@function
def skip_obs_amp(data_path=None,*args,**kwargs):
    varargin = skip_obs_amp.varargin
    nargin = skip_obs_amp.nargin

    fid=fopen(concat([data_path,'filenr.lis']),'rt')
# skip_obs_amp.m:6
    REA_filenames=textscan(fid,'%*6c %s','eof')
# skip_obs_amp.m:7
    
    fclose(fid)
    mkdir(concat([data_path,'/No_OBS_Amp/']))
    for i in arange(1,length(REA_filenames[1])).reshape(-1):
        fit=fopen(concat([data_path,REA_filenames[1][i]]))
# skip_obs_amp.m:13
        fot=fopen(concat([data_path,'/No_OBS_Amp/',REA_filenames[1][i]]),'w')
# skip_obs_amp.m:14
        tline=0
# skip_obs_amp.m:15
        while tline != - 1:

            tline=fgetl(fit)
# skip_obs_amp.m:18
            if tline != - 1:
                if (tline(arange(2,3)) == 'OS'):
                    if (tline(arange(11,13)) == 'AMP'):
                        continue
                    else:
                        fprintf(fot,'%s \n',tline)
                else:
                    fprintf(fot,'%s \n',tline)

        fclose(fit)
        fclose(fot)
    