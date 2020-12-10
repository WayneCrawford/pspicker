# Generated with SMOP  0.41-beta
from smop.libsmop import *
# write_rea.m

    # rea_path='/Volumes/donnees/SE_DATA/REA/CHRIS/2008/06/';
# rea_name='01-0003-53L.S200806';
# rea_path='/Volumes/donnees/SE_DATA/REA/CHRIS/2008/05/';
# rea_name='01-0008-11L.S200805';
    
    
@function
def write_rea(rea_name=None,rea_path=None,pick_cell=None,amp_cell=None,*args,**kwargs):
    varargin = write_rea.varargin
    nargin = write_rea.nargin

    # function to write an S-file with auto picks into directory
    
    # check line of the header to plot those that you really want
# The input is the cell (with stations, components and arrival times) and
# the original REA file name (to copy the header)
    
    Amp_P=copy(amp_cell)
# write_rea.m:15
    F=copy(pick_cell)
# write_rea.m:16
    fis=fopen(concat([rea_path,rea_name]),'rt')
# write_rea.m:17
    current_dir=copy(pwd)
# write_rea.m:18
    sfile_dir=fullfile(current_dir,'/Sfile_directory')
# write_rea.m:18
    if logical_not(exist(sfile_dir,'dir')):
        error('You must create a subdirectory "%s" to store new s-files',sfile_dir)
    
    fwr=fopen(fullfile(sfile_dir,rea_name),'w')
# write_rea.m:22
    
    if fwr == - 1:
        error('Could not create output s-file "%s"!',fullfile(sfile_dir,rea_name))
    
    ## Read all the text lines of S-file
    type_='None'
# write_rea.m:28
    j=1
# write_rea.m:29
    while type_ != '7':

        C[j]=cellarray([fgetl(fis)])
# write_rea.m:31
        type_=C[j](end())
# write_rea.m:32
        j=j + 1
# write_rea.m:33

    
    ## Select only desired lines
    for i in arange(1,length(C)).reshape(-1):
        a=C[i](end())
# write_rea.m:38
        if ' ' == a:
            type_1=C[i]
# write_rea.m:41
            type_1[80]='1'
# write_rea.m:42
        else:
            if '1' == a:
                type_1=C[i]
# write_rea.m:44
            else:
                if '6' == a:
                    type_6=C[i]
# write_rea.m:46
                else:
                    if 'I' == a:
                        type_I=C[i]
# write_rea.m:48
                    else:
                        if '7' == a:
                            type_7=C[i]
# write_rea.m:50
    
    fprintf(fwr,'%s\n',type_1)
    fprintf(fwr,'%s\n',type_6)
    fprintf(fwr,'%s\n',type_7)
    yy,mo,dd,ho,mi,se=datevec(F[5],nargout=6)
# write_rea.m:59
    jour=min(dd)
# write_rea.m:60
    #F{4}=F{4}+1;
    for i in arange(1,size(F[1],1)).reshape(-1):
        if dd(i) == jour + 1:
            ho[i]=ho(i) + 24
# write_rea.m:64
        #     if i==9
#         continue
#     end
        fprintf(fwr,' %-5s%2s %3s  %iA  %2i%2i %5.2f\n',F[1][i],F[2][i],F[3][i],F[4](i),ho(i),mi(i),se(i))
    
    clear('yy','mo','dd','ho','mi','se')
    yy,mo,dd,ho,mi,se=datevec(Amp_P[4],nargout=6)
# write_rea.m:74
    jour2=min(dd)
# write_rea.m:75
    for i in arange(1,size(Amp_P[1],1)).reshape(-1):
        if dd(i) == jour2 + 1:
            ho[i]=ho(i) + 24
# write_rea.m:78
        if Amp_P[5](i) < 1000000.0:
            format=' %-5s%2s  %4s A  %2i%2i %5.2f      %6.2f %4.2f\n'
# write_rea.m:81
        else:
            format=' %-5s%2s  %4s A  %2i%2i %5.2f      %6.1G %4.2f\n'
# write_rea.m:83
        fprintf(fwr,format,Amp_P[1][i],Amp_P[2][i],Amp_P[3][i],ho(i),mi(i),se(i),Amp_P[5](i),Amp_P[6](i))
    
    
    fclose(fis)
    fclose(fwr)
    return
    
if __name__ == '__main__':
    pass
    