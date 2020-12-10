# Generated with SMOP  0.41-beta
from smop.libsmop import *
# write_rea_amplitude.m

    # function made to write the S-file into directory, check line of the
# header to plot those that you really want
# The input is the cell (with stations, components and arrival times) and
# the original REA file name (to copy the header)
    
    # Contrary to the write_rea, it doesn't take phase arrival cell as input
    
    # rea_path='/Volumes/donnees/SE_DATA/REA/CHRIS/2008/06/';
# rea_name='01-0003-53L.S200806';
# rea_path='/Volumes/donnees/SE_DATA/REA/CHRIS/2008/05/';
# rea_name='01-0008-11L.S200805';
    
    
@function
def write_rea_amplitude(rea_name=None,rea_path=None,amp_cell=None,remove_flag=None,*args,**kwargs):
    varargin = write_rea_amplitude.varargin
    nargin = write_rea_amplitude.nargin

    Amp_P=copy(amp_cell)
# write_rea_amplitude.m:15
    fis=fopen(concat([rea_path,rea_name]),'rt')
# write_rea_amplitude.m:16
    current_dir=copy(pwd)
# write_rea_amplitude.m:17
    fwr=fopen(concat([current_dir,'/Sfile_directory/',rea_name]),'w')
# write_rea_amplitude.m:18
    if fwr == - 1:
        disp('wrong path')
    
    ## Read all the text lines of S-file
    type_='None'
# write_rea_amplitude.m:24
    j=1
# write_rea_amplitude.m:25
    while type_ != '7':

        C[j]=cellarray([fgetl(fis)])
# write_rea_amplitude.m:27
        type_=C[j](end())
# write_rea_amplitude.m:28
        j=j + 1
# write_rea_amplitude.m:29

    
    ## Select only desired lines
    for i in arange(1,length(C)).reshape(-1):
        a=C[i](end())
# write_rea_amplitude.m:34
        if ' ' == a:
            type_1=C[i]
# write_rea_amplitude.m:37
            type_1[80]='1'
# write_rea_amplitude.m:38
        else:
            if '1' == a:
                type_1=C[i]
# write_rea_amplitude.m:40
            else:
                if '6' == a:
                    type_6=C[i]
# write_rea_amplitude.m:42
                else:
                    if 'I' == a:
                        type_I=C[i]
# write_rea_amplitude.m:44
                    else:
                        if '7' == a:
                            type_7=C[i]
# write_rea_amplitude.m:46
                            k=copy(i)
# write_rea_amplitude.m:47
    
    fprintf(fwr,'%s\n',type_1)
    fprintf(fwr,'%s\n',type_6)
    fprintf(fwr,'%s\n',type_7)
    k=0
# write_rea_amplitude.m:55
    while logical_not(feof(fis)):

        line=fgetl(fis)
# write_rea_amplitude.m:57
        if remove_flag == 1:
            if strcmp(line(arange(11,12)),'IA') or strcmp(line(11),'A'):
                continue
            else:
                fprintf(fwr,'%s\n',line)
                k=k + 1
# write_rea_amplitude.m:63
                aaa[k]=ftell(fwr)
# write_rea_amplitude.m:64
        else:
            fprintf(fwr,'%s\n',line)
            k=k + 1
# write_rea_amplitude.m:68
            aaa[k]=ftell(fwr)
# write_rea_amplitude.m:69

    
    if logical_not(isempty(Amp_P)):
        fseek(fwr,aaa(end() - 1),'bof')
        clear('yy','mo','dd','ho','mi','se')
        yy,mo,dd,ho,mi,se=datevec(Amp_P[4],nargout=6)
# write_rea_amplitude.m:78
        jour2=min(dd)
# write_rea_amplitude.m:79
        for i in arange(1,size(Amp_P[1],1)).reshape(-1):
            if dd == jour2 + 1:
                ho=ho + 24
# write_rea_amplitude.m:82
            if Amp_P[5](i) < 100000.0:
                format=' %s%s  %4s A  %2i%2i %5.2f      %6.2f %4.2f\n'
# write_rea_amplitude.m:85
            else:
                format=' %s%s  %4s A  %2i%2i %5.2f      %6.1G %4.2f\n'
# write_rea_amplitude.m:87
            fprintf(fwr,format,Amp_P[1][i],Amp_P[2][i],Amp_P[3][i],ho(i),mi(i),se(i),Amp_P[5](i),Amp_P[6](i))
    
    
    fclose(fis)
    fclose(fwr)
    return
    
if __name__ == '__main__':
    pass
    