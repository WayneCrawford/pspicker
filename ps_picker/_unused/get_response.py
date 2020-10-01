# Generated with SMOP  0.41-beta
from smop.libsmop import *
from obspy.core.inventory.response.ResponseStage import PolesZerosStage
# get_response.m

    
@function
def get_response(filename, format=None):
    """
    Read response file and output PoleZeros object
    """

    if strcmp(format,'GSE'):
        Info = read_GSE(filename)
        rsample=Info.sampling
        norm_cons=Info.normal
        dig_gain=Info.amplif
        cali_gain=Info.sens
        f_ref=1
        poles=Info.poles
        zes=Info.zeros
    else:
        foc=fopen(filename,'r')
        i=1
        j=1
        k=1
        m=1
        while 1:
            tline=fgetl(foc)
            if tline == - 1:
                break
            if i <= 4:
                a[m]=str2double(tline)
                m=m + 1
            else:
                if i == 5:
                    p=str2double(tline)
                    poles=zeros(p,1)
                else:
                    if i == 6:
                        z=str2double(tline)
                        zes=zeros(z,1)
                    else:
                        if i >= 7 and i < 7 + p:
                            clear('A')
                            A=textscan(tline,'%s %s')
                            poles[j]=complex(str2double(A[1][1]),str2double(A[2][1]))
                            j += 1
                        else:
                            if i >= 7 + p and i < 7 + p + z:
                                clear('A')
                                A=textscan(tline,'%s %s')
                                zes[k]=complex(str2double(A[1][1]),str2double(A[2][1]))
                                k += 1
                            else:
                                dig_gain=str2double(tline)
            i += 1

        cali_gain=a(1)
        f_ref=1 / a(2)
        norm_cons=a(4)
        rsample=a(3)
        fclose(foc)
    
    return PolesZerosStage(
        pz_transfer_function_type='LAPLACE (HERTZ)',
        stage_gain=dig_gain*cali_gain,
        stage_gain_frequency=f_ref,
        normalization_frequency= f_ref,
        zeros=zes,
        poles=poles,
        normalization_factor= norm_cons,
        decimation_input_sample_rate=rsample)
    
if __name__ == '__main__':
    pass
    