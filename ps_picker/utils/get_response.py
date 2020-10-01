# Generated with SMOP  0.41-beta
# from smop.libsmop import *
import numpy as np
# from obspy.core.inventory.response import PolesZerosResponseStage


def get_response(filename, format=None):
    """
    Read response file and output PoleZeros object
    """
    if format == 'GSE':
        info = read_GSE(filename)
    else:
        info = read_Baillard_PZ(filename)
    return {'gain': info['amplif'] * info['sens'],
            'poles': info['poles'],
            'sensitivity': 1.0,
            'zeros': info['zeros']}
#    return PolesZerosResponseStage(
#        pz_transfer_function_type='LAPLACE (HERTZ)',
#        stage_gain=info['amplif'] * info['sens'],
#        stage_gain_frequency=info['f_ref'],
#        normalization_frequency=info['f_ref'],
#        zeros=info['zeros'],
#        poles=info['poles'],
#        normalization_factor=info['normal'],
#        decimation_input_sample_rate=info['sampling'],
#        stage_sequence_number=1,
#        input_units='m/s',
#        output_units='counts')


def read_Baillard_PZ(filename):
    """
    read funky Baillard PZ format
    """
    info = {}
    with open(filename, 'r') as foc:
        i = 0
        a, info['poles'], info['zeros'] = [], [], []
        for tline in foc:
            if i < 4:
                a.append(float(tline))
            elif i == 4:
                n_poles = float(tline)
            elif i == 5:
                n_zeros = float(tline)
            elif i >= 6 and i < (6 + n_poles):
                A = tline.split()
                info['poles'].append(np.complex(float(A[0]), float(A[1])))
            elif i >= 6 + n_poles and i < 6 + n_poles + n_zeros:
                A = tline.split()
                info['zeros'].append(np.complex(float(A[0]), float(A[1])))
            else:
                info['amplif'] = float(tline)
            i += 1
    info['sens'] = a[0]
    info['f_ref'] = 1 / a[1]
    info['sampling'] = a[2]
    info['normal'] = a[3]
    return info


# read_GSE.m
def read_GSE(name=None):
    """
    read a GSE response file

    :param name:    Name of GSE2 response file
    :returns: Structure with all parameters
    """
# read_GSE.m:6
    info = dict()
    with open(name, 'r') as fic:
        for line in fic:
            if 'CAL' == line[:3]:
                info['station'] = line[6:10]
                info['type'] = line[20:29].trim()
                info['sens'] = float(line[29:40])
                info['sampling'] = float(line[46:56])
            else:
                if 'PAZ' == line[:3]:
                    info['normal'] = float(line[11:25])
                    num_poles = float(line[41:43])
                    num_zeros = float(line[45])
                    info['poles'] = []
                    info['zeros'] = []
                    for p in np.arange(num_poles):
                        A = next(fic).split()
                        info['poles'].extend(np.complex(float(A[0]),
                                                        float(A[1])))
                    for z in np.arange(num_zeros):
                        A = next(fic).split()
                        info['zeros'].extend(np.complex(float(A[0]),
                                                        float(A[1])))
                else:
                    if 'DIG' == line[:3]:
                        info['amplif'] = float(line[9:23])
                        info['digitizer'] = line[36:42].trim()
                    else:
                        continue
    info['f_ref'] = 1
    return info
