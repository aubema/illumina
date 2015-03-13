from MetaArray import MetaArray
import math

fine_number_density = 0.999875
coarse_number_density = 0.000125

size_parameter_increment = 0.010

size_parameter_info = [
        {'name':'Type', 'cols':[
            {'name':'rural'}, 
            {'name':'urban'}, 
            {'name':'maritime'}
            ] },
        {'name':'RH', 'units':'%', 'values':[0,50,70,80,90,95,99,99]},
        {'name':'Stats', 'cols':[
                {'name':'Average radius fine'},
                {'name':'Sigma fine'},
                {'name':'Average radius coarse'},
                {'name':'Sigma coarse'}
            ]}
        ]

size_parameter = MetaArray((3,8,4),dtype=float,
        info=size_parameter_info)

size_parameter['Type':'rural', 'RH':0] = [0.02700, 0.35, 0.4300, 0.4]
size_parameter['Type':'rural', 'RH':1] = [0.02748, 0.35, 0.4377, 0.4]
size_parameter['Type':'rural', 'RH':2] = [0.02846, 0.35, 0.4571, 0.4]
size_parameter['Type':'rural', 'RH':3] = [0.03274, 0.35, 0.5477, 0.4]
size_parameter['Type':'rural', 'RH':4] = [0.03884, 0.35, 0.6462, 0.4]
size_parameter['Type':'rural', 'RH':5] = [0.04238, 0.35, 0.7078, 0.4]
size_parameter['Type':'rural', 'RH':6] = [0.04751, 0.35, 0.9728, 0.4]
size_parameter['Type':'rural', 'RH':7] = [0.05215, 0.35, 1.1755, 0.4]

size_parameter['Type':'urban', 'RH':0] = [0.02500, 0.35, 0.4000, 0.4]
size_parameter['Type':'urban', 'RH':1] = [0.02563, 0.35, 0.4113, 0.4]
size_parameter['Type':'urban', 'RH':2] = [0.02911, 0.35, 0.4777, 0.4]
size_parameter['Type':'urban', 'RH':3] = [0.03514, 0.35, 0.5805, 0.4]
size_parameter['Type':'urban', 'RH':4] = [0.04187, 0.35, 0.7061, 0.4]
size_parameter['Type':'urban', 'RH':5] = [0.04904, 0.35, 0.8634, 0.4]
size_parameter['Type':'urban', 'RH':6] = [0.05996, 0.35, 1.1691, 0.4]
size_parameter['Type':'urban', 'RH':7] = [0.06847, 0.35, 1.4858, 0.4]

size_parameter['Type':'maritime', 'RH':0] = [0.1600, 0.4, 0.1600, 0.4]
size_parameter['Type':'maritime', 'RH':1] = [0.1711, 0.4, 0.1711, 0.4]
size_parameter['Type':'maritime', 'RH':2] = [0.2041, 0.4, 0.2041, 0.4]
size_parameter['Type':'maritime', 'RH':3] = [0.3180, 0.4, 0.3180, 0.4]
size_parameter['Type':'maritime', 'RH':4] = [0.3803, 0.4, 0.3803, 0.4]
size_parameter['Type':'maritime', 'RH':5] = [0.4606, 0.4, 0.4606, 0.4]
size_parameter['Type':'maritime', 'RH':6] = [0.6024, 0.4, 0.6024, 0.4]
size_parameter['Type':'maritime', 'RH':7] = [0.7505, 0.4, 0.7505, 0.4]

def get_size_parameter(type, rh, wv, mod):
    """
    type: one of 'rural', 'urban' or 'maritime'
    rh: relative humidity supported by opac (integer)
    wv: wavelength in micron
    mod: on of 'fine' or 'coarse

    Return a tuple containing average size and sigma
    """
    rh_array = [0,50,70,80,90,95,98,99]

    stats = size_parameter[type, rh_array.index(rh)]
    if mod == 'fine':
        average_radius = stats[0]
        sigma = stats[1]
    elif mod == 'coarse':
        average_radius = stats[2]
        sigma = stats[3]
    else:
        raise Error('Mode must be fine or corse, received %s' % (mod,))

    final_average_radius = 2.0 * math.pi * average_radius / wv

    return (final_average_radius, sigma)
