from isotope_constants import *

def init_230_232Th_a(measured_ratio, age):
    '''Calculate intitial 230/232 activity ratio from measured and age estimate'''

    return measured_ratio * np.e**(-lambda_230 * age)
