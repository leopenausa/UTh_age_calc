from isotope_constants import *
from calc_Th230_U238_act import calc_Th230_U238_act

def age_calc_Th230_U238_CS(T_ini, m_Th230_U238_a, m_d234U, diff = 0.00001):
    ''' Function that iterates times until computed U/Th activity is close
    enough to measured U/Th activity based on preset threshold'''

    # Compute initial 230Th_238U activity estimate
    Th230_238U_a_estim = calc_Th230_U238_act(m_d234U, T_ini)

    # Loop until condition diference in activity ratios measured and estimated
    # is below threshold

    T_cal = T_ini

    while abs(Th230_238U_a_estim - m_Th230_U238_a) > diff:
        if (Th230_238U_a_estim - m_Th230_U238_a) <= 0:
            T_cal = T_cal + 1
        else:
            T_cal = T_cal - 1

        # Compute ratio estimate with new age to change condition
        if T_cal >= 700000:
            return 0
        Th230_238U_a_estim = calc_Th230_U238_act(m_d234U, T_cal)
        #print(T_cal, abs(Th230_238U_a_estim - m_Th230_U238_a))

    # Return age in years
    return T_cal
