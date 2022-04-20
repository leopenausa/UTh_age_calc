####### IMPORTS ##############
import json
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from calc_Th230_U238_act import calc_Th230_U238_act
from age_calc_Th230_U238_CS import age_calc_Th230_U238_CS
from init_230_232Th_a import init_230_232Th_a
from m_Th230_U238_a_from_concs import m_Th230_U238_a_from_concs
from isotope_constants import *
###############################

########################################################################3####
# Use the full page instead of a narrow central column
st.set_page_config(layout="wide")

######## MAIN PAGE LAYOUT #############
st.title('U/Th age calculation')

st.markdown("""
This app is for age calculation from **U238**, **Th230** concentrations and age corrections
based on Th230""")

######## SIDEBAR LAYOUT ##############
st.sidebar.header('Enter isotope concentrations')
st.sidebar.markdown("""
***
Enter the isotope concentrations and their respective
units. Also enter uncertainties
***
""")

smp_name = st.sidebar.text_input('Enter sample name')

col1, col2 = st.sidebar.columns([2,1])

U238_tmp = col1.number_input('238U')
U238_units = col2.selectbox('units', ['ppb', 'pmol'], key = 1)
# Unit conversion
if U238_units == 'ppb':
    U238_C = U238_tmp
elif U238_units == 'pmol':
    U238_C = U238_tmp * 1000 / U238_atom

U238_err = col1.number_input('238U uncertainty')
U238_err_units = col2.selectbox('units', ['rsd', 'abs'], key = 2)
#st.write(U238_C)
st.sidebar.markdown("""
***
""")

col3, col4 = st.sidebar.columns([2,1])

Th230_tmp = col3.number_input('230Th')
Th230_units = col4.selectbox('units', ['pmol', 'ppt'], key = 3)
# Unit conversion
if Th230_units == 'pmol':
    Th230_C = Th230_tmp
elif Th230_units == 'ppt':
    Th230_C = Th230_tmp / Th230_atom

Th230_err = col3.number_input('230Th uncertainty')
Th230_err_units = col4.selectbox('units', ['rsd', 'abs'], key = 4)
#st.write(Th230_C)
st.sidebar.markdown("""
***
""")

col5, col6 = st.sidebar.columns([2,1])

Th232_tmp = col5.number_input('232Th')
Th232_units = col6.selectbox('units', ['pmol', 'ppt'], key = 5)
# Unit conversion
if Th232_units == 'pmol':
    Th232_C = Th232_tmp
elif Th232_units == 'ppt':
    Th232_C = Th232_tmp / Th232_atom

Th232_err = col5.number_input('232Th uncertainty')
Th232_err_units = col6.selectbox('units', ['rsd', 'abs'], key = 6)
#st.write(Th232_C)
st.sidebar.markdown("""
***
""")

col7, col8 = st.sidebar.columns([2,1])

d234U_C = col7.number_input('d234U')
d234U_units = col8.selectbox('units', ['permil', 'pmol'], key = 7)

d234U_err = col7.number_input('d234U uncertainty')
d234U_err_units = col8.selectbox('units', ['rsd', 'abs'], key = 8)

st.sidebar.markdown("""
***
""")


####################################################################3####
# Execute on button press

############################  Load data files for isolines ############################
#@st.cache
def import_isolines(act_iso, d234U_iso, age_iso):

    with open(act_iso, "r") as f:
        res_final = json.load(f)

    with open(d234U_iso, "r") as f:
        delt_final = json.load(f)

    with open(age_iso, "r") as f:
        age_iso = json.load(f)

    return res_final, delt_final, age_iso

res_final, delt_final, age_iso = import_isolines("Th230_U238_act_isolines.txt",
                                                "d234U_measured_isolines.txt",
                                                "age_isolines.txt")

if st.sidebar.button('Calculate corrected age'):

    # Read conc values from number input widgets
    U238_ppb = U238_C
    Th230_pmol = Th230_C
    Th232_pmol = Th232_C
    m_d234 = d234U_C
    ########################################################################

    # Calculate initial 230/238 activity ratio from concs
    ini_230_238_a = m_Th230_U238_a_from_concs(Th230_pmol, U238_ppb)

    # Determine intial age estimate
    ini_Yr = age_calc_Th230_U238_CS(0, ini_230_238_a, m_d234,diff=0.000005)
    ages = [ini_Yr]
    #initial_age = [ini_Yr]

    for i in range(5):
        #st.progress(i)
        # initial 230
        old_Yr = ini_Yr
        Th230_232_now = Th230_Th232_ini * np.e**(-lambda_230 * old_Yr)
        conc_230 = Th230_pmol - (Th230_232_now * Th232_pmol)
        ini_230_238_a = m_Th230_U238_a_from_concs(conc_230, U238_ppb)
        ini_Yr = age_calc_Th230_U238_CS(0, ini_230_238_a, m_d234, diff=0.00001)

        #print(Th230_232_now)
        ages.append(ini_Yr)

    st.write(pd.DataFrame({
     'Sample name': smp_name,
     '[238U] (ppb)': [U238_ppb],
     '[230Th] (pmol)': [Th230_pmol],
     '[232Th] (pmol)': [Th232_pmol],
     'initial 230/238': [ini_230_238_a],
     'corrected age': [ages[-1]]
     }))


    ########################################################################################
    fig = plt.figure(figsize=(20,20))
    ax_dict = fig.subplot_mosaic(
    [['R1', 'R3', 'R3', 'R3'],
     ['R2', 'R3', 'R3', 'R3'],
     ['', 'R3', 'R3', 'R3']
     ],
    )

    ##################### SOME CALCULATIONS###############################################3
    time = np.arange(0, 600000, 1)

    res = []
    res_2 = []
    for i in time:
        res.append(calc_Th230_U238_act(0, i))
        res_2.append(calc_Th230_U238_act(m_d234, i))



    ############################ Gridding for isolines  #############################################
    Th230U238_A = np.arange(0, 3.5, 0.05)
    delta234U_m = np.arange(0, 4000, 10)

    X, Y = np.meshgrid(Th230U238_A, delta234U_m)
    ###############################################################################################

    # Ratio Box 1
    R1_plt = ax_dict['R1'].plot(ages)
    ax_dict['R1'].set_xlabel('Nº iterations', size = 20)
    ax_dict['R1'].set_ylabel('Age correction yrs', size = 20)
    ax_dict['R1'].set_title('232Th correction', size = 20)

    # Ratio Box 2
    R2_plt = ax_dict['R2'].plot(res, 'k')
    R2_plt = ax_dict['R2'].plot(res_2)
    R2_plt = ax_dict['R2'].plot([0, ages[-1]], [ini_230_238_a, ini_230_238_a], ':r')
    R2_plt = ax_dict['R2'].plot([ages[-1], ages[-1]], [0, ini_230_238_a], ':r')
    ax_dict['R2'].set_xlabel('230/238 activity', size = 20)
    ax_dict['R2'].set_ylabel('measured d234U', size = 20)
    ax_dict['R2'].set_title('Age solution', size = 20)

    # Ratio Box 3
    for i in range(8):
        R3_plt = ax_dict['R3'].plot(res_final[i], delt_final[i], 'grey')

    R3_plt = ax_dict['R3'].contour(X,Y,age_iso,[25, 50, 75, 100, 150, 200, 300, 400, 500], colors='k', linestyles='dotted')
    R3_plt = ax_dict['R3'].clabel(R3_plt, inline=True, fontsize=20)
    ax_dict['R3'].set_xlim(0, 3.5)
    ax_dict['R3'].set_ylim(0, 3500)

    ##################### PLOT SAMPLE RESULTS  ##############################

    R3_plt = ax_dict['R3'].scatter(ini_230_238_a, m_d234, s=40, linewidths=10, edgecolors='red')
    R3_plt = ax_dict['R3'].contour(X,Y,age_iso, [ages[-1]/1000], colors='red', linestyles='solid')
    R3_plt = ax_dict['R3'].clabel(R3_plt, inline=True, fontsize=20)

    tmp=[]
    acts=[]
    initial_234 = m_d234 * np.e**(lambda_234 * ages[-1])

    for i in range(ages[-1]):
        intermediate =initial_234 * np.e**(-lambda_234 *i)
        tmp.append(intermediate)
        acts.append(calc_Th230_U238_act(intermediate, i))

    R3_plt = ax_dict['R3'].plot(acts, tmp, 'red')

    st.pyplot(fig)

    # st.write(pd.DataFrame({
    #  'Sample name': ['No name'],
    #  'initial 230/238': [ini_230_238_a],
    #  'initial age': [initial_age],
    #  'corrected age': [ages[-1]]
    #  }))
