import numpy as np

## Calculation of Fine Fuel Moisture Code (FFMC)
def calculate_ffmc(temp, rh, precip, wind, ffmc_prev):
    """
    Parameters:
    temp: Temperature at noon [°C]
    rh: relative humidity at noon [%]
    precip: precipitation in mm
    wind: wind speed at noon [m/s]
    ffmc_prev: previous days FFMC value

    Return:
    ffmc: Fine Moisture Fuel Code
    """
    
    # Fine fuel moisture content (mt−1) from the previous day
    mo = 147.2 * ((101 - ffmc_prev) / (59.5 + ffmc_prev))

    if precip > 0.5:
        rf = precip - 0.5
        if mo <= 150:
            mr = mo + 42.5 * rf * np.exp(-100/ (251-mo)) * (1 - np.exp(-6.93/rf))
        else: #(mo > 150)
            mr = mo + 42.5 * rf * (np.exp(-100/ (251-mo))) * (1 - np.exp(-6.93/rf)) + 0.0015 * (mo - 150)**2 * rf**0.5
        if mr > 250:
            mr = 250
        mo = mr

    Ed = 0.942 * rh**0.679 + 11 * np.exp((rh - 100) /10) + 0.18 * (21.1 - temp) * (1 - np.exp(-0-115*rh))

    if mo > Ed:
        ko = 0.424 * (1 - (rh/100)**1.7) + 0.0694 * wind**0.5 * (1 - (rh/100)**8)
        kd = ko * 0.581 * np.exp(0.0365 * temp)
        m = Ed + (mo - Ed) * 10**(-kd)
    else: #(mo < Ed)
        Ew = 0.618 * rh**0.753 + 10 * np.exp((rh - 100)/10) + 0.18 * (21.1 - temp) * (1 - np.exp(-0-115*rh))
        if mo < Ew:
            k1 = 0.424 * (1 - ((100 - rh)/100)**1.7) + 0.0694 * wind * (1 - ((100 - rh)/100)**8)
            kw = k1 * 0.581 * np.exp(0.0365 * temp)
            m = Ew - (Ew - mo) * 10**(-kw)
        else:
            m = mo

    ffmc = 59.5 * (250 - m) / (147.2 + m)
    if (ffmc > 101.0):
        ffmc = 101.0
    if (ffmc <= 0.0): 
        ffmc = 0.0
    return ffmc


## Calculate Duff Moisture Code
def calculate_dmc(temp, rh, precip, dmc_prev, month):

    if precip <= 1.5:
        Pr = dmc_prev
    else:
        re = 0.92 * precip - 1.27
        # Mo = 20 + np.exp(5.6348 - (dmc_prev/43.43))
        Mo = 20.0 + 280.0/np.exp(0.023*dmc_prev)
        if dmc_prev <= 33:
            b = 100 / (0.5 + 0.3* dmc_prev)
        elif 33 < dmc_prev <= 65:
            b = 14 - 1.3* np.log(dmc_prev)
        else:  # This covers the case where p > 65
            b = 6.2 * np.log(dmc_prev) - 17.2

        Mr = Mo + 1000 * re / (48.77 + b * re)
        Pr = 244.72 - 43.43 * np.log(Mr - 20)
        # Pr = 43.43 * (5.6348 - np.log(Mr-20.0))

    day_length = [6.5, 7.5, 9.0, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 8.0, 7.0, 6.0]
    Le = day_length[month - 1]

    if temp >= -1.1:
        K = 1.894 * (temp + 1.1) * (100 - rh) * Le * 10**(-6)
    else:
        K = 0.0

    if Pr < 0.0:
        Pr = 0.0
        
    dmc = Pr + 100*K # in previous formula it's 100*K
    if dmc <= 0.0:
        dmc = 0.0
    return dmc

##### Drought Code #####
def calculate_dc(temp, precip, dc_prev, month):
    Do = dc_prev
    if precip > 2.8:
        rd = 0.83 * precip - 1.27
        Qo = 800 * np.exp(-Do/400)
        Qr = Qo + 3.937 * rd
        Dr = 400 * np.log(800/Qr)
        if Dr > 0.0:
            Do = Dr
        else:
            Do = 0.0

    # Length Factor array for months 1 to 12
    length_factors = [-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6]
    Lf = length_factors[month - 1]

    if temp < -2.8:
        V = Lf
    else:
        V = (0.36 * (temp + 2.8)) + Lf
    if V <= 0.0:
        V = 0.0
    dc = Do + 0.5*V
    return dc


## Initial speed index
def calculate_isi(wind, ffmc):
    W = np.exp(0.05039*wind)
    m = 147.2 * (101 - ffmc) / (59.5 + ffmc)
    F = (91.9*np.exp(-0.1386*m)) * (1 + (m**5.31)/(4.93*10**7))
    isi = 0.208 * W * F
    return isi


## Buildup Index
def calculate_bui(dmc, dc):
    denominator = dmc + 0.4 * dc
    if denominator == 0:  # Handle division by zero
        return 0.0

    if dmc <= 0.4 * dc:
        bui = 0.8 * (dmc * dc) / denominator
    else:
        bui = dmc - (1 - 0.8 * dc / denominator) * (0.92 + (0.0114 * dmc) ** 1.7)

    return max(bui, 0.0)  # Ensure BUI is non-negative




## Fire Weather Index
## calculate FWI
def calculate_fwi(bui, isi):
    if bui <= 80:
        D = 0.626 * bui**0.809 + 2
    else:
        D = 1000 / (25 + 108.64 * np.exp(-0.023 * bui))

    B = 0.1 * isi * D

    if B > 1:
        fwi = np.exp(2.72 * (0.434 * np.log(B))**0.647)
    else:
        fwi = B
    return fwi
