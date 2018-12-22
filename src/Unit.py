import numpy as np

SI = {
    'mm':0.001, 'cm':0.01, 'm':1.0, 'km':1000., 
    'GHz':10**9, 'Hz':1.0
}

def convert_SI(val_unit, unit_in="cm", unit_out="mm"):
    val = float(val_unit.split(unit_in)[0])
    return val*SI[unit_in]/SI[unit_out]

def convert (unit_in="cm", unit_out="mm"):
    return SI[unit_in]/SI[unit_out]