def convert_um_to_num(concentration):
    """Convert uM concentration to # of molecules, assuming v=1.661e-12 L."""
    return concentration*1e6

def convert_nm_to_num(concentration):
    """Convert uM concentration to # of molecules, assuming v=1.661e-12 L."""
    return concentration*1e3

def convert_um_kf_to_stoch(rate):
    molar_rate = rate*1e6
    return molar_rate*1e-12

def convert_nm_kf_to_stoch(rate):
    molar_rate = rate*1e9
    return molar_rate*1e-12
