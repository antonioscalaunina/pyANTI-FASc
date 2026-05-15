import numpy as np


SCALING_LAWS = {
    # Strasser et al. 2010, subduction interface
    # Mw = c + d log10(A), so A = 10 ** ((Mw - c) / d)
    # or you can use the classic log10(Y) = a + b * Mw
    "Strasser2010_interface": {
        "area": lambda Mw:  10 ** (-3.476 + 0.952 * Mw),                                #10 ** ((Mw - 4.441) / 0.846),
        "length": lambda Mw: 10 ** (-2.477 + 0.585 * Mw)                                #10 ** ((Mw - 4.868) / 1.392),
    },

    "Strasser2010_intraslab": {
    "area": lambda Mw: 10 ** (-3.225 + 0.89 * Mw),
    "aspect_ratio": 2.0,
    },

    #Murotani scaling law. It defines area from seismic moment
    #To estimate length we use the aspect ratio from Strasser
    "Murotani2013": {
        "area": lambda Mw: 1.34e-10 * (10 ** (1.5 * Mw + 9.1)) ** (2.0 / 3.0),
        "aspect_ratio": lambda Mw: (
            (10 ** (-2.477 + 0.585 * Mw))**2 /
            (10 ** (-3.476 + 0.952 * Mw))
        )
    },            

    # Wells & Coppersmith 1994
    # log10(Y) = a + b Mw
    "WC1994_all": {
        "area": lambda Mw: 10 ** (-3.49 + 0.91 * Mw),
        "length": lambda Mw: 10 ** (-2.44 + 0.59 * Mw),
    },
    "WC1994_reverse": {
        "area": lambda Mw: 10 ** (-3.99 + 0.98 * Mw),
        "length": lambda Mw: 10 ** (-2.42 + 0.58 * Mw),
    },
    "WC1994_normal": {
        "area": lambda Mw: 10 ** (-2.87 + 0.82 * Mw),
        "length": lambda Mw: 10 ** (-1.88 + 0.50 * Mw),
    },
    "WC1994_strike_slip": {
        "area": lambda Mw: 10 ** (-3.42 + 0.90 * Mw),
        "length": lambda Mw: 10 ** (-2.57 + 0.62 * Mw),
    },

    # Leonard 2014: OpenQuake implements mainly magnitude-area.
    # Define area aspect ratio.
    "Leonard2014_interplate": {
        "area": lambda Mw: 10 ** (Mw - 3.995),
        "aspect_ratio": 2.0,
    },
    "Leonard2014_SCR_dip_slip": {
        "area": lambda Mw: 10 ** (Mw - 4.185),
        "aspect_ratio": 2.0,
    },

    "Leonard2014_SCR_strike_slip": {
    "area": lambda Mw: 10 ** (Mw - 4.18),
    "aspect_ratio": 2.0,
    },

    "Thingbaijam2017_interface": {
    "area": lambda Mw: 10 ** (-3.292 + 0.949 * Mw),
     "aspect_ratio": lambda Mw: (
            (10 ** (-2.477 + 0.585 * Mw))**2 /
            (10 ** (-3.476 + 0.952 * Mw))
        )
    },

    "Thingbaijam2017_strike_slip": {
    "area": lambda Mw: 10 ** (-3.486 + 0.942 * Mw),
    "aspect_ratio": 2.0,
    },

    # Allen & Hayes 2017 - subduction interface
    # Preferred bilinear area S2 and length L.
    # Units: area in km^2, length in km.
    # Validity suggested by the paper: 7.1 <= Mw <= 9.5.
    "AllenHayes2017_interface": {
        "area": lambda Mw: (
            10 ** (-5.62 + 1.22 * Mw)
            if Mw <= 8.63
            else 10 ** (2.23 + 0.31 * Mw)
        ),
        "length": lambda Mw: 10 ** (-2.90 + 0.63 * Mw),
    },

    # Allen & Hayes 2017 - intraslab
    # Units: area in km^2, length in km.
    # Validity suggested by the paper: 7.3 <= Mw <= 8.3.
    "AllenHayes2017_intraslab": {
        "area": lambda Mw: 10 ** (-3.89 + 0.96 * Mw),
        "length": lambda Mw: 10 ** (-3.03 + 0.63 * Mw),
    },
}
