import matplotlib
# matplotlib.use('TkAgg')
import numpy as np
import csv
import scipy.integrate as integrate
import scipy.constants as const
import scipy.stats as stats
from PyAstronomy import pyasl
from aflare import aflare1
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
from astropy.io import fits,ascii
import florian_flares
from astropy.stats import sigma_clip
from scipy import interpolate
import os
import math
import warnings
from copy import deepcopy
from functools import partial
from scipy.optimize import least_squares