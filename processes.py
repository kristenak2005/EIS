import datetime as dt
import os
import matplotlib.pyplot as plt
from sunpy.map import Map
import numpy as np
import matplotlib.colors as colors
import eispac.net
from eispac.net.attrs import FileType
import glob
from astropy import units as u
from sunpy.time import parse_time
from astropy.io import fits as fits
from sunpy.net import Fido, attrs as a
from aiapy.calibrate import register, update_pointing, correct_degradation
import matplotlib.gridspec as gridspec
import sunkit_image.enhance as enhance
from astropy.coordinates import SkyCoord
from general_routines import closest
import asdf

