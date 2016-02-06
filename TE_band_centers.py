from pymatgen import Composition, Element
from pymatgen.matproj.rest import MPRester

__author__ = 'Saurabh Bajaj'
# __copyright__ = 'Copyright 2014, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Saurabh Bajaj'
__email__ = 'sbajaj@lbl.gov'
__date__ = 'Feb 05, 2016'


def get_band_center(form):
    c = Composition(str(form))
    prod = 1.0
    for el, amt in c.get_el_amt_dict().iteritems():
        prod = prod * (Element(el).X ** amt)

    return -prod ** (1 / sum(c.get_el_amt_dict().values()))


if __name__ == '__main__':
    mpr = MPRester()
