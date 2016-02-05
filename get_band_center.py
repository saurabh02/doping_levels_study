from __future__ import division

from pymatgen import Composition, Element
from pymatgen.matproj.rest import MPRester

__author__ = 'Anubhav Jain'
__copyright__ = 'Copyright 2014, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Anubhav Jain'
__email__ = 'ajain@lbl.gov'
__date__ = 'Oct 01, 2014'


def get_band_center(form):
    c = Composition(str(form))
    prod = 1.0
    for el, amt in c.get_el_amt_dict().iteritems():
        prod = prod * (Element(el).X ** amt)

    return -prod ** (1 / sum(c.get_el_amt_dict().values()))


if __name__ == '__main__':
    mpr = MPRester()
    for i in mpr.query({"elements": {"$in": ["O", "S", "Se", "Te"]}, "nelements": 2, "e_above_hull": {"$lte": 0.05}},
                       properties=['band_gap', 'formula', 'material_id', 'icsd_id']):
        c = Composition(i['formula'])
        center = get_band_center(c)
        print ','.join([str(x) for x in
                        c.formula, c.reduced_formula, center, center + i['band_gap'] / 2, center - i['band_gap'] / 2, i[
                            'material_id'], i['icsd_id']])
