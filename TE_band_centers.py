from pymatgen import Composition, Element
from pymatgen.matproj.rest import MPRester
import csv
import pandas

__author__ = 'Saurabh Bajaj'
# __copyright__ = 'Copyright 2014, The Materials Project'
__version__ = '0.1'
__maintainer__ = 'Saurabh Bajaj'
__email__ = 'sbajaj@lbl.gov'
__date__ = 'Feb 05, 2016'


def get_band_center(form):
    comp = Composition(str(form))
    prod = 1.0
    for el, amt in comp.get_el_amt_dict().iteritems():
        prod = prod * (Element(el).X ** amt)
        return -prod ** (1 / sum(comp.get_el_amt_dict().values()))


if __name__ == '__main__':
    mpr = MPRester()
    f_w = open('Best_n_p_TEs_out.csv', 'wb')
    csv_outfile = csv.writer(f_w)
    csv_outfile.writerow(['formula', 'band gap', '(band gap)/2', 'CBM', 'band center', 'VBM', 'materials_id', 'icsd_id',
                          'Umutlist_type'])
    colnames = ['p-type', 'n-type']
    df = pandas.read_csv('Best_n_p_TEs.csv', names=colnames)
    input_data_ptype = df['p-type'].tolist()[1:]
    y = 0
    compositions = []
    for column in df:
        if column == 'p-type':
            te_type = 'p-type'
        elif column == 'n-type':
            te_type = 'n-type'
        for compound in df[column].tolist()[1:]:
            if not pandas.isnull(compound) and Composition(compound) not in compositions:
                c = Composition(compound)
                compositions.append(c)
                center = get_band_center(c)
                comp_results = mpr.query(c.reduced_formula,
                                         properties=['band_gap', 'material_id', 'icsd_id', 'e_above_hull'])
                try:
                    i = comp_results[0]
                    min_e_above_hull = i['e_above_hull']
                except:
                    y += 1
                    print 'No compounds exists in MP for ' + c.reduced_formula
                    continue
                for x in range(1, len(comp_results)):
                    if comp_results[x]['e_above_hull'] < min_e_above_hull:
                        i = comp_results[x]
                        min_e_above_hull = comp_results[x]['e_above_hull']
                print ','.join([str(x) for x in
                                c.formula, c.reduced_formula, center + i['band_gap'] / 2, center, center, center - i[
                                    'band_gap'] / 2, i['material_id'], i['icsd_id'], column])
                csv_outfile.writerow(
                    [c.reduced_formula, i['band_gap'], i['band_gap'] / 2, center + i['band_gap'] / 2, center,
                     center - i['band_gap'] / 2, i['material_id'], i['icsd_id'], column])
    print 'DONE! Total number of compounds with no data found for in MP = {}'.format(y)
    f_w.close()
