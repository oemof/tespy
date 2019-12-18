# -*- coding: utf-8

from nose.tools import eq_

from tespy.tools.characteristics import (char_line, compressor_map,
                                         load_default_char, load_custom_char)

from tespy.tools.helpers import extend_basic_path
from pkg_resources import resource_filename
import json
import numpy as np
import shutil


class characteristics_tests:

    def setup(self):
        # create data path and write json files into path
        self.path = extend_basic_path('data')

    def test_custom_char_line_import(self):

        # we need to write some data to the path first, using defaults
        data_path = resource_filename('tespy.data', 'char_lines.json')

        with open(path) as f:
            raw_data = json.loads(f.read())

        data = raw_data['heat exchanger']['kA_char1']
        with open(os.path.join(self.path, 'char_maps.json'), 'w') as outfile:
            json.dump(data, outfile)

        char_original = load_default_char('heat exchanger', 'kA_char1',
                                          'EVAPORATING FLUID', char_line)
        char_custom = load_custom_char('EVAPORATING FLUID', char_line)

        x_cond = np.array_equal(char_original.x, char_custom.x)
        y_cond = np.array_equal(char_original.y, char_custom.y)

        msg = ('The x values from the custom characteristic line ' +
               str(char_custom.x) + ' must be identical to the x values from '
               'the default characteristic line ' + str(char_original.x) + ' '
               'as these have been duplicated before load.')
        eq_(True, x_cond, msg)

        msg = ('The y values from the custom characteristic line ' +
               str(char_custom.y) + ' must be identical to the y values from '
               'the default characteristic line ' + str(char_original.y) + ' '
               'as these have been duplicated before load.')
        eq_(True, y_cond, msg)
