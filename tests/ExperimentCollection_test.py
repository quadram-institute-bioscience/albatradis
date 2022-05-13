import os
import unittest

from albatradis.ExperimentCollection import ExperimentCollection
from albatradis.ExperimentMetaData import ExperimentMetaData


class ErrorReadingFile(Exception): pass


class InvalidFileFormat(Exception): pass


test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data', 'experimentcollection')


class TestExperimentCollection(unittest.TestCase):

    def test_vary_mic(self):
        e = [ExperimentMetaData('drug', 'X', 'X', 'X', 1, 0, os.path.join(data_dir, 'condition_1_0_rep1'),
                                os.path.join(data_dir, 'condition_1_0_rep2')),
             ExperimentMetaData('drug', 'X', 'X', 'X', 0.5, 0, os.path.join(data_dir, 'condition_0.5_0_rep1'),
                                os.path.join(data_dir, 'condition_0.5_0_rep2')),
             ExperimentMetaData('drug', 'X', 'X', 'X', 2, 0, os.path.join(data_dir, 'condition_2_0_rep1'),
                                os.path.join(data_dir, 'condition_2_0_rep2')),
             ExperimentMetaData('drug', 'X', 'X', 'X', 0.25, 0, os.path.join(data_dir, 'condition_0.25_0_rep1'),
                                os.path.join(data_dir, 'condition_0.25_0_rep2'))]
        ec = ExperimentCollection(e, [os.path.join(data_dir, 'control')], os.path.join(data_dir, 'reference'))
        self.assertEqual({'drug': 'drug', 'induction': 0}, ec.properties_in_common())
        self.assertEqual(['mic'], ec.properties_not_in_common())
        self.assertEqual('drug_0', ec.project_name)
