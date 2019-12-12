import unittest
import os
import logging
import filecmp
from albatradis.ExperimentCollection import ExperimentCollection
from albatradis.ExperimentMetaData import ExperimentMetaData

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

data_dir = os.path.join('data','experimentcollection')

class TestExperimentCollection(unittest.TestCase):

	def test_vary_mic(self):
		e = []
		e.append(ExperimentMetaData('drug','X','X','X',1,    0, os.path.join(data_dir, 'condition_1_0_rep1'),    os.path.join(data_dir, 'condition_1_0_rep2')))
		e.append(ExperimentMetaData('drug','X','X','X',0.5,  0, os.path.join(data_dir, 'condition_0.5_0_rep1'),  os.path.join(data_dir, 'condition_0.5_0_rep2')))
		e.append(ExperimentMetaData('drug','X','X','X',2,    0, os.path.join(data_dir, 'condition_2_0_rep1'),    os.path.join(data_dir, 'condition_2_0_rep2')))
		e.append(ExperimentMetaData('drug','X','X','X',0.25, 0, os.path.join(data_dir, 'condition_0.25_0_rep1'), os.path.join(data_dir, 'condition_0.25_0_rep2')))
		ec = ExperimentCollection(e, [os.path.join(data_dir, 'control')], os.path.join(data_dir, 'reference'))
		self.assertEqual({'drug': 'drug','induction': 0}, ec.properties_in_common())
		self.assertEqual(['mic'], ec.properties_not_in_common())
		self.assertEqual('drug_0', ec.project_name)
		
		expected_output = """project.drug_0.sequence=data/experimentcollection/reference
project.drug_0.userplot=data/experimentcollection/control data/experimentcollection/condition_0.25_0_rep1 data/experimentcollection/condition_0.25_0_rep2 data/experimentcollection/condition_0.5_0_rep1 data/experimentcollection/condition_0.5_0_rep2 data/experimentcollection/condition_1_0_rep1 data/experimentcollection/condition_1_0_rep2 data/experimentcollection/condition_2_0_rep1 data/experimentcollection/condition_2_0_rep2
"""
		self.assertEqual(expected_output, str(ec))
		
