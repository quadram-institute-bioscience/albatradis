'''Driver class'''
import logging
import os
import sys

from albatradis.GeneReportSets import GeneReportSets

class AlbaTraDISGeneReports:
	def __init__(self, options):
		self.logger            = logging.getLogger(__name__)
		self.genereports       = options.genereports
		self.verbose           = options.verbose
		self.prefix            = options.prefix

		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)
			
		if not os.path.exists(self.prefix ):
			os.makedirs(self.prefix )

	def run(self):
		p = GeneReportSets(self.genereports, self.prefix )
		p.write_union_file()
		return self
