'''Driver class'''
import logging
import os
import sys

from tradistron.PresenceAbsence import PresenceAbsence
from tradistron.GeneReportSets import GeneReportSets

class TradisTronPresenceAbsence:
	def __init__(self, options):
		self.logger            = logging.getLogger(__name__)
		self.genereports       = options.genereports
		self.verbose           = options.verbose
		self.emblfile          = options.emblfile
		self.prefix            = options.prefix

		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)
			
		if not os.path.exists(self.prefix ):
			os.makedirs(self.prefix )

	def run(self):
		p = PresenceAbsence(self.genereports, self.emblfile, self.verbose, self.prefix )
		p.create_output_files()
		
		p = GeneReportSets(self.genereports, self.prefix )
		p.write_union_file()
		p.write_intersection_file()
		return self
