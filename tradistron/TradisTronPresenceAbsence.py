'''Driver class'''
import logging
import os
import sys

from tradistron.PresenceAbsence import PresenceAbsence


class TradisTronPresenceAbsence:
	def __init__(self, options):
		self.logger            = logging.getLogger(__name__)
		self.genereports       = options.genereports
		self.verbose           = options.verbose
		self.emblfile          = options.emblfile
		self.filter_no_data    = options.filter_no_data

		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)

	def run(self):
		PresenceAbsence(self.genereports, self.emblfile, self.filter_no_data ,self.verbose ).run()
		return self
