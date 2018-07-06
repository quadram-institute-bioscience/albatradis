import logging

from tradistron.ArtemisProject import ArtemisProject

class TradisTronArtemisProject:
	def __init__(self, options):
		self.logger               = logging.getLogger(__name__)
		self.outputfile           = options.outputfile
		self.verbose              = options.verbose
		self.experiments_metadata = options.experiments_metadata
		self.reference            = options.reference
		self.controls             = options.control

		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)
	
	def run(self):
		ap = ArtemisProject(self.outputfile, self.verbose, self.experiments_metadata, self.reference, self.controls)
		ap.create_project_file()
		return self
