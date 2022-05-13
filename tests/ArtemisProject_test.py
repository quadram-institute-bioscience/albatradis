import os
import unittest

from albatradis.ArtemisProject import ArtemisProject


class ErrorReadingFile(Exception): pass


class InvalidFileFormat(Exception): pass


test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data', 'experimentcollection')


class TestArtemisProject(unittest.TestCase):

    def test_vary_mic(self):
        outputfile = "project_file"
        ap = ArtemisProject(outputfile, False, os.path.join(data_dir, 'spreadsheet'),
                            os.path.join(data_dir, 'reference'), [os.path.join(data_dir, 'control')])

        self.assertTrue(ap.create_project_file())
        self.assertTrue(os.path.exists(outputfile))

        os.remove(outputfile)
