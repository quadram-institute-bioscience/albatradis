#!/usr/bin/env python3

import argparse
import sys
import os
import pkg_resources
sys.path.append('../')
sys.path.append('./')
from tradistron.TradisTronGeneReports import TradisTronGeneReports


version = ''
try:
	version = pkg_resources.get_distribution("tradistron").version
except pkg_resources.DistributionNotFound:
	version = 'x.y.z'

parser = argparse.ArgumentParser(
	description = 'Manipulate gene_report.csv files',
	usage = 'tradistron-gene_reports [options] gene_report1.csv gene_report2.csv ...', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('genereports', help='Gene report spreadsheets', nargs='+', type=str)
parser.add_argument('--prefix', '-o',  help='Output directory prefix', type=str, default='output')

parser.add_argument('--verbose', '-v', action='store_true', help='Print out more information about the analysis while it runs', default = False)
parser.add_argument('--debug', action='store_true', help='Turn on debugging', default = False)
parser.add_argument('--version', action='version', version=str(version))

options = parser.parse_args()

if options.debug:
	options.verbose = True
	import cProfile, pstats, io
	pr = cProfile.Profile()
	pr.enable()
		
	tradistron = TradisTronGeneReports(options)
	tradistron.run()
		
	pr.disable()
	s = io.StringIO()
	sortby = 'cumulative'
	ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
	ps.print_stats()
	print(s.getvalue())
else:
	tradistron = TradisTronGeneReports(options)
	tradistron.run()
