import unittest
import os
import sys
import logging
from tradistron.BlockInsertions import BlockInsertions
import shutil

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','gent')

class TestGent(unittest.TestCase):
	
	#def test_gent_small_block_13800(self):
	#	self.check_single_region(13800)
	#def test_gent_small_block_121600(self):
	#	self.check_single_region(121600)
	#def test_gent_small_block_156866(self):
	#	self.check_single_region(156866)
	#def test_gent_small_block_211487(self):
	#	self.check_single_region(211487)
	#def test_gent_small_block_223182(self):
	#	self.check_single_region(223182)
	#def test_gent_small_block_445412(self):
	#	self.check_single_region(445412)
	#def test_gent_small_block_634775(self):
	#	self.check_single_region(634775)
	#def test_gent_small_block_652833(self):
	#	self.check_single_region(652833)
	#def test_gent_small_block_691444(self):
	#	self.check_single_region(691444)
	#def test_gent_small_block_769616(self):
	#	self.check_single_region(769616)
	#def test_gent_small_block_1014715(self):
	#	self.check_single_region(1014715)
	#def test_gent_small_block_1026874(self):
	#	self.check_single_region(1026874)
	#def test_gent_small_block_1162192(self):
	#	self.check_single_region(1162192)
	#def test_gent_small_block_1327574(self):
	#	self.check_single_region(1327574)
	#def test_gent_small_block_1751715(self):
	#	self.check_single_region(1751715)
	#def test_gent_small_block_2388671(self):
	#	self.check_single_region(2388671)
	#def test_gent_small_block_2529399(self):
	#	self.check_single_region(2529399)
	#def test_gent_small_block_2651320(self):
	#	self.check_single_region(2651320)
	#def test_gent_small_block_2721782(self):
	#	self.check_single_region(2721782)
	#def test_gent_small_block_2715832(self):
	#	self.check_single_region(2715832)
	#def test_gent_small_block_2747711(self):
	#	self.check_single_region(2747711)
	#def test_gent_small_block_2940889(self):
	#	self.check_single_region(2940889)
	#def test_gent_small_block_2960266(self):
	#	self.check_single_region(2960266)
	#def test_gent_small_block_3045044(self):
	#	self.check_single_region(3045044)
	#def test_gent_small_block_3047108(self):
	#	self.check_single_region(3047108)
	#def test_gent_small_block_3210373(self):
	#	self.check_single_region(3210373)
	#def test_gent_small_block_3302653(self):
	#	self.check_single_region(3302653)
	#def test_gent_small_block_3334569(self):
	#	self.check_single_region(3334569)
	#def test_gent_small_block_3339684(self):
	#	self.check_single_region(3339684)
	#def test_gent_small_block_3344576(self):
	#	self.check_single_region(3344576)
	#def test_gent_small_block_3404735(self):
	#	self.check_single_region(3404735)
	#def test_gent_small_block_3417659(self):
	#	self.check_single_region(3417659)
	#def test_gent_small_block_3539386(self):
	#	self.check_single_region(3539386)
	#def test_gent_small_block_3713438(self):
	#	self.check_single_region(3713438)
	#def test_gent_small_block_3916762(self):
	#	self.check_single_region(3916762)
	#def test_gent_small_block_3935861(self):
	#	self.check_single_region(3935861)
	#def test_gent_small_block_3976025(self):
	#	self.check_single_region(3976025)
	#def test_gent_small_block_4013397(self):
	#	self.check_single_region(4013397)
	#def test_gent_small_block_4020294(self):
	#	self.check_single_region(4020294)
	#def test_gent_small_block_4031450(self):
	#	self.check_single_region(4031450)
	#def test_gent_small_block_4094661(self):
	#	self.check_single_region(4094661)
	#def test_gent_small_block_4159281(self):
	#	self.check_single_region(4159281)
	#def test_gent_small_block_4198762(self):
	#	self.check_single_region(4198762)
	#def test_gent_small_block_4246679(self):
	#	self.check_single_region(4246679)
	#def test_gent_small_block_4532102(self):
	#	self.check_single_region(4532102)
	#def test_gent_small_block_4629721(self):
	#	self.check_single_region(4629721)
	#def test_gent_small_block_4596028(self):
	#	self.check_single_region(4596028)
	#def test_gent_small_block_4415908(self):
	#	self.check_single_region(4415908)
		
	def check_single_region(self, start_coord):
		print("Testing:\t"+ str(start_coord))
		case = os.path.join(data_dir, 'gent-5MIC.'+ str(start_coord) +'.plot.gz')
		control = os.path.join(data_dir, 'controlrep1.'+ str(start_coord) +'.plot.gz')
		emblfile = os.path.join(data_dir, 'annotation.embl')
		
		logger = logging.getLogger(__name__)
		b = BlockInsertions(logger,[case, control], 6, 100, 25, False, 2, 0.05, 'testoutput', 8, 50, 1, emblfile)
		self.assertTrue(b.run())
		self.assertTrue(os.path.exists('testoutput'))

		self.assertTrue(self.is_in_block(b.blocks,10000))
		shutil.rmtree("testoutput")
		
	def is_in_block(self,blocks,coord):
		for block in blocks:
			if block.start <= coord and block.end >= coord:
				return True
		return False
		
		