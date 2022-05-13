import os
import unittest

from albatradis.BlockIdentifier import BlockIdentifier

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data', 'blockidentifier')

def create_block_identifier(basename):
    return BlockIdentifier(
        combined_mask_file=os.path.join(data_dir, basename),
        forward_mask_file=os.path.join(data_dir, basename + '_forward'),
        reverse_mask_file=os.path.join(data_dir, basename + '_reverse'),
        combined_score_file=os.path.join(data_dir, basename + '_scores'),
        window_size=1)

class TestBlockIdentifier(unittest.TestCase):

    def test_block_identifier_normal_case(self):
        b = create_block_identifier('normalcase')
        blocks = b.block_generator()

        self.assertEqual(2, len(blocks))

        self.assertEqual([6, 14], [x.start for x in blocks])
        self.assertEqual([8, 19], [x.end for x in blocks])
        self.assertEqual([2, -3], [x.max_logfc for x in blocks])
        self.assertEqual([3, 6], [x.block_length for x in blocks])
        self.assertEqual(['forward', 'reverse'], [x.direction for x in blocks])
        self.assertEqual(['increased_insertions', 'decreased_insertions'], [x.expression for x in blocks])

    def test_noblocks(self):
        b = create_block_identifier('noblocks')
        blocks = b.block_generator()
        self.assertEqual(0, len(blocks))

    def test_blocks_at_end(self):
        b = create_block_identifier('blocksatends')
        blocks = b.block_generator()
        self.assertEqual(2, len(blocks))

        self.assertEqual([1, 16], [x.start for x in blocks])
        self.assertEqual([3, 21], [x.end for x in blocks])
        self.assertEqual([2, -3], [x.max_logfc for x in blocks])
        self.assertEqual([3, 6], [x.block_length for x in blocks])
        self.assertEqual(['forward', 'reverse'], [x.direction for x in blocks])
        self.assertEqual(['increased_insertions', 'decreased_insertions'], [x.expression for x in blocks])

    def test_peak_from_array(self):
        b = create_block_identifier('')
        self.assertEqual(15, b.peak_from_array([1, 5, 10, 15, 0, 1]))
        self.assertEqual(-15, b.peak_from_array([-1, -5, -10, -15, 0, -1]))
        self.assertEqual(-15, b.peak_from_array([1, -5, 10, -15, 0, 1]))
