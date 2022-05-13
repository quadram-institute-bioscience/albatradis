from albatradis.Gene import Gene
import unittest


class TestGene(unittest.TestCase):
    def test_parse(self):
        line = "caiE	knockout	34780	35371	1.40823784330935	increased_insertions	reverse	rpsT__3prime	1.7376552909545698e-23	3.6220055890654397e-23"

        gene = Gene().parse_line(line)

        self.assertEqual(gene.gene_name, "caiE")
        self.assertEqual(gene.categories, ["knockout"])
        self.assertEqual(gene.feature.location.start, 34780)
        self.assertEqual(gene.feature.location.end, 35371)
        self.assertEqual(gene.max_logfc, 1.40823784330935)
        self.assertEqual(gene.min_pvalue, 1.7376552909545698e-23)
        self.assertEqual(gene.min_qvalue, 3.6220055890654397e-23)
