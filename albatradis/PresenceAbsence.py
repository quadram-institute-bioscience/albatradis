import os

import dendropy
import matplotlib
import numpy
from dendropy.utility.textprocessing import StringIO
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

matplotlib.use('agg')
import matplotlib.pyplot as plt
from graphviz import Digraph

from albatradis.EMBLReader import EMBLReader
from albatradis.GeneReport import GeneReport
from albatradis.HeatMap import HeatMap
from albatradis.ReorderGenes import ReorderGenes


def saneLogFC(gene):
    return gene.max_logfc if gene is not None else 0.0


class PresenceAbsence:
    def __init__(self, genereports, emblfile, verbose, prefix):
        self.genereports = genereports
        self.verbose = verbose
        self.emblfile = emblfile
        self.prefix = prefix

        self.run()

    def run(self):
        self.features = EMBLReader(self.emblfile).features
        self.gene_names = self.generate_gene_names()
        self.reports_to_genes = self.create_reports_to_gene_list(self.gene_names)

        self.filtered_gene_names = ReorderGenes(self.gene_names, self.genereports, self.reports_to_genes).reorder_genes()
        self.filtered_reports_to_genes = self.create_reports_to_gene_list(self.filtered_gene_names)

    def generate_gene_names(self):
        return [f.qualifiers["gene"][0] for f in self.features if "gene" in f.qualifiers]

    def create_reports_to_gene_list(self, gene_names):
        reports_to_gl = {}
        for g in self.genereports:
            reports_to_gl[g] = GeneReport(g).filtered_genes(gene_names)
        return reports_to_gl

    def create_output_files(self):
        self.create_gene_logfc_spreadsheet(os.path.join(self.prefix, 'all_logfc.csv'), self.gene_names,
                                           self.reports_to_genes)
        self.create_gene_logfc_spreadsheet(os.path.join(self.prefix, 'filtered_logfc.csv'), self.filtered_gene_names,
                                           self.filtered_reports_to_genes)

        self.create_dot_graph_genes(os.path.join(self.prefix, 'logfc.dot'), self.filtered_gene_names,
                                    self.filtered_reports_to_genes)

        self.create_network(os.path.join(self.prefix, 'network.csv'), self.filtered_gene_names,
                            self.filtered_reports_to_genes)

        self.plot_distance_matrix(os.path.join(self.prefix, 'distance_matrix_dendrogram.tre'))
        self.create_nj_newick(os.path.join(self.prefix, 'nj_newick_tree.tre'))

        HeatMap(self.reports_to_genes, self.gene_names,
                os.path.join(self.prefix, 'full_heatmap.png')).create_heat_map()
        HeatMap(self.filtered_reports_to_genes, self.filtered_gene_names,
                os.path.join(self.prefix, 'filtered_heatmap.png')).create_heat_map()

        return self

    def create_dot_graph_genes(self, filename, gene_names, reports_to_genes):
        dot = Digraph(comment='Gene Graph')

        # Antibiotics
        for g in self.genereports:
            dot.node(g, g)

        # gene names
        for g in gene_names:
            dot.node(g, g)

        for anti in self.genereports:
            for i, gene in enumerate(reports_to_genes[anti]):
                if numpy.absolute(saneLogFC(gene)) > 0:
                    dot.edge(anti, gene_names[i])

        with open(filename, 'w') as fh:
            fh.write(dot.source)

        return self

    def create_network(self, filename, gene_names, reports_to_genes):

        header = ",".join(
            ["Source", "Target", "Type", "Id", "Label", "timeset", "Weight", "concentration", "antibiotic",
             "qvals", "logFC"])

        # Antibiotics
        # for g in self.genereports:
        #    dot.node(g, g)

        # gene names
        # for g in gene_names:
        #    dot.node(g, g)

        edges = []

        for anti in self.genereports:
            for i, gene in enumerate(reports_to_genes[anti]):
                if numpy.absolute(saneLogFC(gene)) > 0.0:
                    edges.append([anti, gene_names[i], "Directed", "0", "NA", "NA", "1", "?", "?", str(gene.min_qvalue),
                                  str(gene.max_logfc)])

        with open(filename, 'w') as fh:
            fh.write(header + "\n")
            for e in edges:
                fh.write(",".join(e) + "\n")

        return self

    def create_gene_logfc_spreadsheet(self, filename, gene_names, reports_to_genes):
        with open(filename, 'w') as fh:
            # Header
            fh.write("\t".join(['Sample'] + gene_names) + "\n")

            # Body
            for report in self.genereports:
                logfcs = [str(saneLogFC(g)) for g in reports_to_genes[report]]
                fh.write("\t".join([report] + logfcs) + "\n")

        return self


    def filter_genes_with_no_changes(self):
        genes_with_changes = []

        gene_to_freq = {}
        for i, g in enumerate(self.gene_names):
            for report_file in self.genereports:
                cell_logfc = numpy.absolute(saneLogFC(self.reports_to_genes[report_file][i]))
                if cell_logfc > 0.0:
                    if g in gene_to_freq:
                        gene_to_freq[g] += 1
                    else:
                        gene_to_freq[g] = 1

        # sort by value descending, then by gene name asc
        sorted_genes = sorted(gene_to_freq.items(), key=lambda x: (-x[1], x[0]))
        gene_name_index = {g: i for i, g in enumerate(self.gene_names)}
        sorted_gene_index = []
        for s in sorted_genes:
            ordered_gene = s[0]
            sorted_gene_index.append(gene_name_index[ordered_gene])

        for i in sorted_gene_index:
            g = self.gene_names[i]
            for report_file in self.genereports:
                if numpy.absolute(saneLogFC(self.reports_to_genes[report_file][i])) > 0.0:
                    genes_with_changes.append(g)
                    break
        return genes_with_changes


    # assumption is that any changes (+-) on a gene is counted the same.
    def pair_wise_distance(self, file_a, file_b):
        distance = 0
        for i in range(len(self.filtered_gene_names)):
            a_abs = numpy.absolute(saneLogFC(self.filtered_reports_to_genes[file_a][i]))
            b_abs = numpy.absolute(saneLogFC(self.filtered_reports_to_genes[file_b][i]))

            if (a_abs == 0.0 and b_abs == 0.0) or (a_abs > 0.0 and b_abs > 0.0):
                continue
            else:
                distance += 1
        return distance


    def distance_matrix(self):
        distances = []
        for a in self.genereports:
            dist_row = []
            for b in self.genereports:
                dist_row.append(self.pair_wise_distance(a, b))
            distances.append(dist_row)
        return distances


    def plot_distance_matrix(self, outputfile):
        mat = numpy.array(self.distance_matrix())
        dists = squareform(mat)
        linkage_matrix = linkage(dists, "single")
        plt.title("Distance matrix")
        tree = hierarchy.to_tree(linkage_matrix, False)
        ntree = self.get_newick(tree, "", tree.dist, self.genereports)
        with open(outputfile, 'w') as fh:
            fh.write(ntree)


    def nj_distance_matrix_str(self):
        output = "\t".join(['.'] + self.genereports) + "\n"

        dm = self.distance_matrix()
        for i, g in enumerate(self.genereports):
            output += "\t".join([g] + [str(d) for d in dm[i]]) + "\n"

        return output


    def create_nj_newick(self, outputfile):
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
            src=StringIO(self.nj_distance_matrix_str()),
            delimiter="\t")
        nj_tree = pdm.nj_tree()

        with open(outputfile, 'w') as fh:
            fh.write(nj_tree.as_string("newick"))

        return self


    def get_newick(self, node, newick, parentdist, leaf_names):
        if node.is_leaf():
            return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = "):%.2f%s" % (parentdist - node.dist, newick)
            else:
                newick = ");"
            newick = self.get_newick(node.get_left(), newick, node.dist, leaf_names)
            newick = self.get_newick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
            newick = "(%s" % (newick)
            return newick
