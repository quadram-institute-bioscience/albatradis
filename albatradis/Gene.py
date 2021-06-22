import numpy


class Gene:
    def __init__(self, feature, blocks):
        self.feature = feature
        self.blocks = blocks
        self.categories = []
        self.upstream = []
        self.five_prime = None
        self.three_prime = None
        self.max_logfc = 0.0
        self.min_pvalue = 1.0
        self.min_qvalue = 1.0
        self.gene_name = self.calc_gene_name()

    def calc_gene_name(self):
        if self.feature:
            gene_name_val = str(self.feature.location.start) + "_" + str(self.feature.location.end)
            if "gene" in self.feature.qualifiers:
                gene_name_val = self.feature.qualifiers["gene"][0]
                return gene_name_val
        else:
            return "unknown"

    def category(self):
        if 'knockout' in self.categories:
            return 'knockout'
        elif 'upregulated' in self.categories:
            return 'upregulated'
        elif 'downregulated' in self.categories:
            return 'downregulated'
        else:
            return "/".join(list(set(self.categories)))

    def upstream_gene(self):
        return "/".join(list(set(self.upstream)))

    def max_logfc_from_blocks(self):
        if self.blocks:
            all_logfc = [b.max_logfc for b in self.blocks]
            highest_logfc = numpy.max(numpy.absolute(all_logfc))
            for a in all_logfc:
                if a == highest_logfc:
                    return highest_logfc
                elif a == highest_logfc * -1.0:
                    return highest_logfc * -1.0
        else:
            return 0

    def min_pvalue_from_blocks(self):
        if self.blocks:
            all_pvals = [b.min_pvalue for b in self.blocks]
            return numpy.min(all_pvals)
        else:
            return 0.0

    def min_qvalue_from_blocks(self):
        if self.blocks:
            all_qvals = [b.min_qvalue for b in self.blocks]
            return numpy.min(all_qvals)
        else:
            return 0.0

    def max_logfc_from_category(self):
        l = self.max_logfc_from_blocks()

        if self.category() == 'upregulated':
            if l < 0:
                return l * -1
            else:
                return l
        elif self.category() == 'downregulated':
            if l > 0:
                return l * -1
            else:
                return l

        if self.expression_from_blocks() == 'increased_insertions':
            if l < 0:
                return l * -1
            else:
                return l
        elif self.expression_from_blocks() == 'decreased_insertions':
            if l > 0:
                return l * -1
            else:
                return l

        return l

    def expression_from_blocks(self):
        if self.blocks:
            return "/".join(list(set([b.expression for b in self.blocks])))
        else:
            return ""

    def direction_from_blocks(self):
        if self.blocks:
            return "/".join(list(set([b.direction for b in self.blocks])))
        else:
            return ""

    def update_feature(self):
        self.feature = """FT   CDS             {window_start}..{window_end}
    		FT                   /gene="{gene_name}"
    		FT                   /locus_tag="{gene_name}"
    		FT                   /product=product
    		""".format(gene_name=self.gene_name, window_start=str(self.start), window_end=str(self.end))

    def window_string(self):
        return "\t".join(
            [str(self.start) + '_' + str(self.end), str(self.category()), str(self.start), str(self.end),
             str(self.max_logfc_from_category()), str(self.expression_from_blocks()), str(self.direction_from_blocks()),
             str(self.upstream_gene()), str(self.min_pvalue_from_blocks()), str(self.min_qvalue_from_blocks())])

    def __str__(self):
        try:
            teststring = "\t".join([str(self.gene_name), str(self.category()), str(self.feature.location.start),
                                    str(self.feature.location.end), str(self.max_logfc_from_category()),
                                    str(self.expression_from_blocks()), str(self.direction_from_blocks()),
                                    str(self.upstream_gene()), str(self.min_pvalue_from_blocks()),
                                    str(self.min_qvalue_from_blocks())])
            return teststring
        except:
            raise ValueError("some problem")

    def header(self):
        return "\t".join(['Gene', 'Category', 'Start', 'End', 'MaxLogFC', 'Expression', 'Direction', 'Upstream', 'P-Value', 'Q-Value'])
