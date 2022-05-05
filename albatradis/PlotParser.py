import numpy
import pandas


class ErrorReadingFile(Exception): pass


class InvalidFileFormat(Exception): pass


class PlotParser:

    def __init__(self, filename, minimum_threshold=0):
        self.filename = filename
        self.minimum_threshold = minimum_threshold
        self.forward = []
        self.reverse = []
        self.combined = []
        self.genome_length = 0

        self.split_lines()

        self.total_reads = sum(self.combined)
        self.total_insertions = sum([1 for a in self.combined if a > 0])

    def split_lines(self):
        insert_site_array = pandas.read_csv(self.filename, delim_whitespace=True, dtype=float, engine='c',
                                            header=None).values

        self.genome_length = len(insert_site_array)

        self.forward = insert_site_array[:, 0]
        self.reverse = insert_site_array[:, 1]

        if self.minimum_threshold != 0:
            self.forward = self.filter_column(self.forward, self.genome_length)
            self.reverse = self.filter_column(self.reverse, self.genome_length)

        self.combined = [self.forward[i] + self.reverse[i] for i in range(0, self.genome_length)]

        return self

    def filter_column(self, ins_array, genome_length):
        abs_ins_values = numpy.absolute(ins_array)
        return [ins_array[i] if abs_ins_values[i] >= self.minimum_threshold else 0 for i in range(0, genome_length)]

    def __str__(self):
        return "\t".join((self.filename, str(len(self.combined)), str(self.total_reads), str(self.total_insertions),
                          str(self.total_reads / self.total_insertions)))


class ScoreParser:
    def __init__(self, filename, minimum_threshold=0):
        self.filename = filename
        self.minimum_threshold = minimum_threshold
        self.pvals = []
        self.qvals = []
        self.genome_length = 0

        self.split_lines()

    def split_lines(self):
        insert_site_array = pandas.read_csv(self.filename, delim_whitespace=True, dtype=float, engine='c',
                                            header=None).values

        self.genome_length = len(insert_site_array)

        self.pvals = insert_site_array[:, 0]
        self.qvals = insert_site_array[:, 1]

        if self.minimum_threshold != 0:
            self.pvals = self.filter_column(self.pvals, self.genome_length)
            self.qvals = self.filter_column(self.qvals, self.genome_length)

        return self

    def filter_column(self, ins_array, genome_length):
        abs_ins_values = numpy.absolute(ins_array)
        return [ins_array[i] if abs_ins_values[i] >= self.minimum_threshold else 0 for i in range(0, genome_length)]
