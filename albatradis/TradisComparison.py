import subprocess
from tempfile import mkstemp
import os
import csv
import shutil


class GeneEssentiality:
    def __init__(self):
        self.condition = 0
        self.control = 0
        self.number_of_reps = 1

    def status(self):

        if self.condition == self.number_of_reps and self.control == 0:
            return 'conditionally essential'
        elif self.condition / self.number_of_reps < 0.5 and self.condition / self.number_of_reps > 0.5:
            return 'probably conditionally essential'
        elif self.condition == 0 and self.control == self.number_of_reps:
            return 'essential in control'
        elif self.control / self.number_of_reps > 0.5 and self.condition / self.number_of_reps < 0.5:
            return 'probably essential in control'
        elif self.condition == self.control and self.condition == self.number_of_reps:
            return 'always essential'
        elif (self.control / self.number_of_reps > 0.5) and (self.condition / self.number_of_reps > 0.5):
            return 'probably always essential'
        elif self.condition == 0 and self.control == 0:
            return 'always non-essential'
        elif (self.condition / self.number_of_reps < 0.5) and (self.control / self.number_of_reps < 0.5):
            return 'probably always non-essential'
        else:
            return 'inconsistent replicates'

class TradisComparison:
    def __init__(self, condition_files, control_files, verbose, minimum_block, only_ess_files_condition,
                 only_ess_files_control, analysis_type, prefix, exec="tradis_comparison.R"):
        self.condition_files = condition_files
        self.control_files = control_files
        self.exec = exec
        self.verbose = verbose
        self.minimum_block = minimum_block
        self.only_ess_files_condition = only_ess_files_condition
        self.only_ess_files_control = only_ess_files_control
        self.analysis_type = analysis_type
        self.prefix = prefix

        fd, self.output_filename = mkstemp()
        fd, self.conditions_fofn = mkstemp()
        fd, self.controls_fofn = mkstemp()

    def gene_names_from_essentiality_file(self, filename):

        with open(filename, 'r') as fileh:
            reader = csv.reader(fileh, delimiter=',', quotechar='"')
            gene_names = [r[1] for r in reader if r[1] != 'gene_name']
            print("Number of all genes:" + str(len(gene_names)))

        return gene_names

    def get_all_gene_names(self):

        all_gene_names = []

        for filename in self.condition_files:
            with open(filename, 'r') as fileh:
                reader = csv.reader(fileh, delimiter='\t', quotechar='"')
                gene_names1 = [r[1] for r in reader if r[1] != 'gene_name']
                all_gene_names = list(set(all_gene_names) | set(gene_names1))

        for f in self.control_files:
            with open(filename, 'r') as fileh:
                reader = csv.reader(fileh, delimiter='\t', quotechar='"')
                gene_names2 = [r[1] for r in reader if r[1] != 'gene_name']
                all_gene_names = list(set(all_gene_names) | set(gene_names2))

        return all_gene_names


    def all_gene_essentiality(self, input_filename):

        all_gene_names = self.get_all_gene_names()
        print("all_gene_names: " + str(len(all_gene_names)))
        genes_ess = {g: GeneEssentiality() for g in all_gene_names}
        if self.analysis_type == "original":
            for f in self.only_ess_files_condition:
                ess_gene_names = self.gene_names_from_essentiality_file(f)
                print("ess_gene_names condition: " + str(len(ess_gene_names)))
                print("genes_ess: " + str(len(genes_ess)))
                for e in genes_ess:
                    if e in ess_gene_names:
                        genes_ess[e].condition += 1
                    genes_ess[e].number_of_reps = len(self.only_ess_files_condition)
            for f in self.only_ess_files_control:
                ess_gene_names = self.gene_names_from_essentiality_file(f)
                print("ess_gene_names control: " + str(len(ess_gene_names)))
                print("genes_ess: " + str(len(genes_ess)))
                for e in genes_ess:
                    if e in ess_gene_names:
                        genes_ess[e].control += 1
                    genes_ess[e].number_of_reps = len(self.only_ess_files_control)
        else:
            for e in genes_ess:
                genes_ess[e].control = 0
                genes_ess[e].condition = 0
                genes_ess[e].number_of_reps = len(self.only_ess_files_control)

        return genes_ess

    def create_fofn(self):
        with open(self.conditions_fofn, 'w') as fileh:
            if len(self.condition_files) == 1:
                fileh.write(self.condition_files[0] + "\n")
                fileh.write(self.condition_files[0] + "\n")
            else:
                for i in self.condition_files:
                    fileh.write(i + "\n")

        with open(self.controls_fofn, 'w') as fileh:
            if len(self.control_files) == 1:
                fileh.write(self.control_files[0] + "\n")
                fileh.write(self.control_files[0] + "\n")
            else:
                for i in self.control_files:
                    fileh.write(i + "\n")
        return self

    def add_gene_essentiality_to_file(self, input_filename, output_filename, genes_ess):

        # We can add information on gene essentiality to the comparison output,
        # but this will not reflect full set of essential genes as the output does not contain all genes
        # In order to prevent confusion this is "switched" off, but could be used if uncommented


        with open(input_filename, 'r') as inputfh:
            output_content = []

            reader = csv.reader(inputfh, delimiter=',', quotechar='"')
            input_content = [r for r in reader]
            # if self.analysis_type == "original":
            #     print("Number of cells: " + len(input_content))
            #     for i, cells in enumerate(input_content):
            #         if i == 0:
            #             cells.append("Essentiality")
            #         elif cells[1] in genes_ess and not ("3prime" in cells[1] or "5prime" in cells[1]):
            #             cells.append(genes_ess[cells[1]].status())
            #         else:
            #             cells.append('N/A')
            #         output_content.append(cells)
            # else:
            #     output_content = input_content

            output_content = input_content

            with open(output_filename, 'w') as outputfh:
                for line in output_content:
                    outputfh.write(",".join(line) + "\n")
        return self

    def construct_command(self):
        return " ".join(
            [self.exec, '-f', '-t', str(self.minimum_block), '--controls', self.controls_fofn, '--conditions',
             self.conditions_fofn])

    def run(self):
        self.create_fofn()
        cmd = self.construct_command()
        if self.verbose:
            print(cmd)
        if self.analysis_type != "original":
            subprocess.check_output(cmd, shell=True)
        genes_ess = self.all_gene_essentiality(".output.csv")

        ess = open(os.path.join(self.prefix, "Essentiality.txt"), "w+")
        ess.write("Gene, Essentiality, Control, Condition, Replicates\n")
        for e in genes_ess:
            ess.write(e + ", " + str(genes_ess[e].status()) + ", " + str(genes_ess[e].control) + ", " + str(genes_ess[e].condition) + ", " + str(genes_ess[e].number_of_reps) + "\n")
        ess.close()

        self.add_gene_essentiality_to_file(".output.csv", self.output_filename, genes_ess)

        self.cleanup()

        return self

    def cleanup(self):
        os.remove(self.controls_fofn)
        os.remove(self.conditions_fofn)

        for f in self.only_ess_files_condition:
            os.remove(f)

        for f in self.only_ess_files_control:
            os.remove(f)
