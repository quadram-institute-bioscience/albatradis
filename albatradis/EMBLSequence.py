# https://github.com/sanger-pathogens/gff3toembl/blob/master/gff3toembl/EMBLContig.py
# ca3687ba5bb6bc81c1df9ef6fdb7a425fb5a1e26

class EMBLSequence(object):

  def __init__(self, sequence_string):
    nucleotide_counts = self.calculate_nucleotide_counts(sequence_string)
    self.header = self.format_header(nucleotide_counts)
    self.body = self.format_sequence_body(sequence_string)
    self.length = len(sequence_string)

  def format(self):
    return self.header + '\n' + self.body + '//'+ '\n'

  def calculate_nucleotide_counts(self, sequence):
    sequence = sequence.lower()
    counts = {}
    counts['a'] = sequence.count('a')
    counts['c'] = sequence.count('c')
    counts['g'] = sequence.count('g')
    counts['t'] = sequence.count('t')
    count_of_acgt = sum(counts.values())
    counts['other'] = len(sequence) - count_of_acgt
    return counts

  def format_header(self, nucleotide_counts):
    # The SQ line can exceed 80 characters
    template = "XX\nSQ   Sequence {total} BP; {a} A; {c} C; {g} G; {t} T; {other} other;"
    total_counts = sum(nucleotide_counts.values())
    nucleotide_counts['total'] = total_counts
    return template.format(**nucleotide_counts)

  def format_sequence_body(self, sequence_string):
    sequence_string = sequence_string.lower()
    lines = self.split_sequence(sequence_string)
    def format_a_line(line):
      # a line looks like:
      # (["1234567890", "12345", '', '', '', ''], 15)
      # and should look like
      # "     1234567890 12345                                                         15"
      blocks_of_sequence, end_of_line = line
      format_arguments = blocks_of_sequence + [end_of_line]
      return "     {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:>9}".format(*format_arguments)
    formatted_lines = map(format_a_line, lines)
    return '\n'.join(formatted_lines) + '\n'

  def split_line_of_sequence(self, line_of_sequence):
    # Turns "123456789012345" into ["1234567890", "12345", '', '', '', '']
    splits = []
    line_breaks = range(0, 60, 10)
    for line_break in line_breaks:
      split = line_of_sequence[line_break:line_break+10]
      splits.append(split)
    return splits

  def split_sequence(self, sequence_string):
    splits = []
    sequence_length = len(sequence_string)
    for start_of_line in range(0, sequence_length, 60):
      # might not actually be the end of the line if the line isn't long enough
      end_of_line = start_of_line + 60
      line_of_sequence = sequence_string[start_of_line:end_of_line]
      length_of_line = len(line_of_sequence)
      end_of_line = start_of_line + length_of_line # actually end of the line
      splits.append((self.split_line_of_sequence(line_of_sequence), end_of_line))
    return splits