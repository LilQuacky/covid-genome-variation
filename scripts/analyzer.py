# analyzer.py

import collections
import os

from Bio import SeqIO


class GenomeVariationAnalyzer:
    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
        self.records = list(SeqIO.parse(fasta_file, "fasta"))
        self.reference_id = self.records[0].id
        self.reference_seq = str(self.records[0].seq).upper()

    def check_input_file(self):
        if not os.path.exists(self.fasta_file):
            print(f"Error: The file {self.fasta_file} does not exists.")
            return False
        if os.path.getsize(self.fasta_file) == 0:
            print(f"Error: The file {self.fasta_file} is empty.")
            return False
        print("The input file is valid.")
        return True

    def validate_sequences(self):
        reference_length = len(self.reference_seq)
        for record in self.records[1:]:
            if len(record.seq) != reference_length:
                raise ValueError(f"Sequence {record.id} has a different length compared to the reference sequence.")
            if set(record.seq.upper()) - set('ACGTN-'):
                raise ValueError(f"Sequence {record.id} contains invalid bases.")
        print("All sequences have been successfully validated.")

    def identify_variations(self, reference, sequence):
        variations = []
        for pos, (ref_base, seq_base) in enumerate(zip(reference, sequence), start=1):
            if ref_base != seq_base:
                if seq_base == '-':
                    variations.append((pos, 'deletion', ref_base))
                elif ref_base == '-':
                    variations.append((pos, 'insertion', seq_base))
                elif ref_base != 'N' and seq_base != 'N':
                    variations.append((pos, 'substitution', f"{ref_base}->{seq_base}"))
        return variations

    def compare_sequences(self, record, output_file):
        compare_id = record.id
        compare_seq = str(record.seq).upper()
        result = f"\nComparing {self.reference_id} with {compare_id}:\n"
        variations = self.identify_variations(self.reference_seq, compare_seq)
        with open(output_file, "a") as f:
            f.write(result)
            if variations:
                for var in variations:
                    f.write(f"Position: {var[0]}, Detail: {var[2]}, Reference Character: {self.reference_seq[var[0]-1]}\n")
            else:
                f.write("No variations found.\n")
        return compare_id, len(variations)

    def generate_final_report(self, output_comparisons_file="comparisons.txt",
                              output_final_report_file="final_report.txt"):
        print("Generating final report...")
        genome_variations_count = {}
        all_variations = collections.defaultdict(list)
        reference_positions = []

        with open(output_comparisons_file, "w") as f:
            for record in self.records[1:]:
                compare_id, variations_count = self.compare_sequences(record, output_comparisons_file)
                genome_variations_count[compare_id] = variations_count
                for var in self.identify_variations(self.reference_seq, str(record.seq).upper()):
                    all_variations[var].append(compare_id)

        with open(output_final_report_file, "w") as f:
            f.write("Final Report:\n")
            max_variations_genome = max(genome_variations_count, key=genome_variations_count.get)
            min_variations_genome = min(genome_variations_count, key=genome_variations_count.get)
            f.write(
                f"Genome with the most variations: {max_variations_genome} ({genome_variations_count[max_variations_genome]} variations)\n")
            f.write(
                f"Genome with the least variations: {min_variations_genome} ({genome_variations_count[min_variations_genome]} variations)\n")

            common_positions = [var[0] for var, genomes in all_variations.items() if
                                len(genomes) == len(self.records) - 1]
            f.write("\nPosition of the reference sequence that change in all comparisons:\n")
            f.write(f"Total positions: {len(common_positions)}\n")
            f.write("Characters at these positions in the reference sequence:\n")
            for pos in common_positions:
                reference_positions.append((pos, self.reference_seq[pos - 1]))
            f.write(f"{reference_positions}\n")

            same_variations_positions = {'substitution': collections.Counter(), 'insertion': collections.Counter(),
                                         'deletion': collections.Counter()}
            for var, genomes in all_variations.items():
                if len(set(genomes)) == 1:
                    same_variations_positions[var[1]][(var[0], var[2])] += 1  # Including reference position
            f.write("\nPosition of the reference sequence that change in the same way in all comparisons:\n")
            for var_type, counts in same_variations_positions.items():
                f.write(f"\n{var_type.capitalize()} variations:\n")
                f.write(f"Total positions: {sum(counts.values())}\n")
                for (position, detail), count in counts.items():
                    f.write(f"Position: {position}, Detail: {detail}\n")

        print("Report Generated!")
