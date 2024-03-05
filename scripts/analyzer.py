import collections
import os

from typing import List, Tuple
from Bio import SeqIO


class GenomeVariationAnalyzer:
    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
        # TODO: sequence parser with regex and model for record
        # TODO: .upper() during parsing
        self.records = list(SeqIO.parse(fasta_file, "fasta"))
        self.reference_record = self.records.pop(0)

    def check_input_file(self) -> None:
        """
        Checks if the input file is a valida fasta file
        :return: True if valid, False otherwise
        """
        # Check if file exists
        if not os.path.exists(self.fasta_file):
            raise ValueError(f"Error: The file {self.fasta_file} does not exists.")
        # Check if file is not empty
        if os.path.getsize(self.fasta_file) == 0:
            raise ValueError(f"Error: The file {self.fasta_file} is empty.")


    def validate_sequences(self) -> None:
        """
        Checks if the sequences are valid, by their length and their bases.
        """
        reference_length = len(self.reference_record.seq)

        for record in self.records[1:]:
            # Check length
            if len(record.seq) != reference_length:
                raise ValueError(f"Sequence {record.id} has a different length compared to the reference sequence.")
            # Check bases
            if set(record.seq.upper()) - set('ACGTN-'):
                raise ValueError(f"Sequence {record.id} contains invalid bases.")


    def identify_variations(self, sequence) -> List[Tuple[int, str, str]]:
        """
        Identify variations of the given sequence compared to the reference sequence
        :param sequence: String representing the sequence to be compared to the reference
        :return: List of variations of the given sequence. Each variation is a triple containing:
        position, type of variation and the base variation
        """
        # Skipping first and last gap
        first_gap = 0
        while sequence[first_gap] == '-':
            first_gap += 1
        last_gap = len(sequence) - 1
        while sequence[last_gap] == '-':
            last_gap -= 1

        variations = []
        for pos, (ref_base, seq_base) in enumerate(zip(self.reference_record.seq[first_gap:last_gap + 1], 
                                                       sequence[first_gap:last_gap + 1]), 
                                                    start=first_gap+1):
            if ref_base != seq_base:
                # Deletion
                if seq_base == '-':
                    variations.append((pos, 'deletion', ref_base))
                # Insertion
                elif ref_base == '-':
                    variations.append((pos, 'insertion', seq_base))
                # Substitutions, Ignore the ones with base N
                elif ref_base != 'N' and seq_base != 'N':
                    variations.append((pos, 'substitution', f"{ref_base}->{seq_base}"))

        return variations

    def compare_sequences(self, record, output_file) -> Tuple[int, int]:
        """
        Compare given record to the reference record and writes the report to file
        :param record: record to be compared
        :param output_file: file to save the report to
        :return: Tuple containing id of the compared record and the number of variations
        """
        compare_id = record.id
        compare_seq = str(record.seq).upper()
        result = f"\nComparing {self.reference_record.id} with {compare_id}:\n"
        variations = self.identify_variations(compare_seq)
        
        output_file.write(result)
        if variations:
            for var in variations:
                output_file.write(f"Position: {var[0]}: {var[1]}: {var[2]}\n")
        else:
            output_file.write("No variations found.\n")

        return compare_id, len(variations)

    def generate_final_report(self, 
                              output_comparisons_file="comparisons.txt",
                              output_final_report_file="final_report.txt"):
        
        # Doing needed checks to avoid errors during execution
        self.check_input_file()
        self.validate_sequences()

        print("Generating final report...")
        genome_variations_count = {}
        all_variations = collections.defaultdict(list)
        reference_positions = []

        # Writing comparison file
        with open(output_comparisons_file, "w") as f:
            for record in self.records:
                compare_id, variations_count = self.compare_sequences(record, f)
                genome_variations_count[compare_id] = variations_count
                # TODO: Remove .upper() when model parsing with regex is implemented
                for var in self.identify_variations(str(record.seq).upper()):
                    all_variations[var].append(compare_id)

        # Writing final report file
        with open(output_final_report_file, "w") as f:
            f.write("Final Report:\n")
            max_variations_genome = max(genome_variations_count)
            min_variations_genome = min(genome_variations_count)
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
                reference_positions.append((pos, self.reference_record.seq[pos - 1]))
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
