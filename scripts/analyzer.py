import collections
import os

from typing import List, Tuple
from Bio import SeqIO


class GenomeVariationAnalyzer:
    def __init__(self, fasta_file: str):
        self.fasta_file = fasta_file
        # TODO: sequence parser with regex and model for record
        # TODO: .upper() during parsing
        self.records = list(SeqIO.parse(fasta_file, "fasta"))
        self.reference_record = self.records.pop(0)

    def check_input_file(self) -> None:
        """
        Checks if the input file is a valida fasta file
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


    def identify_variations(self, sequence: str) -> List[Tuple[int, str, str]]:
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

        return variations, len(variations)

    def generate_final_report(self, 
                              output_comparisons_file: str = "comparisons.txt",
                              output_final_report_file: str = "final_report.txt"):
        
        # Doing needed checks to avoid errors during execution
        self.check_input_file()
        self.validate_sequences()

        print("Generating final report...")
        
        # Saving for each genome its variations
        variations = {}
        variations_counter = {}
        for r in self.records:
            variations[r.id], variations_counter[r.id] = self.identify_variations(r.seq)

        # Writing comparison file
        with open(output_comparisons_file, "w") as f:
            for r in self.records:
                f.write(f"\nComparing {self.reference_record.id} with {r.id}:\n")
                if variations[r.id]:
                    for var in variations[r.id]:
                        f.write(f"Position: {var[0]}: {var[1]}: {var[2]}\n")
                else:
                    f.write("No variations found.\n")


        # Writing final report file
        with open(output_final_report_file, "w") as f:
            f.write("Final Report:\n")

            # Writing genome with most and least variations
            max_variations_genome = max(variations_counter, key=variations_counter.get)
            min_variations_genome = min(variations_counter, key=variations_counter.get)
            f.write(
                f"Genome with the most variations: {max_variations_genome} ({variations_counter[max_variations_genome]} variations)\n")
            f.write(
                f"Genome with the least variations: {min_variations_genome} ({variations_counter[min_variations_genome]} variations)\n")



            # Writing postions that change in all records
            common_positions = list(set.intersection(*[set(map(lambda x: x[0], variations[genome])) for genome in variations.keys()]))
            common_positions.sort()
            f.write("\nPosition of the reference sequence that change in all comparisons:\n")
            f.write(f"Total positions: {len(common_positions)}\n")
            for pos in common_positions:
                f.write(f"Position: {pos}, Base: {self.reference_record.seq[pos - 1]}\n")


            # Writing postions that change in the same way in all records
            common_variations = list(set.intersection(*[set(variations[genome]) for genome in variations.keys()]))
            common_variations.sort()
            f.write("\nPosition of the reference sequence that change in the same way in all comparisons:\n")
            f.write(f"Total positions: {len(common_variations)}\n")
            for var in common_variations:
                f.write(f"Position: {var[0]}: {var[1]}: {var[2]}\n")

        print("Report Generated!")
