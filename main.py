from scripts.analyzer import GenomeVariationAnalyzer

def main():
    fasta_file = "fasta/aligned-covid-sequences.fasta"
    output_comparisons_file = "report/comparisons.txt"
    output_final_report_file = "report/final_report.txt"

    analyzer = GenomeVariationAnalyzer(fasta_file)
    analyzer.check_input_file()
    analyzer.validate_sequences()
    analyzer.generate_final_report(output_comparisons_file, output_final_report_file)

if __name__ == "__main__":
    main()