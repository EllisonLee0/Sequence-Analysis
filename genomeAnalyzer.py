from sequenceAnalysis import FastAreader, NucParams

def main(fileName=None):
    """
    Program overview:
    The genomeAnalyzer program takes a fastA file containing the genetic information of some organism as input and outputs its total genome
    sequence length, GC content percentage, and codon preference. To achieve this, classes from module sequenceAnalysis, FastAreader and 
    NucParams were imported. A FastAreader object is instantiated to format input fastA file data into a form that can then be used to
    instantiate a NucParams object. This main() function's code then makes calls to NucParams methods on the NucParams object to calculate and 
    display desired final output.

    Note:
    Of course, this program should be able to produce a similar corresponding output for any other input fastA file.

    Assumption:
    It was assumed that initial input fastA files contained DNA (ACGTN) sequences, but NucParams is coded such that even if it contained RNA 
    (ACGUN) sequences were contained in the initial input fastA file instead, the correct final output should still be produced.

    Assumption:
    It was assumed that regardless of if initial input fastA files contain DNA and/or RNA sequences, that these sequences code directly for the
    resulting amino acids. Therefore, nothing in any imported modules or code of this program translates input sequences to mRNA first.

    Assumption:
    It was also assumed that the first codon starts at reading frame 1. 
    """

    # initialize a NucParams object using for input the information in a fastA file formatted for NucParams using FastAreader class
    myReader = FastAreader(fileName) 
    myNuc = NucParams()
    for head, seq in myReader.readFasta() :
        myNuc.addSequence(seq)

    # compute necessary class object's attribute values for final output
    nucleotideComposition = myNuc.nucComposition()
    codonComposition = myNuc.codonComposition()
    aminoAcidComposition = myNuc.aaComposition()
    sequenceLength = myNuc.nucCount()

    # compute the GC content percentage for final output
    gcContentPercent = ( nucleotideComposition["G"] + nucleotideComposition["C"] ) * 100 / sequenceLength

    # create a copy of NucParams.rnaCodonTable, sorted alphabetically by Amino Acid
    sortItems = sorted(NucParams.rnaCodonTable.items(), key=lambda item: (item[1], item[0]))
    sortedRnaCodonTable = dict(sortItems)

    # print final output
    print(f"sequence length = {sequenceLength/1000000:.2f} Mb")
    print()
    print(f"GC content = {gcContentPercent:.1f}%")
    print()
    # calculate relative codon usage for each codon and print
    for codon, aa in sortedRnaCodonTable.items():
        print(f"{codon} : {aa} {codonComposition[codon]*100/aminoAcidComposition[aa]:5.1f} ({codonComposition[codon]:6d})")

if __name__ == "__main__":
    main()