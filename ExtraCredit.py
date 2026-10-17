from sequenceAnalysis import FastAreader, NucParams

def getOrganismDataForComparing(fileName=None):
    ''' Runs the exact same code as main() but instead of printing desired output, returns data for codon usage, gc content, and amino acid 
    count needed to compare between organisms. '''
    
    # initialize a NucParams object using for input the information in a fastA file formatted for NucParams using FastAreader class
    myReader = FastAreader(fileName) 
    myNuc = NucParams()
    for head, seq in myReader.readFasta() :
        myNuc.addSequence(seq)

    # compute necessary class object's attribute values
    nucleotideComposition = myNuc.nucComposition()
    codonComposition = myNuc.codonComposition()
    aminoAcidComposition = myNuc.aaComposition()
    sequenceLength = myNuc.nucCount()

    # compute the GC content percentage for final output
    gcContentPercent = ( nucleotideComposition["G"] + nucleotideComposition["C"] ) * 100 / sequenceLength

    # create a copy of NucParams.rnaCodonTable, sorted alphabetically by Amino Acid
    sortItems = sorted(NucParams.rnaCodonTable.items(), key=lambda item: (item[1], item[0]))
    sortedRnaCodonTable = dict(sortItems)
    
    codonUsageList = dict() # for saving codon use percents as list of floats for comparing organisms later
    for codon, aa in sortedRnaCodonTable.items(): # calculate relative codon usage for each codon
        codonUsageList[codon] = codonComposition[codon]*100/aminoAcidComposition[aa]

    return gcContentPercent, dict(sorted(codonUsageList.items())), dict(sorted(aminoAcidComposition.items())) # return consistently sorted data needed to compare organisms
        

def compareOrganisms(organism1Data, organism2Data):
    """ Compares both organism 1 and 2's gc content, codon usage, and amino acid composition. """

    # values returned by getOrganismDataForComparing() were tuple, separates that tuple's elements into their own names for simplification
    organism1GC = organism1Data[0]
    organism1CodonCount = organism1Data[1]
    organism1AaCount = organism1Data[2]

    organism2GC = organism2Data[0]
    organism2CodonCount = organism2Data[1]
    organism2AaCount = organism2Data[2]
    
    # compute  differences
    GC_difference = abs(organism1GC - organism2GC)
    codonDifferenceDict = dict()
    for codon in organism1CodonCount.keys():
        codonDifferenceDict[codon] = abs(organism1CodonCount[codon] - organism2CodonCount[codon])
    aaDifferenceDict = dict()
    for aa in organism1AaCount.keys():
        aaDifferenceDict[aa] = abs(organism1AaCount[aa] - organism2AaCount[aa])

    # return differences computed
    return GC_difference, dict(sorted(codonDifferenceDict.items())), dict(sorted(aaDifferenceDict.items()))

if __name__ == "__main__":
    """ Prints final desired output of differences, creates a few variables to achieve this. """

    # again, get the tuple of data returned by compareOrganisms inot separate names for simplification
    differences = compareOrganisms(getOrganismDataForComparing('testGenome.fa'), getOrganismDataForComparing('haloVolc1_1-genes.fa'))
    GC_difference = differences[0]
    codonDifference = dict(sorted(differences[1].items()))
    aaDifference = dict(sorted(differences[2].items()))
    
    # use the newly made variables to output final information
    print("organisms compared: testGenome.fa, haloVolc1_1-genes.fa")
    print()
    print(f"GC difference = {GC_difference} %")
    print()
    for codon, usage in codonDifference.items():
        print(f"{codon} usage difference = {usage} %")
    print()
    for aa, count in codonDifference.items():
        print(f"{aa} usage difference = {count} %")