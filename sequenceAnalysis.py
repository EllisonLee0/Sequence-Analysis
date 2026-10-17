import sys

class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

class NucParams:
    """ NucParams objects are instantiated by taking some string of text as input. From there, methods can be called on a NucParams to
    determine its DNA sequence, sequence length, count of each nucleotide, count of each codon, and count of each amino acid. Currently only 
    takes DNA sequences as input, doesn't work for RNA sequences. Note: if you use addSequence, you have to call all other methods on the 
    object to update their values. Also, to correctly work, the methods need to be called in this order each time: 1. nucComposition, 
    2. nucCount, 3. codonComposition, 4. aaComposition. """

    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        ''' Initializes attributes that will store mRNA sequence, nucleotide counts, codon counts, amino acid counts. '''
        self.seq = inString.upper().replace("T","U") # str that replaces all T's with U's and uppercases input sequence, removes spaces
        self.aaComp = { item:0 for item in set(NucParams.rnaCodonTable.values()) } # dict that stores amino acid counts (assume start at frame 1)
        self.nucComp = {"A":0, "C":0, "G":0, "U":0, "N":0} # dict that stores nucleotide counts
        self.codonComp = { key:0 for key in NucParams.rnaCodonTable.keys() } # dict that stores codon counts (assume start at frame 1)
    
    def addSequence (self, inSeq):
        ''' Updates self.seq by concatenating the value of str inSeq to the end of self.seq. '''
        inSeq = inSeq.replace("T","U")
        self.seq += inSeq
        return self.seq
      
    def aaComposition(self):
        ''' Updates key:value pairs of self.aaComp with number of times each amino acid will appear in protein sequence resulting from value of self.seq. '''
        # for each codon:count item in self.codonComp, figure out which amino acid corresponds to the codon key, then add the value at that self.codonComp codon key to the value of the corresponding amino acid key in self.aaComp
        # I wonder if there is a way to do this using less if elif statements?
        for codon in self.codonComp.keys():
            if NucParams.rnaCodonTable[codon] == "L":
                self.aaComp["L"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "R":
                self.aaComp["R"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "K":
                self.aaComp["K"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "C":
                self.aaComp["C"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "G":
                self.aaComp["G"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "W":
                self.aaComp["W"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "-":
                self.aaComp["-"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "D":
                self.aaComp["D"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "E":
                self.aaComp["E"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "Y":
                self.aaComp["Y"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "M":
                self.aaComp["M"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "Q":
                self.aaComp["Q"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "N":
                self.aaComp["N"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "I":
                self.aaComp["I"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "P":
                self.aaComp["P"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "S":
                self.aaComp["S"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "T":
                self.aaComp["T"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "F":
                self.aaComp["F"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "H":
                self.aaComp["H"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "A":
                self.aaComp["A"] += self.codonComp[codon]
            elif NucParams.rnaCodonTable[codon] == "V":
                self.aaComp["V"] += self.codonComp[codon]
        return self.aaComp
        
    def nucComposition(self):
        ''' Updates key:value pairs of self.nucComp with count of each nucleotide in self.seq. '''
        for nuc in self.seq:
            self.nucComp[nuc] += 1
        return self.nucComp
        
    def codonComposition(self):
        ''' Updates key:value pairs of self.codonComp with count of each codon in self.seq. '''
        for i in range(0, len(self.seq), 3):
            codon = self.seq[i:i+3]
            self.codonComp[codon] += 1
        return self.codonComp
      
    def nucCount(self):
        ''' Returns length of self.seq by summing counts of each nucleotide in self.nucComp. '''
        return sum(self.nucComp.values())

class ProteinParam :
    """ A ProteinParam object has only 1 attribute: proteinDict, a dictionary mapping every valid amino acid to their respective number of
    occurences in the input protein. ProteinParam has several following class variables, explained just below. Using this initial data, the
    class's methods calculate the physical-chemical properties of the input protein, and return the results. """
    
    # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini
    #  Feel free to move these to appropriate methods as you like
    
    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34
    
    def __init__ (self, protein):
        ''' ProteinParam objects have only 1 attribute: a dictionary recording count of each type of amino acid present in input protein. '''
        
        fixedProtein = protein.upper() # user input protein sequence might not be all uppercase but needs to be for this program to work, this takes care of that
        self.proteinDict = { aa:0 for aa in ProteinParam.aa2mw.keys() } # initialize attribute as dict w/ keys being same as classwide dict aa2mw all mapping to value (int) 0
        
        # update values of self.proteinDict to hold count of each amino acid in protein sequence
        for a in fixedProtein:
            if a in self.proteinDict.keys(): # make sure that no extra keys and counts for characters that aren't valid amino acids are added
                self.proteinDict[a] = fixedProtein.count(a)
                # final dict format: { "amino acid 1 letter abreviation" : int(number of times that amino acid appears in input protein sequence) }

    def aaCount (self):
        ''' Returns total number of valid amino acids there are in input protein sequence as integer. '''
        return sum(self.proteinDict.values())

    def pI (self, precision=0.01):
        ''' Calculates the theoretical pH value that will neutralize the net charge of the protein, to 2 decimal places. Uses binary search
        method, not brute force of calling _charge_() at a pH value in the middle of the search range, and adjusting the search bounds
        according to the result. '''
        
        high = 14.0 # initialize high end of search range
        low = 0.0 # initialize low end of search range
        
        while high - low > precision: # until the pI value is found to a precision of 2 decimal places (or otherwise specified precision)
            
            # calculate the net charge of the protein at the pH value in the middle of the search range
            mid = (low + high) / 2
            charge = self._charge_(mid)

            # if the net charge is positive at the middle pH, in the next loop iteration, repeat the search using the larger half of the range used in the previous iteration
            if charge > 0:
                low = mid
                
            # if the net charge is negative at the middle pH, in the next loop iteration, repeat the search using the smaller half of the range used in the previous iteration
            else:
                high = mid
                
        return (low + high) / 2 # return the pI that will eventually be found       

    def aaComposition (self) :
        ''' Return a dictionary with keys being every valid amino acid 1-letter abbreviation (str) and their associated values being the
        number of times that amino acid appears in the input protein sequence (int). If an amino acid never appears in input, key still
        exists and has value 0. '''
        return self.proteinDict # such a dictionary is in fact the exact value of self.proteinDict, so just return that
    
    def _charge_ (self, pH):
        ''' Calculates and returns the net charge of input protein at some given pH, used in the pI method. Formula for this is given in
        instructions. '''
        
        aaPosTerm = 0
        
        # calculate sum of charges for groups with positive default charge (K, R, H, N-terminus)
        for aa in ProteinParam.aa2chargePos.keys():
            aaPosTerm += self.proteinDict[aa] * (1 / (1 + 10**(pH - ProteinParam.aa2chargePos[aa])))
        aaPosTerm += 1 / (1 + 10**(pH - ProteinParam.aaNterm))
        aaNegTerm = 0
        
        # calculate sum of charges for groups with negative default charge (D, E, C, Y, C-terminus)
        for aa in ProteinParam.aa2chargeNeg.keys():
            aaNegTerm += self.proteinDict[aa] * (1 / (1 + 10**(ProteinParam.aa2chargeNeg[aa] - pH)))
        aaNegTerm += 1 / (1 + 10**(ProteinParam.aaCterm - pH))
        
        # Plug in both terms into the final expression for net charge in a given pH and return result
        netCharge = aaPosTerm - aaNegTerm
        return netCharge

    def molarExtinction (self):
        ''' Calculates and returns how much light input protein absorbs at wavelength 280nnm. Formula given in instructions. '''
        return self.proteinDict["Y"] * ProteinParam.aa2abs280["Y"] + self.proteinDict["W"] * ProteinParam.aa2abs280["W"] + self.proteinDict["C"] * ProteinParam.aa2abs280["C"]

    def massExtinction (self):
        ''' Returns the mass extinction coefficient of the input protein. '''
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        ''' Calculates and returns molecular weight of input protein as float. '''
        molWeightSum = 0 # will store combined weight of every amino acid in the protein
        
        # this loop calculates the correct value for molWeightSum
        for aa in self.proteinDict.keys(): # for every amino acid type present in the input protein, increment molWeightSum by the following amount:
            molWeightSum += self.proteinDict[aa] * (ProteinParam.aa2mw[aa] - ProteinParam.mwH2O) # number of this particular type of amino acid present in the protein * (molecular weight of that amino acid - molecular weight of water)
        
        return (ProteinParam.mwH2O + molWeightSum) # this is exactly equivalent to the formula given, where the first term is molecular weight of water etc