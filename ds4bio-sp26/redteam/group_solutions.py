def group_1(template_sequence):
    try:
        # 1. Validation Logic: Ensure only a, g, c, t and no spaces/numbers
        valid_bases = set('agct')
        if not isinstance(template_sequence, str):
            return "ERROR"
        
        if not all(char in valid_bases for char in template_sequence.lower()) or \
           any(char.isdigit() or char.isspace() for char in template_sequence):
            return "ERROR"

        # 2. Transcription: Reverse and match with base_pairs dictionary
        base_pairs = {"g": "c", "c": "g", "a": "u", "t": "a"}
        reversed_template = template_sequence[::-1].lower()
        mrna = "".join([base_pairs[nuc] for nuc in reversed_template])

        # 3. Open Reading Frame (ORF): Must contain 'aug'
        if 'aug' not in mrna:
            return "ERROR"
        
        start_idx = mrna.find('aug')
        sliced_mrna = mrna[start_idx:]
        
        # Slicing mRNA into codons (triplets)
        codons = []
        for i in range(0, len(sliced_mrna), 3):
            chunk = sliced_mrna[i:i+3]
            # Ensure codon length is 3
            if len(chunk) == 3:
                codons.append(chunk.upper())
            else:
                break 

        # 4. Stop Codon Logic: Must find a stop codon or raise ERROR
        stop_codons = {'UAA', 'UAG', 'UGA'}
        final_codons = []
        stop_found = False
        
        for c in codons:
            final_codons.append(c)
            if c in stop_codons:
                stop_found = True
                break
        
        # NEW RULE: If no stop codon appears within codons, an error should arise
        if not stop_found:
            return "ERROR"
                
        # 5. Translation using codon_table
        codon_table = {
            'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
            'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
            'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', 'UGA': '*',
            'UGU': 'C', 'UGC': 'C', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L',
            'CUA': 'L', 'CUG': 'L', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P',
            'CCG': 'P', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AUU': 'I',
            'AUC': 'I', 'AUA': 'I', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T',
            'ACG': 'T', 'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V',
            'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'GCU': 'A', 'GCC': 'A',
            'GCA': 'A', 'GCG': 'A', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E',
            'GAG': 'E', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }

        amino_acid_array = []
        for codon in final_codons:
            # Verification: individual codon length check
            if len(codon) != 3:
                return "ERROR"
            amino_acid_array.append(codon_table[codon])

        aa_seq = "".join(amino_acid_array)
        return aa_seq

    except:
        return "ERROR"

def group_2(template_sequence):
    """
    Performs transcription and translation on a DNA template strand.
    Handles non-string inputs and invalid sequences by returning "ERROR".
    """
    
    # 1. Input Type Validation: Ensure input is a non-empty string
    if not isinstance(template_sequence, str) or not template_sequence:
        return "ERROR"

    # Define the codon table
    codon_table = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S', 'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA':'', 'UAG':'', 'UGC':'C', 'UGU':'C', 'UGA':'', 'UGG':'W'
    }

    # Internal helper: Transcribe DNA template to mRNA
    def transcribe(seq):
        base_pairs = {"g": "c", "c": "g", "a": "u", "t": "a"}
        # Process in lowercase as per lab standards
        # Template is 5'->3', transcribe reads 3'->5' to produce mRNA 5'->3'
        reversed_seq = seq.lower()[::-1]
        mrna = ""
        for nucleotide in reversed_seq:
            complement = base_pairs.get(nucleotide)
            if complement:
                mrna += complement
            else:
                # Return None if non-DNA characters are found
                return None 
        return mrna

    # Internal helper: Find the Open Reading Frame (ORF)
    def find_orf(mrna_seq):
        if mrna_seq is None: return ""
        start_codon = 'aug'
        stop_codons = ['uaa', 'uag', 'uga']
        
        start_index = mrna_seq.find(start_codon)
        if start_index == -1:
            return "" # No start codon found
        
        # Slicing from the start codon
        coding_seq = mrna_seq[start_index:]
        orf_sequence = ""
        stop_found = False
        
        # Split into triplets
        for i in range(0, len(coding_seq) - 2, 3):
            codon = coding_seq[i:i+3]
            orf_sequence += codon
            if codon in stop_codons:
                stop_found = True
                break
        
        return orf_sequence if stop_found else ""

    # Internal helper: Translate mRNA to Amino Acids
    def translate(orf_seq, table):
        if not orf_seq: return ""
        aa_seq = ""
        for i in range(0, len(orf_seq) - 2, 3):
            codon = orf_seq[i:i+3].upper()
            aa_seq += table.get(codon, "")
        return aa_seq

    # Execute workflow with safety checks
    mrna_seq = transcribe(template_sequence)
    if mrna_seq is None:
        return "ERROR" # Triggered by invalid nucleotides
        
    mrna_orf = find_orf(mrna_seq)
    aa_seq = translate(mrna_orf, codon_table)
    
    # Requirement: Return "ERROR" if no amino acid output is produced
    return aa_seq if aa_seq != "" else "ERROR"

def group_3(template):
    """
    Translates a DNA template strand to an amino acid sequence.
    
    Args:
        template: DNA template strand (string)
    
    Returns:
        Amino acid sequence (string) or "ERROR" if translation fails
    """
    # Genetic code dictionary (mRNA codon to amino acid)
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    # Input validation
    if not template or not isinstance(template, str):
        return "ERROR"
    
    # Convert to uppercase for consistency
    template = template.upper()
    
    # Check for valid DNA characters
    valid_bases = set('ATGC')
    if not all(base in valid_bases for base in template):
        return "ERROR"
    
    # Check length is multiple of 3
    if len(template) % 3 != 0:
        return "ERROR"
    
    # Convert template strand to coding strand (complement)
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    coding_strand = ''.join(complement_map[base] for base in template)
    
    # Transcribe coding strand to mRNA (T -> U)
    mrna = coding_strand.replace('T', 'U')
    
    # Check for start codon (AUG)
    if not mrna.startswith('AUG'):
        return "ERROR"
    
    # Translate mRNA to amino acids
    amino_acids = []
    for i in range(0, len(mrna), 3):
        codon = mrna[i:i+3]
        if codon in genetic_code:
            amino_acid = genetic_code[codon]
            if amino_acid == '*':  # Stop codon
                break
            amino_acids.append(amino_acid)
        else:
            return "ERROR"
    
    # Return the amino acid sequence
    return ''.join(amino_acids)

def group_4(dna_template_strand):
  
   # 1. Input Validation
   # Ensure input is a string and contains only valid DNA characters
   if not isinstance(dna_template_strand, str):
       return "ERROR"
  
   sequence = dna_template_strand.upper().strip()
   valid_bases = {'A', 'T', 'C', 'G'}
  
   if any(char not in valid_bases for char in sequence):
       return "ERROR"
  
   complement_map = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
  
   # Reverse the sequence first
   #reversed_seq = sequence[::-1]
  
   # Create mRNA by complementing
   mrna_sequence = "".join(complement_map[base] for base in sequence)
  
   # 3. Determine Open Reading Frame (ORF)
   # Find the first start codon 'AUG'
   start_codon = "AUG"
   start_index = mrna_sequence.find(start_codon)
  
   if start_index == -1:
       return "ERROR" # No ORF found
      
   # Slice the mRNA starting from the start codon
   orf_mrna = mrna_sequence[start_index:]
  
   # 4. Translation
   # Standard Genetic Code Dictionary
   codon_table = {
       'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
       'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
       'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
       'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
       'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
       'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
       'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
       'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
       'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
       'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
       'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
       'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
       'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
       'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
       'UAC': 'Y', 'UAU': 'Y', 'UAA': '_', 'UAG': '_',
       'UGC': 'C', 'UGU': 'C', 'UGA': '_', 'UGG': 'W',
   }
  
   amino_acids = []
  
   # Read codons in triplets
   for i in range(0, len(orf_mrna), 3):
       codon = orf_mrna[i:i+3]
      
       # If we have a partial codon at the end, stop
       if len(codon) < 3:
           break
          
       aa = codon_table.get(codon, '?')
      
       if aa == '_': # Stop codon
           break
          
       amino_acids.append(aa)
      
   return "".join(amino_acids)

def group_5(template_sequence):
    """
    Gemini-generated code for translate_from_template
    Translates a DNA template strand into an amino acid sequence.
    """

    # ─────────────────────────────────────────────
    # Helper function: complement a single DNA base to its RNA counterpart
    # ─────────────────────────────────────────────
    def complement_base(base):
        """
        Returns the RNA complement of a single DNA template base.
        Template base → mRNA base:  T→A, A→U, C→G, G→C
        Returns None if the base is not a valid DNA nucleotide.
        """
        complement_map = {'T': 'A', 'A': 'U', 'C': 'G', 'G': 'C'}
        return complement_map.get(base)


    # ─────────────────────────────────────────────
    # Helper function: transcribe a DNA template strand to mRNA
    # ─────────────────────────────────────────────
    def transcribe(template):
        """
        Transcribes a DNA template strand into an mRNA string by
        complementing each base (T→A, A→U, C→G, G→C).
        Returns the mRNA string, or None if any base is invalid.
        """
        mrna = ""
        for base in template:
            comp = complement_base(base)
            if comp is None:
                return None
            mrna += comp
        return mrna


    # ─────────────────────────────────────────────
    # Standard genetic code: mRNA codon → amino acid (single-letter)
    # ─────────────────────────────────────────────
    codon_table = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGG': 'W',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    stop_codons = {'UAA', 'UAG', 'UGA'}


    # ─────────────────────────────────────────────
    # Helper function: translate mRNA to protein
    # ─────────────────────────────────────────────
    def translate(mrna):
        """
        Translates an mRNA sequence to a protein (amino acid) string.
        - Scans for the first AUG start codon.
        - Reads in-frame codons (groups of 3) from the position after AUG.
        - Stops at the first in-frame stop codon (UAA, UAG, UGA).
        - Returns the amino acid sequence (the initial Met from AUG is excluded).
        - Returns "ERROR" if no start codon is found, no in-frame stop codon
        is found, or an unrecognised codon is encountered.
        """
        # Find the first AUG start codon at any position
        start_index = mrna.find('AUG')
        if start_index == -1:
            return "ERROR"

        # Begin reading codons right after the start codon
        protein = ""
        i = start_index + 3
        found_stop = False

        while i + 3 <= len(mrna):
            codon = mrna[i:i+3]

            # Check for stop codon
            if codon in stop_codons:
                found_stop = True
                break

            # Look up amino acid
            if codon in codon_table:
                protein += codon_table[codon]
            else:
                return "ERROR"

            i += 3

        if not found_stop:
            return "ERROR"

        return protein


    # ─────────────────────────────────────────────
    # Main function: translate_from_template
    # ─────────────────────────────────────────────
    def translate_from_template(template):
        """
        Takes a DNA template strand and returns the translated amino acid sequence.

        Process:
        1. Validates input is a non-empty string.
        2. Converts to uppercase (case-insensitive handling).
        3. Checks that every character is a valid DNA base (A, T, C, G).
            Spaces, numbers, and special characters cause an "ERROR".
        4. Transcribes the template strand to mRNA by complementing each base
            (T→A, A→U, C→G, G→C).
        5. Locates the first AUG start codon in the mRNA.
        6. Reads codons in the reading frame established by that AUG.
        7. Stops translation at the first in-frame stop codon (UAA, UAG, UGA).
        8. Returns the amino acid sequence as a string of single-letter codes
            (the initial Met encoded by the start codon is NOT included).
        9. Returns the string "ERROR" if translation cannot be completed for
            any reason (invalid input, missing start codon, missing stop codon).

        Parameters:
            template (str): A DNA template strand sequence.

        Returns:
            str: The amino acid sequence, or "ERROR".
        """
        # --- Input validation ---
        if not isinstance(template, str) or len(template) == 0:
            return "ERROR"

        # Case-insensitive: convert to uppercase
        template = template.upper()

        # Only A, T, C, G are valid bases (no spaces, digits, or special chars)
        valid_bases = set('ATCG')
        for base in template:
            if base not in valid_bases:
                return "ERROR"

        # --- Transcription ---
        mrna = transcribe(template)
        if mrna is None:
            return "ERROR"

        # --- Translation ---
        return translate(mrna)
    
    return translate_from_template(template_sequence)

def group_6(template_sequence):
    codon_table = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA':'', 'UAG':'',
        'UGC':'C', 'UGU':'C', 'UGA':'', 'UGG':'W'}

    def is_valid_dna(sequence):
        if not isinstance(sequence, str) or len(sequence) == 0:
            return False
        return all(base in 'ATCG' for base in sequence.upper())

    def transcribe(dna_seq):
        # Standard biological transcription: 5'-3' Template -> 5'-3' mRNA
        complement = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
        mrna_comp = "".join([complement.get(base, 'N') for base in dna_seq.upper()])
        return mrna_comp[::-1]

    def orf(mrna_seq):
        # Find ALL start codon indices (AUG)
        start_indices = [i for i in range(len(mrna_seq)) if mrna_seq[i:i+3] == 'AUG']
        
        candidates = []
        
        for start in start_indices:
            # For each start, search for a Stop Codon in-frame (steps of 3)
            for i in range(start, len(mrna_seq), 3):
                codon = mrna_seq[i:i+3]
                if codon in ['UAA', 'UAG', 'UGA']:
                    # Found a valid ORF! Save it.
                    candidates.append(mrna_seq[start:i])
                    break
                    
        # If we found any candidates, return the longest one.
        if candidates:
            return max(candidates, key=len)
            
        return None

    def translate(mrna_orf, table):
        if not mrna_orf:
            return None
        
        aa_seq = ""
        for i in range(0, len(mrna_orf), 3):
            codon = mrna_orf[i:i+3]
            if len(codon) < 3: break
            aa = table.get(codon)
            # If codon is valid and NOT a stop codon, add to sequence
            if aa and aa != '_':
                aa_seq += aa
            elif aa == '_': # If we hit a stop unexpectedly
                break
            else:
                return None
        return aa_seq

    def translate_from_template(template_sequence):
        try:
            if not is_valid_dna(template_sequence):
                return "ERROR, invalid dna"
                
            mrna_seq = transcribe(template_sequence)
            mrna_orf = orf(mrna_seq)
            
            if mrna_orf is None:
                return "ERROR, no orf"
                
            aa_seq = translate(mrna_orf, codon_table)
            
            # Final check: if aa_seq is empty string (e.g. only a start and stop), return ERROR or result?
            # Usually, a valid ORF must produce at least the Start Methionine.
            return aa_seq if aa_seq else "ERROR, not valid sequence"
            
        except Exception:
            return "ERROR, unknown"
    
    return translate_from_template(template_sequence)

def group_7(dna_sequence):
    """
    Processes a 5' to 3' DNA template strand and returns ONLY the 
    resulting amino acid sequence.
    """
    # --- STEP 0: VALIDATION LAYER ---
    if not isinstance(dna_sequence, str):
        return "ERROR"
    
    dna_sequence = dna_sequence.upper().strip()
    valid_nucleotides = set("ATCG")
    if not dna_sequence or not set(dna_sequence).issubset(valid_nucleotides):
        return "ERROR"
        
    if len(dna_sequence) < 3:
        return "ERROR"

    # --- STEP 1: REVERSE COMPLEMENT ---
    # Convert template to coding strand (reverse complement)
    complement_table = str.maketrans("ATCG", "TAGC")
    coding_strand = dna_sequence.translate(complement_table)[::-1]

    # --- STEP 2: TRANSCRIPTION ---
    mrna = coding_strand.replace('T', 'U')
    
    # --- STEP 3: ORF DETECTION ---
    start_index = mrna.find("AUG")
    if start_index == -1:
        return "ERROR"
    
    orf_mrna = mrna[start_index:]
    
    # --- STEP 4: TRANSLATION TABLE ---
    codon_map = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_', 'UGA':'_',
        'UGC':'C', 'UGU':'C', 'UGG':'W',
    }
    
    # --- STEP 5: TRANSLATION LOGIC ---
    protein = []
    has_stop_codon = False
    
    for i in range(0, len(orf_mrna) - 2, 3):
        codon = orf_mrna[i:i+3]
        amino_acid = codon_map.get(codon, '?')
        if amino_acid == '_':
            has_stop_codon = True
            break
        protein.append(amino_acid)
        
    if not has_stop_codon:
         return "ERROR"
        
    # Return ONLY the amino acid sequence string
    return "".join(protein)

def group_9(template_sequence):
    def get_dna_to_rna_map():
        """
        Returns the dictionary mapping DNA template bases to mRNA bases.
        Template -> mRNA (Complementary)
        A -> U
        T -> A
        C -> G
        G -> C
        """
        return {
            'A': 'U', 
            'T': 'A', 
            'C': 'G', 
            'G': 'C'
        }

    def get_codon_map():
        """
        Returns the standard genetic code dictionary (mRNA to Amino Acid).
        """
        return {
            # Phenylalanine
            'UUU': 'F', 'UUC': 'F',
            # Leucine
            'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
            # Isoleucine
            'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
            # Methionine (Start)
            'AUG': 'M',
            # Valine
            'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
            # Serine
            'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S',
            # Proline
            'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            # Threonine
            'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            # Alanine
            'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            # Tyrosine
            'UAU': 'Y', 'UAC': 'Y',
            # Histidine
            'CAU': 'H', 'CAC': 'H',
            # Glutamine
            'CAA': 'Q', 'CAG': 'Q',
            # Asparagine
            'AAU': 'N', 'AAC': 'N',
            # Lysine
            'AAA': 'K', 'AAG': 'K',
            # Aspartic Acid
            'GAU': 'D', 'GAC': 'D',
            # Glutamic Acid
            'GAA': 'E', 'GAG': 'E',
            # Cysteine
            'UGU': 'C', 'UGC': 'C',
            # Tryptophan
            'UGG': 'W',
            # Arginine
            'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
            # Glycine
            'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
            # Stop Codons
            'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP'
        }

    def translate_from_template(template_dna):
        """
        Takes a DNA template string, transcribes it to mRNA, finds the start codon,
        and translates it to an amino acid sequence.
        Returns "ERROR" if validation fails.
        """
        
        # 1. Validation: Check if input is empty or contains invalid characters
        if not template_dna:
            return "ERROR"
            
        dna_map = get_dna_to_rna_map()
        valid_bases = set(dna_map.keys())
        
        # Check for invalid characters (case-insensitive check handled by upper())
        template_dna = template_dna.upper()
        if any(base not in valid_bases for base in template_dna):
            return "ERROR"

        # 2. Transcription: DNA Template -> mRNA
        # We iterate through the DNA string and swap bases based on the map
        mrna = []
        for base in template_dna:
            mrna.append(dna_map[base])
        mrna_seq = "".join(mrna)

        # 3. Translation Initialization
        codon_map = get_codon_map()
        amino_acids = []
        
        # Find the start codon 'AUG'
        start_index = mrna_seq.find('AUG')
        
        # If no start codon is found, we cannot begin translation
        if start_index == -1:
            return "ERROR"

        # 4. Translation Loop
        # Start reading from the start codon index
        # We step by 3 to read triplets
        found_stop = False
        
        for i in range(start_index, len(mrna_seq), 3):
            codon = mrna_seq[i:i+3]
            
            # If we run out of bases (incomplete codon at the end), it's an error
            # unless we already hit STOP (which breaks the loop earlier)
            if len(codon) < 3:
                break

            # Get amino acid from map
            # If a weird codon appears not in map (unlikely if validation passed), return ERROR
            if codon not in codon_map:
                return "ERROR"
                
            aa = codon_map[codon]

            if aa == 'STOP':
                found_stop = True
                break
            else:
                amino_acids.append(aa)

        # 5. Final Validation
        # In strict biological translation, if we don't hit a STOP codon, 
        # the protein is often considered incomplete or the reading frame invalid.
        if not found_stop:
            return "ERROR"

        return "".join(amino_acids)
    
    return translate_from_template(template_sequence)