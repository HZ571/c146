import pandas as pd

############################## Group 1 ##############################

def group_1_count_repeats_in_string(string: str, substring: str) -> int:
    max_repeats = 0
    sub_len = len(substring)
    str_len = len(string)

    if sub_len == 0:
        return 0

    for i in range(str_len - sub_len + 1):
        count = 0
        current_pos = i
        
        while current_pos + sub_len <= str_len:
            if string[current_pos : current_pos + sub_len] == substring:
                count += 1
                current_pos += sub_len
            else:
                break
        
        if count > max_repeats:
            max_repeats = count
            
    return max_repeats

def group_1_search_database(df: pd.DataFrame, STR_counts: list) -> list:
    matches = []
    
    for name, row in df.iterrows():
        if list(row) == STR_counts:
            matches.append(name)
            
    return matches

def group_1_STR_search(database_df: pd.DataFrame, dna_seq: str) -> list:
    valid_chars = {'A', 'C', 'G', 'T'}
    dna_seq = dna_seq.upper()
    if not all(char in valid_chars for char in dna_seq):
        raise ValueError("DNA sequence contains invalid characters.")

    str_counts = []

    for str_sequence in database_df.columns:
        str_sequence = str_sequence.upper()
        
        if not all(char in valid_chars for char in str_sequence):
            raise ValueError(f"STR sequence '{str_sequence}' contains invalid characters.")
        
        count = group_1_count_repeats_in_string(dna_seq, str_sequence)
        str_counts.append(count)
        
    return group_1_search_database(database_df, str_counts)

############################## Group 2 ##############################

def group_2_count_repeats_in_string(string: str, substring: str) -> int:
    """
    Blue-Team Version: Computes the longest run of consecutive repeats.
    Optimized for robustness and efficiency without using 'in' or regex.
    
    Args:
        string (str): The sequence to search through (e.g., DNA sequence).
        substring (str): The specific sequence to count (e.g., "AATG").
        
    Returns:
        int: The maximum number of consecutive occurrences found.
    """
    # 1. Robustness: Handle empty substring or sequence shorter than STR
    sub_len = len(substring)
    str_len = len(string)
    
    if sub_len == 0 or str_len < sub_len:
        return 0

    max_repeats = 0

    # 2. Iterate through the string to check every potential starting position
    # This ensures we don't miss a run that starts at any index
    for i in range(str_len - sub_len + 1):
        count = 0
        
        # 3. Check for consecutive matches starting exactly at index i
        while True:
            # Calculate the window for the next potential repeat in the run
            start = i + (count * sub_len)
            end = start + sub_len
            
            # Boundary check: Ensure we don't slice past the end of the string
            if end <= str_len:
                # Use equality comparison (==) instead of 'in'
                if string[start:end] == substring:
                    count += 1
                else:
                    break # Run ended
            else:
                break # Not enough characters left for another repeat
        
        # Update the global maximum if this specific run was the longest found so far
        if count > max_repeats:
            max_repeats = count
            
    return max_repeats

def group_2_search_database(df: pd.DataFrame, STR_counts: list) -> list:
    """
    Blue-Team Version: Searches a pandas DataFrame of STR profiles for an exact match.
    Optimized for robustness by validating input dimensions.
    
    Args:
        df (pd.DataFrame): The database where the index contains names 
                           and columns represent different STR sequences.
        STR_counts (list): A list of integers representing the STR counts.
                           
    Returns:
        list: A list of names for all individuals who match the profile.
              Returns an empty list if no matches are found or if inputs are invalid.
    """
    # 1. Robustness: If the number of STRs provided doesn't match the database columns,
    # it's impossible to have an exact match. Returning [] prevents a Pandas ValueError.
    if len(STR_counts) != len(df.columns):
        return []

    # 2. Element-wise comparison and row-wise reduction
    # This checks every row against the list and ensures ALL columns in that row match.
    matches_mask = (df == STR_counts).all(axis=1)
    
    # 3. Filter index and convert to standard Python list
    return df.index[matches_mask].tolist()

def group_2_STR_search(database_df: pd.DataFrame, dna_seq: str) -> list:
    
    def is_valid_dna(seq: str) -> bool:
        """Helper to ensure string only contains A, T, C, G after cleaning."""
        # Check for empty sequences or non-string inputs
        if not isinstance(seq, str) or len(seq) == 0:
            return False
        valid_bases = {'A', 'T', 'C', 'G'}
        for base in seq:
            if base not in valid_bases:
                return False
        return True

    # 1. Cleaning and Normalization
    # We strip whitespace and newlines (common in text files) and uppercase everything
    dna_seq_clean = dna_seq.strip().upper()
    
    # 2. Input Validation (Defense against Red-team malformed DNA)
    if not is_valid_dna(dna_seq_clean):
        # We return an empty list or could raise an error; 
        # For a search function, [] indicates "no valid match found"
        return []

    # 3. Extract STR patterns from DataFrame columns
    str_patterns = list(database_df.columns)
    
    # 4. Compute the STR profile for the DNA sequence
    query_profile = []
    for pattern in str_patterns:
        pattern_upper = pattern.strip().upper()
        
        # Validation of database columns
        if not is_valid_dna(pattern_upper):
            continue # Skip invalid STR columns
        
        # Use our robust counting function
        count = group_2_count_repeats_in_string(dna_seq_clean, pattern_upper)
        query_profile.append(count)

    # 5. Search the database using the calculated profile
    return group_2_search_database(database_df, query_profile)

############################## Group 3 ##############################

def group_3_count_repeats_in_string(string: str, substring: str) -> int:
    if not substring or not string:
        return 0

    max_count = 0
    i = 0
    n = len(string)
    k = len(substring)

    while i <= n - k:
        if string[i:i + k] == substring:
            # Found the start of a run; count consecutive repeats
            count = 0
            while string[i:i + k] == substring:
                count += 1
                i += k
            max_count = max(max_count, count)
        else:
            i += 1

    return max_count

def group_3_search_database(df: pd.DataFrame, STR_counts: list) -> list:
    if len(STR_counts) != len(df.columns):
        raise ValueError(
            f"STR_counts length ({len(STR_counts)}) does not match "
            f"number of database columns ({len(df.columns)})."
        )

    profile = pd.Series(STR_counts, index=df.columns)
    matches = df[(df == profile).all(axis=1)]
    return list(matches.index)

def group_3_STR_search(database_df: pd.DataFrame, dna_seq: str) -> list:
    STR_counts = [
        group_3_count_repeats_in_string(dna_seq, str(motif))
        for motif in database_df.columns
    ]
    return group_3_search_database(database_df, STR_counts)

############################## Group 4 ##############################

def group_4_count_repeats_in_string(string: str, substring: str) -> int:
    string = string.upper()
    substring = substring.upper()
    max_repeats = 0
    substring_len = len(substring)
    string_len = len(string)

    if substring_len == 0:
        return 0

    for i in range(string_len):
        current_repeats = 0
        current_pos = i
        
        while current_pos + substring_len <= string_len:
            if string[current_pos : current_pos + substring_len] \
            == substring:
                current_repeats += 1
                current_pos += substring_len
            else:
                break
        
        if current_repeats > max_repeats:
            max_repeats = current_repeats
            
    return max_repeats

def group_4_search_database(df: pd.DataFrame, STR_counts: list) -> list:
    matching_names = []
    
    for name, row in df.iterrows():
        row_list = list(row)
        
        is_match = True
        if len(row_list) != len(STR_counts):
            is_match = False
        else:
            for i in range(len(row_list)):
                if row_list[i] != STR_counts[i]:
                    is_match = False
                    break
        
        if is_match:
            matching_names.append(name)
            
    return matching_names

def group_4_STR_search(database_df: pd.DataFrame, dna_seq: str) -> list:
    dna_seq = dna_seq.upper()
    
    valid_bases = "ACGT"
    for char in dna_seq:
        is_valid = False
        for base in valid_bases:
            if char == base:
                is_valid = True
                break
        if not is_valid:
            return []

    str_markers = list(database_df.columns)
    
    query_profile = []
    for marker in str_markers:
        marker_upper = marker.upper()
        
        for char in marker_upper:
            is_valid = False
            for base in valid_bases:
                if char == base:
                    is_valid = True
                    break
            if not is_valid:
                raise []
        
        count = group_4_count_repeats_in_string(dna_seq, marker_upper)
        query_profile.append(count)
    
    results = group_4_search_database(database_df, query_profile)
    
    return results

############################## Group 5 ##############################

def group_5_count_repeats_in_string(string: str, substring: str) -> int:
    """
    Computes the longest run of consecutive repeats of a given substring 
    within a target string without using the 'in' operator or regex.
    """
    max_repeats = 0
    sub_len = len(substring)
    str_len = len(string)
    
    # Handle the edge case of an empty substring to avoid infinite loops
    if sub_len == 0:
        return 0

    # We use a while loop to avoid the 'in' keyword found in 'for i in range'
    i = 0
    while i <= (str_len - sub_len):
        current_run = 0
        
        # Check how many times the substring repeats consecutively starting at index i
        # We slice the string and compare it directly to our target substring
        while string[i + (current_run * sub_len) : i + ((current_run + 1) * sub_len)] == substring:
            current_run += 1
            
            # Stop if the next potential repeat would go past the end of the string
            if i + ((current_run + 1) * sub_len) > str_len:
                break
        
        # If this run is longer than our previous record, update max_repeats
        if current_run > max_repeats:
            max_repeats = current_run
        
        # Move to the next character in the string to check for a new run
        i += 1
        
    return max_repeats

def group_5_search_database(df: pd.DataFrame, STR_counts: list) -> list:
    """
    Searches a pandas DataFrame for rows that exactly match a list of STR counts.
    Returns a list of names (from the DataFrame index) that match.
    """
    
    # GUARD CLAUSE: If the number of STRs in the query doesn't match 
    # the number of columns in our database, there can be no match.
    if len(STR_counts) != len(df.columns):
        return []
        
    # 1. Compare every row in the DataFrame to the STR_counts list.
    # This creates a DataFrame of True/False values for every single cell.
    matches_mask = (df == STR_counts)
    
    # 2. Use .all(axis=1) to find rows where EVERY column matched (True across the whole row).
    # This results in a single Boolean Series.
    perfect_matches = matches_mask.all(axis=1)
    
    # 3. Use that Boolean Series to filter the index and convert the result to a list.
    match_names = df.index[perfect_matches].tolist()
    
    return match_names

def group_5_is_valid_dna(sequence: str) -> bool:
    """Helper to check if a string contains only valid DNA bases."""
    valid_bases = {'A', 'T', 'C', 'G'}
    # Standardize to upper for the check
    return all(base.upper() in valid_bases for base in sequence)
    
def group_5_STR_search(database_df: pd.DataFrame, dna_seq: str) -> list:
    """
    Searches a database for STR matches, with specific formatting 
    to satisfy the lab's test runner.
    """
    # 1. VALIDATION: The runner wants a RETURNED string, not a raised Error
    if not group_5_is_valid_dna(dna_seq):
        return "Non-valid DNA."
    
    # Standardize DNA to upper
    dna_seq = dna_seq.upper()
    
    str_targets = []
    for col in database_df.columns:
        if not group_5_is_valid_dna(col):
            return "Non-valid DNA."
        str_targets.append(col.upper())
        
    # 2. COUNTING
    found_counts = []
    for str_item in str_targets:
        count = group_5_count_repeats_in_string(dna_seq, str_item)
        found_counts.append(count)
        
    # 3. SEARCHING
    results = group_5_search_database(database_df, found_counts)
    
    # 4. FORMATTING FOR THE RUNNER:
    # If there is exactly one match, return it as a string instead of a list.
    if len(results) == 1:
        return results[0]
    
    return results

############################## Group 6 ##############################

def group_6_count_repeats_in_string(string: str, substring: str) -> int:
    def clean_manually(raw_text: str) -> str:
        cleaned = ""
        index = 0
        while index < len(raw_text):
            char = raw_text[index]
            if char.isalnum():
                cleaned += char.upper()
            index += 1
        return cleaned

    clean_string = clean_manually(string)
    clean_substring = clean_manually(substring)

    max_repeats = 0
    substring_len = len(clean_substring)
    string_len = len(clean_string)

    if substring_len == 0:
        return 0

    i = 0
    while i < string_len:
        current_run = 0
        j = i
        
        while clean_string[j : j + substring_len] == clean_substring:
            current_run += 1
            j += substring_len

            if j > string_len:
                break
        
        if current_run > max_repeats:
            max_repeats = current_run
        
        i += 1 

    return max_repeats

def group_6_search_database(database, query): # query is a list of STR count we are checking
    # 1. Strict Type Check: Ensure the query is exactly a list
    if not isinstance(query, list):
        return []
        
    # 2. Empty Checks: Return empty if the database is empty or the query is empty
    if database.empty or len(query) == 0:
        return []
        
    # 3. Shape Alignment Check: Ensure query length exactly matches the number of columns
    if len(query) != len(database.columns):
        return ["Insufficient STR counts in the list"]

    match_mask = (database == query).all(axis=1)

    matching_names = database.index[match_mask].tolist()
    return matching_names

def group_6_STR_search(database_df: pd.DataFrame, dna_seq: str) -> list:
    """
    Searches a database of STR profiles for any matches to a given DNA sequence.
    """
    valid_nucleotides = {'A', 'C', 'G', 'T'}
    
    # 1. Coerce the DNA sequence to uppercase
    dna_seq = dna_seq.upper()
    
    # 2. Validate the DNA sequence
    if not set(dna_seq).issubset(valid_nucleotides):
        raise ValueError("Invalid DNA sequence. Must only contain A, C, G, or T.")
        
    str_counts = []
    
    # 3. Iterate through each STR (the column headers in the database)
    for str_seq in database_df.columns:
        # Coerce STR to uppercase
        str_seq_upper = str(str_seq).upper()
        
        # Validate the STR sequence
        if not set(str_seq_upper).issubset(valid_nucleotides):
            raise ValueError(f"Invalid STR sequence in database column: {str_seq}")
            
        # 4. Count the repeats using the previously defined function
        count = group_6_count_repeats_in_string(dna_seq, str_seq_upper)
        str_counts.append(count)
        
    # 5. Search the database for matches using the previously defined function
    return group_6_search_database(database_df, str_counts)

############################## Group 7 ##############################

def group_7_count_repeats_in_string(string: str, substring: str) -> int:
    if not substring:
        return 0
        
    max_run = 0
    sub_len = len(substring)
    str_len = len(string)
    i = 0
    
    while i <= str_len - sub_len:
        # Check if we have a match at the current position
        if string[i : i + sub_len] == substring:
            current_run = 0
            
            # Count consecutive matches and JUMP the index
            while string[i : i + sub_len] == substring:
                current_run += 1
                i += sub_len
            
            if current_run > max_run:
                max_run = current_run
        else:
            # No match? Just move to the next character
            i += 1
            
    return max_run

def group_7_search_database(df: pd.DataFrame, STR_counts: list) -> list:
    matches = []
# Iterate through each person in the database
    for name, row in df.iterrows():
# Convert row values to a list to compare with query list
        if list(row) == STR_counts:
            matches.append(name)
    return matches

def group_7_STR_search(database_df: pd.DataFrame, dna_seq: str) -> list:
    # 1. Coerce DNA sequence to a consistent case (upper)
    dna_seq = dna_seq.upper()
    # 2. Calculate the longest run for each STR in the database columns
    current_str_counts = []
    for str_sequence in database_df.columns:
    # Coerce the STR sequence to the same case
        count = group_7_count_repeats_in_string(dna_seq, str_sequence.upper())
        current_str_counts.append(count)
    # 3. Search the database for these specific counts
    results = group_7_search_database(database_df, current_str_counts)
    return results

############################## Group 9 ##############################

def group_9_count_repeats_in_string(string: str, substring: str) -> int:
    """
    Returns the maximum number of times 'substring' repeats consecutively 
    within 'string'.
    """
    max_repeats = 0
    sub_len = len(substring)
    str_len = len(string)

    # Iterate through the string to find potential starting points
    for i in range(str_len):
        count = 0
        
        # While the substring matches at the current position, keep counting
        # and jumping forward by the length of the substring
        while True:
            start = i + count * sub_len
            end = start + sub_len
            
            if string[start:end] == substring:
                count += 1
            else:
                break
        
        # Update the maximum count found so far
        max_repeats = max(max_repeats, count)

    return max_repeats

def group_9_search_database(df: pd.DataFrame, STR_counts: list) -> list:
    """
    Compares the calculated STR counts against each person in the database.
    Returns a list containing the names of all individuals who match 
    the counts exactly.
    """
    matches = []
    
    # Iterate through the rows of the database
    # 'name' is the index, 'row' contains the STR counts for that person
    for name, row in df.iterrows():
        # Convert the row data to a list of integers for comparison
        # We ensure they are ints to match the output of count_repeats_in_string
        db_counts = list(row.astype(int))
        
        # Check if the calculated counts match the database row exactly
        if db_counts == STR_counts:
            matches.append(name)
            
    return matches

def group_9_STR_search(database_df: pd.DataFrame, dna_seq: str) -> list:
    """
    Orchestrates the DNA profiling process:
    1. Identifies which STRs to look for from the database columns.
    2. Calculates the longest consecutive repeat for each STR in the DNA sequence.
    3. Searches the database for individuals matching those counts.
    """
    # Get the list of STRs from the column names of the DataFrame
    # (Assuming the name is the index, the columns are the STR sequences)
    strs_to_search = list(database_df.columns)
    
    # Calculate the longest run for each STR found in the columns
    calculated_counts = []
    for seq in strs_to_search:
        count = group_9_count_repeats_in_string(dna_seq, seq)
        calculated_counts.append(count)
    
    # Use the previously defined search_database to find matches
    matches = group_9_search_database(database_df, calculated_counts)
    
    return matches


############################## Final Dicts ##############################

count_repeats_in_string = {
    "group_1": group_1_count_repeats_in_string,
    "group_2": group_2_count_repeats_in_string,
    "group_3": group_3_count_repeats_in_string,
    "group_4": group_4_count_repeats_in_string,
    "group_5": group_5_count_repeats_in_string,
    "group_6": group_6_count_repeats_in_string,
    "group_7": group_7_count_repeats_in_string,
    "group_9": group_9_count_repeats_in_string,
}

search_database = {
    "group_1": group_1_search_database,
    "group_2": group_2_search_database,
    "group_3": group_3_search_database,
    "group_4": group_4_search_database,
    "group_5": group_5_search_database,
    "group_6": group_6_search_database,
    "group_7": group_7_search_database,
    "group_9": group_9_search_database,
}

STR_search = {
    "group_1": group_1_STR_search,
    "group_2": group_2_STR_search,
    "group_3": group_3_STR_search,
    "group_4": group_4_STR_search,
    "group_5": group_5_STR_search,
    "group_6": group_6_STR_search,
    "group_7": group_7_STR_search,
    "group_9": group_9_STR_search,
}
