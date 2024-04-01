##Takes an input .csv file created from the mut_counter.py function and creates a filtered output .csv file that only contains
##the AA substitution counts that should theoretically be possible through PCR mutagenesis, ie those up and down on the codon chart
##from the original WT AA, and those that are at a high enough frequency to be considered valid (threshold is 0.0005).



import csv
import os
import math

#output_cleaner("C:\\Users\\joris\\OneDrive\\Desktop\\SCDS project\\Mutation Counting (code & output)\\Actual output\\initAllFiltered_AA_counts.csv")
def output_cleaner(input_csv):
    #replace input_csv with path to csv file

    #input_filename = os.path.basename(input_csv)
    output_directory = "C:/Users/joris/OneDrive/Desktop/SCDS project/Mutation Counting (code & output)" #change at your discretion
    output_filename = input_csv.replace(".csv", "_onlypossiblethruPCR.csv") 
    output_file = os.path.join(output_directory, output_filename) #creates output file

    altered_csv = [] #initializes altered csv output list

    A = 1
    R = 2 #initializes each of the amino acids as numbers, so that it can be applied to the input csv file which has columns of each AA
    D = 3
    N = 4
    C = 5
    E = 6
    Q = 7
    G = 8
    H = 9
    I = 10
    L = 11
    K = 12
    M = 13
    F = 14
    P = 15
    S = 16
    T = 17
    W = 18
    Y = 19
    V = 20
    nonsense = 21



    with open(input_csv, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        header_row = next(csv_reader) #skips header row
        
        for row in csv_reader:
            # Convert numerical values in the row to integers for comparison
            numeric_row = [int(value) for value in row[1:]]
            
            # Finds the index (header/name) of the maximum value in the numeric row
            WT_AA = numeric_row.index(max(numeric_row))
            
            # Get the corresponding column name from the header row
            AA = header_row[WT_AA + 1]


            print(AA)
            current_aa_counts = {    #resets dictionary of each codon AA count
                        'A': 0,
                        'R': 0,
                        'D': 0,
                        'N': 0,
                        'C': 0,
                        'E': 0,
                        'Q': 0,
                        'G': 0,
                        'H': 0,
                        'I': 0,
                        'L': 0,
                        'K': 0,
                        'M': 0,
                        'F': 0,
                        'P': 0,
                        'S': 0,
                        'T': 0,
                        'W': 0,
                        'Y': 0,
                        'V': 0,
                        'nonsense': 0
            }               
            if AA == "A": #appends the counts of all the amino acids, but only if they are counts that theoretically possible thru PCR mutagenesis
                current_aa_counts['A'] = row[A]
                
                current_aa_counts['V'] = row[V]
                current_aa_counts['D'] = row[D]
                current_aa_counts['E'] = row[E]
                current_aa_counts['G'] = row[G] #those which are up and down on the codon chart from alanine, and all other AA's 
                current_aa_counts['S'] = row[S]
                current_aa_counts['P'] = row[P]
                current_aa_counts['T'] = row[T]
                
            if AA == "R":
                current_aa_counts['R'] = row[R]
                
                current_aa_counts['N'] = row[N]
                current_aa_counts['C'] = row[C]
                current_aa_counts['Q'] = row[Q]
                current_aa_counts['G'] = row[G]
                current_aa_counts['H'] = row[H]
                current_aa_counts['I'] = row[I]
                current_aa_counts['L'] = row[L]
                current_aa_counts['K'] = row[K]
                current_aa_counts['M'] = row[M]
                current_aa_counts['P'] = row[P]
                current_aa_counts['S'] = row[S]
                current_aa_counts['T'] = row[T]
                current_aa_counts['W'] = row[W]
                current_aa_counts['nonsense'] = row[nonsense]
                
                
            if AA == "D":
                current_aa_counts['D'] = row[D]
                
                current_aa_counts['A'] = row[A]
                current_aa_counts['N'] = row[N]
                current_aa_counts['E'] = row[E]
                current_aa_counts['Q'] = row[Q]
                current_aa_counts['G'] = row[G]
                current_aa_counts['H'] = row[H]
                current_aa_counts['K'] = row[K]
                current_aa_counts['Y'] = row[Y]
                current_aa_counts['V'] = row[V]
                current_aa_counts['nonsense'] = row[nonsense]
                
            if AA == "N":
                current_aa_counts['N'] = row[N]
                
                current_aa_counts['R'] = row[R]
                current_aa_counts['D'] = row[D]
                current_aa_counts['E'] = row[E]
                current_aa_counts['Q'] = row[Q]
                current_aa_counts['H'] = row[H]
                current_aa_counts['I'] = row[I]
                current_aa_counts['K'] = row[K]
                current_aa_counts['M'] = row[M]
                current_aa_counts['S'] = row[S]
                current_aa_counts['T'] = row[T]
                current_aa_counts['Y'] = row[Y]
                current_aa_counts['nonsense'] = row[nonsense]
                
            if AA == "C":
                current_aa_counts['C'] = row[C]

                current_aa_counts['R'] = row[R]
                current_aa_counts['G'] = row[G]
                current_aa_counts['L'] = row[L]
                current_aa_counts['F'] = row[F]
                current_aa_counts['S'] = row[S]
                current_aa_counts['W'] = row[W]
                current_aa_counts['Y'] = row[Y]
                current_aa_counts['nonsense'] = row[nonsense]
                
            if AA == "E":
                current_aa_counts['E'] = row[E]

                current_aa_counts['A'] = row[A]
                current_aa_counts['D'] = row[D]
                current_aa_counts['N'] = row[N]
                current_aa_counts['Q'] = row[Q]
                current_aa_counts['G'] = row[G]
                current_aa_counts['H'] = row[H]
                current_aa_counts['K'] = row[K]
                current_aa_counts['Y'] = row[Y]
                current_aa_counts['V'] = row[V]
                current_aa_counts['nonsense'] = row[nonsense]
                
            if AA == "Q":
                current_aa_counts['Q'] = row[Q]
                
                current_aa_counts['R'] = row[R]
                current_aa_counts['D'] = row[D]
                current_aa_counts['N'] = row[N]
                current_aa_counts['E'] = row[E]
                current_aa_counts['H'] = row[H]
                current_aa_counts['L'] = row[L]
                current_aa_counts['K'] = row[K]
                current_aa_counts['P'] = row[P]
                current_aa_counts['Y'] = row[Y]
                current_aa_counts['nonsense'] = row[nonsense]
                
            if AA == "G":
                current_aa_counts['G'] = row[G]
                
                current_aa_counts['A'] = row[A]
                current_aa_counts['R'] = row[R]
                current_aa_counts['D'] = row[D]
                current_aa_counts['C'] = row[C]
                current_aa_counts['E'] = row[E]
                current_aa_counts['S'] = row[S]
                current_aa_counts['W'] = row[W]
                current_aa_counts['V'] = row[V]
                current_aa_counts['nonsense'] = row[nonsense]
                
            if AA == "H":
                current_aa_counts['H'] = row[H]
                
                current_aa_counts['R'] = row[R]
                current_aa_counts['D'] = row[D]
                current_aa_counts['N'] = row[N]
                current_aa_counts['E'] = row[E]
                current_aa_counts['Q'] = row[Q]
                current_aa_counts['L'] = row[L]
                current_aa_counts['K'] = row[K]
                current_aa_counts['P'] = row[P]
                current_aa_counts['Y'] = row[Y]
                current_aa_counts['nonsense'] = row[nonsense]
                
            if AA == "I":
                current_aa_counts['I'] = row[I]
                
                current_aa_counts['R'] = row[R]
                current_aa_counts['N'] = row[N]
                current_aa_counts['L'] = row[L]
                current_aa_counts['K'] = row[K]
                current_aa_counts['M'] = row[M]
                current_aa_counts['F'] = row[F]
                current_aa_counts['S'] = row[S]
                current_aa_counts['T'] = row[T]
                current_aa_counts['V'] = row[V]
                
            if AA == "L":
                current_aa_counts['L'] = row[L]
                
                current_aa_counts['R'] = row[R]
                current_aa_counts['C'] = row[C]
                current_aa_counts['Q'] = row[Q]
                current_aa_counts['H'] = row[H]
                current_aa_counts['I'] = row[I]
                current_aa_counts['M'] = row[M]
                current_aa_counts['F'] = row[F]
                current_aa_counts['P'] = row[P]
                current_aa_counts['S'] = row[S]
                current_aa_counts['W'] = row[W]
                current_aa_counts['Y'] = row[Y]
                current_aa_counts['V'] = row[V]
                current_aa_counts['nonsense'] = row[nonsense]
                
            if AA == "K":
                current_aa_counts['K'] = row[K]
                
                current_aa_counts['R'] = row[R]
                current_aa_counts['D'] = row[D]
                current_aa_counts['N'] = row[N]
                current_aa_counts['E'] = row[E]
                current_aa_counts['Q'] = row[Q]
                current_aa_counts['H'] = row[H]
                current_aa_counts['I'] = row[I]
                current_aa_counts['M'] = row[M]
                current_aa_counts['S'] = row[S]
                current_aa_counts['T'] = row[T]
                current_aa_counts['Y'] = row[Y]
                current_aa_counts['nonsense'] = row[nonsense]
                
            if AA == "M":
                current_aa_counts['M'] = row[M]
                
                current_aa_counts['R'] = row[R]
                current_aa_counts['N'] = row[N]
                current_aa_counts['I'] = row[I]
                current_aa_counts['L'] = row[L]
                current_aa_counts['K'] = row[K]
                current_aa_counts['F'] = row[F]
                current_aa_counts['S'] = row[S]
                current_aa_counts['T'] = row[T]
                current_aa_counts['V'] = row[V]
                
            if AA == "F":
                current_aa_counts['F'] = row[F]
                
                current_aa_counts['C'] = row[C]
                current_aa_counts['I'] = row[I]
                current_aa_counts['L'] = row[L]
                current_aa_counts['M'] = row[M]
                current_aa_counts['S'] = row[S]
                current_aa_counts['W'] = row[W]
                current_aa_counts['Y'] = row[Y]
                current_aa_counts['V'] = row[V]
                current_aa_counts['nonsense'] = row[nonsense]
                
                
            if AA == "P":
                current_aa_counts['P'] = row[P]
                
                current_aa_counts['A'] = row[A]
                current_aa_counts['R'] = row[R]
                current_aa_counts['Q'] = row[Q]
                current_aa_counts['H'] = row[H]
                current_aa_counts['L'] = row[L]
                current_aa_counts['S'] = row[S]
                current_aa_counts['T'] = row[T]
                
            if AA == "S":
                current_aa_counts['S'] = row[S]
                
                current_aa_counts['A'] = row[A]
                current_aa_counts['R'] = row[R]
                current_aa_counts['N'] = row[N]
                current_aa_counts['C'] = row[C]
                current_aa_counts['G'] = row[G]
                current_aa_counts['I'] = row[I]
                current_aa_counts['L'] = row[L]
                current_aa_counts['K'] = row[K]
                current_aa_counts['M'] = row[M]
                current_aa_counts['F'] = row[F]
                current_aa_counts['P'] = row[P]
                current_aa_counts['T'] = row[T]
                current_aa_counts['W'] = row[W]
                current_aa_counts['Y'] = row[Y]
                current_aa_counts['nonsense'] = row[nonsense]
                
            if AA == "T":
                current_aa_counts['T'] = row[T]
                
                current_aa_counts['A'] = row[A]
                current_aa_counts['R'] = row[R]
                current_aa_counts['N'] = row[N]
                current_aa_counts['I'] = row[I]
                current_aa_counts['K'] = row[K]
                current_aa_counts['M'] = row[M]
                current_aa_counts['P'] = row[P]
                current_aa_counts['S'] = row[S]
                
            if AA == "W":
                current_aa_counts['W'] = row[W]
                
                current_aa_counts['R'] = row[R]
                current_aa_counts['C'] = row[C]
                current_aa_counts['G'] = row[G]
                current_aa_counts['L'] = row[L]
                current_aa_counts['F'] = row[F]
                current_aa_counts['S'] = row[S]
                current_aa_counts['Y'] = row[Y]
                current_aa_counts['nonsense'] = row[nonsense]
                
            if AA == "Y":
                current_aa_counts['Y'] = row[Y]
                
                current_aa_counts['D'] = row[D]
                current_aa_counts['N'] = row[N]
                current_aa_counts['C'] = row[C]
                current_aa_counts['E'] = row[E]
                current_aa_counts['Q'] = row[Q]
                current_aa_counts['H'] = row[H]
                current_aa_counts['L'] = row[L]
                current_aa_counts['K'] = row[K]
                current_aa_counts['F'] = row[F]
                current_aa_counts['S'] = row[S]
                current_aa_counts['W'] = row[W]
                current_aa_counts['nonsense'] = row[nonsense]
                
                
            if AA == "V":
                current_aa_counts['V'] = row[V]
                
                current_aa_counts['A'] = row[A]
                current_aa_counts['D'] = row[D]
                current_aa_counts['E'] = row[E]
                current_aa_counts['G'] = row[G]
                current_aa_counts['I'] = row[I]
                current_aa_counts['L'] = row[L]
                current_aa_counts['M'] = row[M]
                current_aa_counts['F'] = row[F]

            current_aa_counts = { #turns all the remaining aa counts in the current_aa_counts dictionary into integers, so that low frequency counts can be removed
                'A': int(row[A]),
                'R': int(row[R]),
                'D': int(row[D]),
                'N': int(row[N]),
                'C': int(row[C]),
                'E': int(row[E]),
                'Q': int(row[Q]),
                'G': int(row[G]),
                'H': int(row[H]),
                'I': int(row[I]),
                'L': int(row[L]),
                'K': int(row[K]),
                'M': int(row[M]),
                'F': int(row[F]),
                'P': int(row[P]),
                'S': int(row[S]),
                'T': int(row[T]),
                'W': int(row[W]),
                'Y': int(row[Y]),
                'V': int(row[V]),
                'nonsense': int(row[nonsense])
            }

                
            total_sum = sum(current_aa_counts.values()) #sums up the total number of value in the current AA counts dictionary
            
            for key, value in list(current_aa_counts.items()): #iterates over the dictionary of only pcr mutagenesis possible counts, and removes those at a frequency lower than 0.0001
                if value / total_sum < 0.0005: #if a value occurs at a frequency less than 0.0005...
                    current_aa_counts[key] = 0 #its turned into zero
            
            altered_csv.append(current_aa_counts) #appends only valid counts to altered_csv


    with open(output_file, 'w', newline='') as csvfile:
            fieldnames = ['Codon', 'A', 'R', 'D', 'N', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'nonsense'] #column names
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            for row in altered_csv:
                writer.writerow(row)

    print(f"Output file: {output_file} Output directory: {output_directory}")




