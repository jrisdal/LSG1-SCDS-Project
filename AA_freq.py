##turns the AA counts in each row into frequencies/their proportion compared to the total counts in each row.
##use this after using output_cleaner.py


import os
import csv

def AA_freq(input_csv): #replace input_csv with path to input file. (output file from last function)

#AA_freq("C:\\Users\\joris\\OneDrive\\Desktop\\SCDS project\\Mutation Counting (code & output)\\Actual output\\initAllFiltered_AA_counts_onlypossiblethruPCR.csv")
    
    with open(input_csv, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file) #reads csv file
        #header_row = next(csv_reader) #skips first row


        output_directory = "C:/Users/joris/OneDrive/Desktop/SCDS project/Mutation Counting (code & output)/Actual output" #change at your discretion
        output_filename = input_csv.replace("AllFiltered_AA_counts_onlypossiblethruPCR.csv", "AAfreqs.csv") 
        output_file = os.path.join(output_directory, output_filename) #creates output file
    
        freq_csv = [] #initializes list to store final freqs

        
        for row in csv_reader: #iterates over each row in the csv file
            

            
            current_codon = {} #initializes/resets the dictionary for each codon/row
            
            row_sum = sum(int(value) for value in list(row.values())[1:]) #sums all the values in each row as integers
            row_max = max(int(value) for value in list(row.values())[1:]) #calculates maximum of each row, which is WT AA 

            for aa, count in list(row.items())[1:]: #iterates over each row, but skips the first column
                if int(count) == row_max:
                    aa_freq = "WT"
                    current_codon[aa] = aa_freq
                    continue
                else:
                    aa_freq = int(count) / row_sum #calculates AA freq as its count over the total row sum
                    current_codon[aa] = aa_freq #adds that frequency to the dictionary at the correct AA column
   

                
            freq_csv.append(current_codon) #appends the current codon dicitonary to the total freq list



                
    with open(output_file, 'w', newline='') as csvfile:
            fieldnames = ['Codon', 'A', 'R', 'D', 'N', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'nonsense'] #column names
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            for row in freq_csv:
                writer.writerow(row)

    print(f"Output file: {output_file}, Output directory: {output_directory}")
