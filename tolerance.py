#Takes two input csv files, produced with output_cleaner.py, and creates 'tolerance scores', which measure how resistant any residue is to
#mutation, by calculating the frequency of the WT AA at each codon in the init vs the prod, and then calculating prod/init. Values higher than 1
#are less resistant to change, values equal to or less than 1 are resistant to change and don't tolerate mutation



import os
import csv

#input files are counts

#tolerance("C:\\Users\\joris\\OneDrive\\Desktop\\SCDS project\\Mutation Counting (code & output)\\Actual output\\initAllFiltered_AA_counts_onlypossiblethruPCR.csv","C:\\Users\\joris\\OneDrive\\Desktop\\SCDS project\\Mutation Counting (code & output)\\Actual output\\prodAllFiltered_AA_counts_onlypossiblethruPCR.csv")



def tolerance(init_csv, prod_csv):

    output_directory = "C:/Users/joris/OneDrive/Desktop/SCDS project/Mutation Counting (code & output)" #change at your discretion
    output_filename = "Tolerance Scores.csv"
    output_file = os.path.join(output_directory, output_filename) #creates output file

    with open(init_csv, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        header_row = next(csv_reader) #skips header row


        init_tolerance = [] #initializes list that contains init freq's of WT allele
        

        for row in csv_reader: #iterates over each row in init_csv
            
            numeric_row = [int(value) for value in row[1:]] #converts counts to ints
            
            WT_AA = max(numeric_row) #finds the max value in the row

            row_sum = sum(int(value) for value in row[1:]) #sums all the values in the row

            init_tolerance.append(WT_AA/row_sum) #adds the frequency of the WT allele at that codon to a list for all the init codons

            

    
    with open(prod_csv, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        header_row = next(csv_reader) #skips header row

        prod_tolerance = [] #initializes list that contains prod freq's of WT allele

        for row in csv_reader: #iterates over each row in prod_csv
            
            numeric_row = [int(value) for value in row[1:]] #converts counts to ints
            
            WT_AA = max(numeric_row) #finds the max value in the row

            row_sum = sum(int(value) for value in row[1:]) #sums all the values in the row

            prod_tolerance.append(WT_AA/row_sum) #adds the freuqency of the WT allele at that codon to a list for all the prod codons




    tolerance_scores = [prod/init for prod,init in zip(prod_tolerance, init_tolerance)] #calculates the tolerance scores for each residue, by doing prod/init for each value in each list


    with open(output_file, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)

            for score in tolerance_scores: #writes scores to csv file
                csv_writer.writerow([score])








                
