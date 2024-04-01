##Creates fitness scores by doing prod freq/init freq
##Run after AA_freq.py


import csv
import os


#fitness("C:\\Users\\joris\\OneDrive\\Desktop\\SCDS project\\Mutation Counting (code & output)\\Actual output\\initAAfreqs.csv","C:\\Users\\joris\\OneDrive\\Desktop\\SCDS project\\Mutation Counting (code & output)\\Actual output\\prodAAfreqs.csv")

def fitness(init_csv, prod_csv):
    output_directory = "C:/Users/joris/OneDrive/Desktop/SCDS project/Mutation Counting (code & output)/Actual output"
    output_filename = "Fitness ScoresTEST.csv"
    output_file = os.path.join(output_directory, output_filename)

    with open(init_csv, 'r') as init_file, open(prod_csv, 'r') as prod_file, open(output_file, 'w', newline='') as output:
        init_reader = csv.reader(init_file)
        prod_reader = csv.reader(prod_file)
        output_writer = csv.writer(output)

        next(init_reader) #skips the first row in both the input and csv files, which contains what AA the row is
        next(prod_reader)

        for init_row, prod_row in zip(init_reader, prod_reader):
            

            # Assuming codon number is in the first column (index 0)
            
            init_values = []
            for i, value in enumerate(init_row):
                if value == "WT": #checks for if a value is a WT value
                    init_values.append(value) 
                elif i != 0:
                    init_values.append(float(value)) #checks to make sure a value isn't zero and app
                else:
                    init_values.append(0.0)

            prod_values = []
            for i, value in enumerate(prod_row):
                if value == "WT":
                    prod_values.append(value)
                elif i != 0: 
                    prod_values.append(float(value)) #checks to make sure a value isn't zero to avoid division by zero
                else:
                    prod_values.append(0.0)


            
            #init_values = [float(value) if i != 0 else 0.0 for i, value in enumerate(init_row)]
            #prod_values = [float(value) if i != 0 else 0.0 for i, value in enumerate(prod_row)]

            # Avoid division by zero
            fitness_scores = []
            for prod, init in zip(prod_values, init_values):
                if init == "WT":
                    fitness_scores.append(init)
                elif init != 0:
                    fitness_scores.append(prod / init)
                else:
                    fitness_scores.append(0.0)


            
            #fitness_scores = [prod / init if init != 0 else 0.0 for prod, init in zip(prod_values, init_values)] #avoids dividing by zero, which is common

            # Write the fitness scores to the output CSV file
            output_writer.writerow(fitness_scores) #writes fitness scores to output csv file




