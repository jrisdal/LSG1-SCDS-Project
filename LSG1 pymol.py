##This script takes the input file of the Tolerance scores and maps it onto a .pdb file in PyMOL, and shows which residues are more or less
##tolerant to mutation.

#if you have trouble running this script in PyMOL, try just pasting each command into the command line, and then when you get to the for loop,
#run just that piece of code in a script form, and then change the spectrum by pasting it in as well.


#open up LSG1 structure in PyMOL before running. You can load it with commands but it's just easier with the GUI


import pandas as pd

scores_file = "C:/Users/joris/OneDrive/Desktop/SCDS project/Mutation Counting (code & output)/Actual output/Tolerance Scores.csv"
#path to scores file



data = pd.read_csv(scores_file) #reads scores file as a pandas df


cmd.alter("LSG1", "b=0.0")  # sets b factors to zero for all residues of LSG1

for index, row in data.iterrows(): #iterates over each row 
    residue_number = int(row['Residue']) #residue number
    cmd.alter(f"LSG1 and resi {residue_number}", f"b={row['Score']}") #turns the b factor into the score

#cmd.iterate(f"{'LSG1'}", "print(f'Residue {resi}: B-factor {b}')")
    #this checks the b factors to make sure they were actually altered


spectrum b, white_red #makes the coloring of the structure based off of the b factors, and on the scale from white to red
#white is closer to 1, red is closer to the max (1.034). Red is less tolerant to mutation

