##Takes an input .sam alignment file and counts the mutations in the reads one codon at a time, line by line.##
#notes: 
#       -COI refers to codon of interest in comments
#       -read length in this experiment is 75nt, LSG1 is 640AA
#       -on this computer, analyzing initL1, it takes ~10sec for one codon. So for all 640AA's, it should take ~2hrs


import os
import csv
import time #just for tracking how long this takes


def mut_counter(input_sam): #change input_sam to path of input sam file
                            #future arguments: read length, sequence length in AA, output directory

    start_time = time.time()
    
    
    output_directory = "C:/Users/joris/OneDrive/Desktop/SCDS project/Mutation Counting" #change at your discretion

    input_filename = os.path.basename(input_sam)
    output_filename = input_filename.replace(".sam", "_AA_counts.csv") 

    output_file = os.path.join(output_directory, output_filename) #creates output file

    #reference LSG1 sequence and AA sequence. Just for my own reference, ended up not being necessary as counting WT codons is valuable info in itself.
    #reference = "ATGCCACCAAAAGAAGCTCCCAAGAAATGGAAGGCGCCAAAAGGGCCAAAACCTACCCACCGTAAAAATAAAAATAAGCTTGAATTAGGCAGAGCTATTAAATATGCACGTCAAAAAGAAAATGCCATCGAGTATTTACCTGATGGTGAAATGAGGTTCACTACCGATAAGCATGAGGCCAACTGGGTTAAATTAAGATCTGTAACTCAAGAATCAGCTTTAGATGAATTCTTGAGTACAGCTGCACTGGCAGACAAAGATTTCACGGCCGATAGACATTCAAATGTTAAAATTATTAGAATGGATAGCGGTAATGATTCTGCGACATCTCAAGGGTTTTCTATGACTAATGAGCAGCGTGGAAATCTTAATGCGAAGCAAAGAGCGCTTGCTAAGGATTTGATTGTTCCAAGGAGGCCTGAATGGAACGAGGGCATGTCCAAGTTTCAGCTTGATAGGCAAGAAAAGGAAGCGTTTTTAGAATGGAGAAGAAAATTGGCACATTTACAAGAAAGCAATGAAGACTTGTTGTTAACACCGTTTGAAAGAAATATCGAAGTTTGGAAACAGTTATGGAGAGTTGTTGAAAGATCAGATTTAGTTGTTCAAATTGTAGATGCGAGGAATCCGTTGCTGTTTAGATCTGTCGATTTAGAAAGATATGTAAAAGAGTCAGATGACAGAAAAGCAAACTTACTGCTAGTTAATAAAGCAGATTTATTGACCAAAAAGCAACGTATCGCTTGGGCAAAGTACTTTATCTCCAAGAATATTTCGTTCACGTTTTACTCTGCATTGAGAGCTAATCAATTATTGGAGAAACAAAAGGAAATGGGGGAAGATTATAGAGAACAAGATTTCGAGGAAGCTGATAAAGAAGGGTTCGATGCTGATGAAAAAGTTATGGAAAAAGTTAAAATTCTGTCCATTGACCAACTGGAAGAATTGTTTTTATCAAAAGCTCCAAACGAGCCTTTATTGCCACCTCTGCCCGGTCAACCTCCACTGATTAATATTGGTTTGGTTGGTTATCCAAATGTAGGTAAATCCTCCACTATTAATTCGCTCGTGGGTGCCAAGAAAGTTTCTGTTTCATCCACGCCTGGTAAAACAAAACACTTCCAAACTATTAAGTTATCTGATTCTGTCATGCTTTGTGACTGTCCCGGTCTTGTCTTCCCAAACTTTGCATATAACAAGGGTGAGCTCGTGTGTAATGGTGTTTTACCTATTGATCAATTACGTGATTATATTGGTCCAGCAGGTTTAGTAGCAGAAAGAATACCAAAGTATTACATTGAAGCGATTTATGGTATCCATATTCAAACCAAATCAAGGGATGAAGGTGGGAATGGGGATATACCAACTGCTCAAGAATTGTTAGTCGCATACGCCAGAGCTCGTGGTTATATGACTCAAGGTTACGGTTCTGCTGATGAACCAAGAGCAAGTCGTTATATATTGAAAGATTACGTCAATGGGAAATTACTGTATGTCAACCCTCCGCCCCATCTAGAGGATGATACACCTTACACTAGAGAAGAGTGTGAAGAATTTAACAAAGATTTATATGTGTTCGACAGATTACCGGACACCAGAAAGGAGCAAGTGCAAAATGCTGCTAAGGCTAAAGGCATTGATATCGTGGATTTAGCTCGTGATTTGAATCAGCTAACGTTTTCAGCTCACACTGGTGGTGACACACAAAAAGAAGCCAAATCTGTTACGCACGGTGGTAAACAAGCTGCATTGTACAATGCCGCTGAGGACTTGGATAGAGACTTCTTCAAGATGAACAACGTCGAAGGTAGATTAAGTACACCATTCCACAAAGTTCAAAATAGTTCAGCTGGTAAGAGACATAACAAAAAAAACAAAAGTAAAAATGCGAAAAGCAAAGTTTTTAGCATTGAAAATAATTAG"
    #ref_AA = "MPPKEAPKKWKAPKGPKPTHRKNKNKLELGRAIKYARQKENAIEYLPDGEMRFTTDKHEANWVKLRSVTQESALDEFLSTAALADKDFTADRHSNVKIIRMDSGNDSATSQGFSMTNEQRGNLNAKQRALAKDLIVPRRPEWNEGMSKFQLDRQEKEAFLEWRRKLAHLQESNEDLLLTPFERNIEVWKQLWRVVERSDLVVQIVDARNPLLFRSVDLERYVKESDDRKANLLLVNKADLLTKKQRIAWAKYFISKNISFTFYSALRANQLLEKQKEMGEDYREQDFEEADKEGFDADEKVMEKVKILSIDQLEELFLSKAPNEPLLPPLPGQPPLINIGLVGYPNVGKSSTINSLVGAKKVSVSSTPGKTKHFQTIKLSDSVMLCDCPGLVFPNFAYNKGELVCNGVLPIDQLRDYIGPAGLVAERIPKYYIEAIYGIHIQTKSRDEGGNGDIPTAQELLVAYARARGYMTQGYGSADEPRASRYILKDYVNGKLLYVNPPPHLEDDTPYTREECEEFNKDLYVFDRLPDTRKEQVQNAAKAKGIDIVDLARDLNQLTFSAHTGGDTQKEAKSVTHGGKQAALYNAAEDLDRDFFKMNNVEGRLSTPFHKVQNSSAGKRHNKKNKSKNAKSKVFSIENN"

    
    
    codon_to_amino_acid = {
                'GCT': 'A',
                'GCC': 'A',
                'GCA': 'A',
                'GCG': 'A',
                
                'CGT': 'R',
                'CGC': 'R',
                'CGA': 'R',
                'CGG': 'R',
                'AGA': 'R',
                'AGG': 'R',
                
                'GAT': 'D',
                'GAC': 'D',
                
                'AAT': 'N',
                'AAC': 'N',
                
                'TGT': 'C',
                'TGC': 'C',
                
                'GAA': 'E',
                'GAG': 'E',
                
                'CAA': 'Q',
                'CAG': 'Q',
                
                'GGT': 'G',
                'GGC': 'G',
                'GGA': 'G',
                'GGG': 'G',
                
                'CAT': 'H',
                'CAC': 'H',
                
                'ATT': 'I',
                'ATC': 'I',
                'ATA': 'I',
                
                'TTA': 'L',
                'TTG': 'L',
                'CTT': 'L',
                'CTC': 'L',
                'CTA': 'L',
                'CTG': 'L',
                
                'AAA': 'K',
                'AAG': 'K',
                
                'ATG': 'M',
                
                'TTT': 'F',
                'TTC': 'F',
                
                'CCT': 'P',
                'CCC': 'P',
                'CCA': 'P',
                'CCG': 'P',
                
                'TCT': 'S',
                'TCC': 'S',
                'TCA': 'S',
                'TCG': 'S',
                'AGT': 'S',
                'AGC': 'S',
                
                'ACT': 'T',
                'ACC': 'T',
                'ACA': 'T',
                'ACG': 'T',
                
                'TGG': 'W',
                
                'TAT': 'Y',
                'TAC': 'Y',
                
                'GTT': 'V',
                'GTC': 'V',
                'GTA': 'V',
                'GTG': 'V',
                
                'TAA': 'nonsense',
                'TAG': 'nonsense',
                'TGA': 'nonsense',
}
      
    

    all_aa_counts = [] #initializes all the aa counts as a list

   
    with open(input_sam, 'r') as sam_file:
        
        codon_count = 2 #sets codon count, don't consider codon 1 since thats the start codon and it won't change
        codon_of_interest = "" #initializes COI
        
        
        while codon_count <= 640: #specific to LSG1. could in the future make another argument in main that takes into account how many codons there are
                                #change at your discretion
            
            sam_file.seek(0) #resets file pointer
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
            
            for line in sam_file:
                if line.startswith('@'): #if sam file has header line, skips it. not necessary in this case but often is
                    continue
                
                fields = line.strip().split('\t') #splits lines into fields separated by tabs

                #if int(fields[4]) < 40: #skips reads with a low mapping quality or no mapping. Change at your own discretion. 
                    #continue
                ##Already trimmed low quality mappings using samtools, unecessary now ^. 
                #Recommend that you do this before running this program to cut down on runtime

               
                start_pos = int(fields[3]) #initializes starting position of each read               
                
                
                if start_pos > codon_count*3 or start_pos + 75 < codon_count * 3: #read doesn't contain COI because its too early/COI is after read length
                   continue

                if start_pos == codon_count*3 + 1 or start_pos == codon_count*3 + 2 or start_pos == codon_count*3 - 1 or start_pos == codon_count*3 - 2:
                    continue #omits reads that contain only 1 or 2 nts of the codon of interest and so are out of frame


                #elif start_pos * 3 == codon_count:  
                 #  codon_of_interest = fields[9][:3]  # extracts first three nts from sequence line if the start pos is where the codon of interest is & is in frame
                   #print("it worked")
                
                for i in range(0,76): #could probably do this in a more computationally efficient manner but couldn't think of how
                        if (start_pos + i) == codon_count*3:  #checks to see where in the read the codon of interest is
                            codon_of_interest = fields[9][(i-2):(i+1)]#not honestly sure why this works but it does, don't touch
                            break
                            



                #not functional but also not necessary tbh       
                #if codon_of_interest == reference[codon_count*3: codon_count*3 + 3]: #skips codon if its identical to ref
                   #print("codon skipped")
                   #continue
                    

                
                if codon_of_interest in codon_to_amino_acid:
                   current_aa_counts[codon_to_amino_acid[codon_of_interest]] += 1 #writes what the codon of interest is into a dictionary of AA counts
                else:
                   continue

            current_aa_counts['Codon'] = codon_count  #adds codon information to the dictionary 
            all_aa_counts.append(current_aa_counts)
            #adds the current aa counts for that specific codon into the list with all aa counts
            codon_count += 1

            if codon_count <= 640:
                print("AA content of codon " + str(codon_count) + " of 640 counted") #just for viewing progress as program runs
            
    #at the end, combines all AA counts into one csv file

        
        with open(output_file, 'w', newline='') as csvfile:
            fieldnames = ['Codon', 'A', 'R', 'D', 'N', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'nonsense'] #column names
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            for row in all_aa_counts:
                writer.writerow(row)

        print(f"Output file: {output_file} Output directory: {output_directory}")
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Script execution time: {elapsed_time} seconds") #shows how much time elapsed during runtime
        

       
