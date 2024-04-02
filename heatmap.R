##creating a 640x20 tile with geom_tile() to represent the fitness scores of LSG1 AA substitutions

library(ggplot2)
library(colorspace)
library(reshape2)
library(dplyr)

#reference AA seq: MPPKEAPKKWKAPKGPKPTHRKNKNKLELGRAIKYARQKENAIEYLPDGEMRFTTDKHEANWVKLRSVTQESALDEFLSTAALADKDFTADRHSNVKIIRMDSGNDSATSQGFSMTNEQRGNLNAKQRALAKDLIVPRRPEWNEGMSKFQLDRQEKEAFLEWRRKLAHLQESNEDLLLTPFERNIEVWKQLWRVVERSDLVVQIVDARNPLLFRSVDLERYVKESDDRKANLLLVNKADLLTKKQRIAWAKYFISKNISFTFYSALRANQLLEKQKEMGEDYREQDFEEADKEGFDADEKVMEKVKILSIDQLEELFLSKAPNEPLLPPLPGQPPLINIGLVGYPNVGKSSTINSLVGAKKVSVSSTPGKTKHFQTIKLSDSVMLCDCPGLVFPNFAYNKGELVCNGVLPIDQLRDYIGPAGLVAERIPKYYIEAIYGIHIQTKSRDEGGNGDIPTAQELLVAYARARGYMTQGYGSADEPRASRYILKDYVNGKLLYVNPPPHLEDDTPYTREECEEFNKDLYVFDRLPDTRKEQVQNAAKAKGIDIVDLARDLNQLTFSAHTGGDTQKEAKSVTHGGKQAALYNAAEDLDRDFFKMNNVEGRLSTPFHKVQNSSAGKRHNKKNKSKNAKSKVFSIENN

setwd("C:/Users/joris/OneDrive/Desktop/SCDS project/Mutation Counting (code & output)/Actual output")

fitness_csv <- read.csv("Fitness Scores.csv") #reads fitness .csv file
fitness_csv[fitness_csv==0] <- NA #turns all the zeroes in the .csv file into NA
fitness_csv[fitness_csv=="WT"] <- 1 #turns all WT into 1 because those theoretically have neutral fitness
melted_data <- reshape2::melt(fitness_csv, id.vars = "Codon") #melts data for easier plotting
View(melted_data) 



WT_csv <- read.csv("Fitness Scores w WT.csv")
View(WT_csv)

#reorderd variable column in order to have y axis properly ordered
desired_order = c("A","R","D","N","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","nonsense")
WT_csv$variable <- factor(WT_csv$variable, levels = desired_order) 


wtseq <- "PPKEAPKKWKAPKGPKPTHRKNKNKLELGRAIKYARQKENAIEYLPDGEMRFTTDKHEANWVKLRSVTQESALDEFLSTAALADKDFTADRHSNVKIIRMDSGNDSATSQGFSMTNEQRGNLNAKQRALAKDLIVPRRPEWNEGMSKFQLDRQEKEAFLEWRRKLAHLQESNEDLLLTPFERNIEVWKQLWRVVERSDLVVQIVDARNPLLFRSVDLERYVKESDDRKANLLLVNKADLLTKKQRIAWAKYFISKNISFTFYSALRANQLLEKQKEMGEDYREQDFEEADKEGFDADEKVMEKVKILSIDQLEELFLSKAPNEPLLPPLPGQPPLINIGLVGYPNVGKSSTINSLVGAKKVSVSSTPGKTKHFQTIKLSDSVMLCDCPGLVFPNFAYNKGELVCNGVLPIDQLRDYIGPAGLVAERIPKYYIEAIYGIHIQTKSRDEGGNGDIPTAQELLVAYARARGYMTQGYGSADEPRASRYILKDYVNGKLLYVNPPPHLEDDTPYTREECEEFNKDLYVFDRLPDTRKEQVQNAAKAKGIDIVDLARDLNQLTFSAHTGGDTQKEAKSVTHGGKQAALYNAAEDLDRDFFKMNNVEGRLSTPFHKVQNSSAGKRHNKKNKSKNAKSKVFSIENN"
wtseq <- strsplit(wtseq, "")[[1]] 
WT_csv$wt_sequence <- wtseq #creates WT seq column for labeling of x axis

View(WT_csv)

##creates the first heatmap from codon 2 to 80, as it's only 79 codons long and so needs to be dealt with separately

ggplot() +
  geom_tile(data = WT_csv, aes(x = factor(Codon), y = factor(variable), fill= as.numeric(score))) +
  scale_fill_gradient2(low= "brown", mid = "lightyellow", high= "darkgreen", midpoint = 1, name = "Fitness Scores") +
  coord_fixed(ratio = 1, xlim = c(1,79)) +
  scale_x_discrete(breaks = seq(2, 80, by=1), labels = paste(seq(2, 80, by=1), wtseq[1:79])) +
  labs(x = "Position & WT AA", y = "Amino Acid", title = paste("Heatmap of Substitution Fitness Scores in LSG1 (Codons 2 to 80)")) +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#function to create all other heatmaps
generate_heatmap <- function(start_codon) {
  ggplot() +
    geom_tile(data = WT_csv, aes(x = factor(Codon), y = factor(variable), fill= as.numeric(score))) +
    scale_fill_gradient2(low= "brown", mid = "lightyellow", high= "darkgreen", midpoint = 1, name = "Fitness Scores") +
    coord_fixed(ratio = 1, xlim = c(start_codon-1, start_codon + 79)) +
    scale_x_discrete(breaks = seq(start_codon, start_codon + 80, by=1), labels = paste(seq(start_codon, start_codon + 80, by=1), wtseq[((start_codon)-1):(start_codon + 79)])) +
    labs(x = "Position & WT AA", y = "Amino Acid", title = paste("Heatmap of Substitution Fitness Scores in LSG1 (Codons", start_codon, "to", start_codon + 80, ")")) +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1), #rotates x axis so that you can easily see AA and pos
          axis.title = element_text(size=10), #changes size of axes labels
          plot.title = element_text(size=12))        
}


#creates all other heatmaps and outputs them
for (start_codon in seq(80, 560, by = 80)) {
  filename <- paste("heatmap_codons_", start_codon, "_to_", start_codon + 80, ".png", sep = "")
  heatmap <- generate_heatmap(start_codon)  # Generate heatmap
  ggsave(filename, heatmap, width = 12, height = 6)  # Save the generated heatmap
}
