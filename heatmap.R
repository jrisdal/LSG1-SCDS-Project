##creating a 640x20 tile with geom_tile() to represent the fitness scores of LSG1 AA substitutions

library(ggplot2)
library(colorspace)
library(reshape2)


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


ggplot() +
  geom_tile(data = WT_csv, aes(x = factor(Codon), y = factor(variable), fill= as.numeric(score))) +
  geom_point(data = subset(WT_csv, value == TRUE), aes(x = factor(Codon), y = factor(variable), color = "WT AA")) + # Add dots for WT sequences
  scale_fill_gradient2(low= "brown", mid = "lightyellow", high= "darkgreen", midpoint = 1, name = "Fitness Scores") + #adds diverging color scale
  scale_color_manual(name= "", values= "black") + # Add legend for dots
  coord_fixed(ratio = 1, xlim = c(319, 479)) + #the codons that the figure includes (add 1 to each of them). Change depending on how big you want figure to be
  scale_x_discrete(breaks = c(320, 340, 360, 380, 400, 420, 440, 460, 480)) + #adds the tickmarks (change every time)
  labs(x = "Codon", y = "Amino Acid Substitution", title = "Heatmap of Amino Acid Substitution Fitness Scores in LSG1") +
  theme(axis.text.y = element_text(size = 6))  #changes size of text on y axis 


