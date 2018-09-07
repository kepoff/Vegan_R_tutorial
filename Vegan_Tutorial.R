rm(list = ls())

# This tutorial will walk you through the basic R skills we cover in class:
#         1.  Overview of diversity and ordination with vegan
#         2.  Loading Data into R - Proper format
#         3.  Common analyses - Rarefactions, Normalization, Diversity Measures
#         4.  Visualizing your diversity metrics and ordinations

# First thing is to install and load the vegan package

# Install vegan, if you haven't already
# install.packages("vegan")

# Load the vegan package (and some other useful stuff)
library(vegan)
library(ggplot2)
library(rgl)
library(plyr)
library(dplyr)


# Now, before we worry about loading our own data or making it look pretty, let's take a basic look
# at how ordinations work in vegan using some simple random data . . . 

################################################################################
#             An overview of metrics and ordination with vegan                 #
################################################################################

# Set random seed so results are reproducible
set.seed(55)

# Generate a random community matrix, similar to the type of data in an OTU Table
prob = c(0.99,rep(0.05,1000)) # sets probability of each number being randomly chosen

community_matrix=matrix(
  sample(0:1000,300,replace=T, prob = prob),nrow=10,
  dimnames=list(paste("Sample_",1:10,sep=""),paste("OTU_",1:30,sep="")))

# Add a treatment vector, assigning samples to one of two groups
treat=c(rep("Treatment_1",5),rep("Treatment_2",5))
trea <- c(rep("Treatment_1",5),rep("Treatment_2",5))

# Take a quick look at your "OTU Table"
community_matrix  # should have 10 samples as rows, and 30 OTUs as columns

# Check to see how even your sampling effort for each community is
barchart(rowSums(community_matrix))
min_depth = min(rowSums(community_matrix)) # this gives us the minimum number of reads in a given sample

# Look at rarefaction (species accumulation) curve
rarecurve(community_matrix, step = 100, sample = min_depth) # vegan's built-in S.A.C. function

# We should normalize our data to account for variable sequencing depth.  One way is to rarefy...
set.seed(55) # set random seed so it's reproducible

rare_community_matrix = rrarefy(community_matrix, min_depth) # randomly subsamples data to given depth and makes new OTU table

# How many species were found in each sample?
specnumber(community_matrix)
specnumber(rare_community_matrix)

# How many species were found in each treatment group?
specnumber(community_matrix, groups = treat)
specnumber(rare_community_matrix, groups = treat)

# Calculate Shannon diversity for each sample
diversity(community_matrix, index = "shannon")
diversity(rare_community_matrix, index = "shannon")
?diversity

# Beta diversity (basic Whittaker index)
beta_div = betadiver(community_matrix, method = "w")
beta_disp = betadisper(beta_div, treat)
plot(beta_disp)

# Run metaMDS(), which automatically calculates a community distance matrix, transforms it,
#                                               and calls monoMDS() which performs the NMDS

example_NMDS = metaMDS(rare_community_matrix, # Our community-by-species matrix
                     k=2) # The number of reduced dimensions



# Take a look at the distances between each pair of communities against their original dissimilarities
# This let's us see how well-preserved community distances are, given that they have been reduced to 2-D

a= stressplot(example_NMDS) # Lots of scatter means poor preservation
?stressplot

?metaMDS()

b=goodness(example_NMDS)
# Plot the ordination
ordiplot(example_NMDS, type = "n") # Sets up the plotting area
orditorp(example_NMDS,display="sites",col="red",air=0.01) # Adds the samples in ordination space
ordiellipse(example_NMDS, groups = treat, label = TRUE) # Calculates the centroid and 95% C.I. of each treatment group

# What you are looking at is a plot of community similarity between your samples that has been reduced to 2D
# Samples that are closer together in the plot have more similar community compositions
# The ellipses overlap, meaning there appears to be little to distinguish communities based on Treatment.
# This is not unexpected, given that our data are fake, and randomly drawn from the same distribution.


# PCOA

princomp(rare_community_matrix)

biplot(princomp(t(rare_community_matrix)))






library(ggplot2)
library(scatterplot3d)


z = example_NMDS$points[,1]
x = example_NMDS$points[,2]
y = example_NMDS$points[,3]

scatterplot3d(x,y,z, box = TRUE, )
?scatterplot3d


################################################################################
#                   Loading data into R - Proper Format                        #
################################################################################

# It is possible to directly import the biom file provided by QIIME, but this depends on the "biomformat" 
# and "Matrix" packages, and versions and dependencies can get confusing, so we will start with a tab-separated OTU table.


# To prepare your biom file you need to convert it to tab-separated format and clean up the first line
# On the command line, you can do the following:

# biom convert --to-tsv --table-type "OTU table" -i otu_table_mc2_w_tax.biom -o otu_table_mc2_w_tax.tsv
# tail -n +2 otu_table_mc2_w_tax.tsv | sed 's/^#OTU ID/OTU_ID/'| head -1 > otu_table.tsv

# those two lines will convert your biom to tsv and clean up the first 2 lines of the file so R can read it easily
# the other file you will need is the sample metadata.  Typically, you will make this in excel, and it has sample IDs as rows
# and information about each sample as columns.

# For this tutorial we're going to start with an example dataset


# To read a tsv file into R we can do the following:
otu_table = read.delim("~/Desktop/vegan_example_otu_table.tsv", sep = "\t", row.names = 1)
metadata = read.delim("~/Desktop/vegan_example_metadata.tsv", sep = "\t", row.names = 1)



# The OTU table from QIIME, by default, has samples/sites as columns and OTUs/Species as rows.
# This is fine, but we will need to transpose it in order to use vegan.
# First, look at the dimensions of the two tables:

dim(otu_table)
dim(metadata)

# There appears to be one additional Sample in the OTU table that isn't in the Metadata... this is the Taxonomic Assignment
# We can pull this off to save it for later:

# Look at the names of the columns in the OTU table (sample names)
names(otu_table) # that last one is the taxonomy from QIIME

# We can build a separate data frame to store our taxonomic assignments
Taxonomy = data.frame(Taxonomy = otu_table$Tax_Name, row.names = rownames(otu_table))

# And now, get rid of them from the otu table for now so we can use vegan
otu_table = otu_table[  ,which(colnames(otu_table) != "Tax_Name")]

# Now, transpose the otu table so rows are samples, and OTUs are columns
otu_table = as.data.frame(t(otu_table))

# Now, let's check the dimensions again
dim(otu_table)
dim(metadata)

# ...And make sure they're in the right order
identical(rownames(otu_table), rownames(metadata)) # if "TRUE", then they are exactly the same, and in the same order

class(rownames(otu_table))

otu_table[order(rownames(otu_table), rownames(metadata)),]

################################################################################
#                Common Analyses - Diversity and Ordination                    #
################################################################################

# Just like the example with random data, we will first take a look at sampling effort
barchart(rowSums(otu_table))

# What is the minimum number of reads in a given sample?
min_depth2 = min(rowSums(otu_table)) # should be 2000 for this data set

# Look at rarefaction curves and draw vertical line at the minimum depth (our probable rarefication level)
rarecurve(otu_table, sample = 2000, step = 100)

# Normalize the data by rarefying
set.seed(1)
rare_otu_table = rrarefy(otu_table, min_depth2)

# How many OTUs are found in each sample?
specnumber(otu_table)
specnumber(rare_otu_table)

# compare before and after rarefication
plot(specnumber(otu_table),specnumber(rare_otu_table)) # there was some loss of species richness during rarification

# Look at it by groups (Ecosystem)
specnumber(otu_table, groups = metadata$Ecosystem) # refers to the metadata table
specnumber(rare_otu_table, groups = metadata$Ecosystem) # refers to the metadata table

# Look at Shannon Diversity
diversity(otu_table, "shannon")
diversity(rare_otu_table, "shannon")

# compare before and after rarefication
plot(diversity(otu_table, "shannon"),
     diversity(rare_otu_table, "shannon"))

# Find the number of samples in which each OTU is found
species_presence = apply(rare_otu_table > 0,2,sum)
barchart(species_presence)

# Calculate a distance matrix for your samples (determine dissimilarity)
# Which method to use?

# We can have some help deciding which distance index to use (Bray-Curtis is the default)
rank_otus = rankindex(metadata$Ecosystem, rare_otu_table, indices = 
            c("bray", "euclid", "manhattan", "horn"), method = "spearman")

print(paste("The highest rank was given by the", names(sort(rank_otus, decreasing = TRUE)[1]), "method."))

# Since rankindex() told us that bray is probably best, we will pass that method to the metaMDS command

# run NMDS (this calculates the distance using "bray" method, transforms the data, and runs isoMDS 20 times to find best
# solution)
MDS = metaMDS(rare_otu_table, distance = "bray", try = 50, trymax = 100)

# Check out the stress plot
stressplot(MDS)  # looks pretty good - has linear fit R2 of 0.946


# Plot it
ordiplot(MDS, type = "p")
ordiellipse(MDS, groups = metadata$Ecosystem, label = TRUE) # Calculates the centroid and 95% C.I. of each treatment group


# Beta Diversity

rare_otu_table.bray = vegdist(rare_otu_table, method = "bray") # Generates the distance matrix, as in metaMDS
beta_div = betadisper(rare_otu_table.bray, group = metadata$Ecosystem) # Calculate homogeneity of multivariate dispersions

boxplot(beta_div) # show boxplot of beta diversity by ecosystem



# Want to make it look pretty?

# We call pull out the cartesian coordinates and build a data frame that is easier to work with
MDS1 = MDS$points[,1] # gives vector of X coordinates
MDS2 = MDS$points[,2] # gives vector of Y coordinates


NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, Ecosystem = metadata$Ecosystem, 
                  Lat = metadata$Latitude, Lon = metadata$Longitude, Host = metadata$Host)




# ggplot2 is its own can of worms, but it can label, color, and make a legend automatically, which is nice
ggplot(NMDS, mapping = aes(x = MDS1, y = MDS2, col = Ecosystem)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  ggtitle("NMDS")

?ggsave()







# PCA

pc = princomp(rare_otu_table, cor = TRUE, scores = TRUE)
plot(pc, type = "lines")
biplot(pc)
biplot(pc, xlim = c(-.2,.2), ylim = c(-.2,.2), cex = .5)

### look at it in 3D

# First, make a vector of colors in your metadata

levels(metadata$Ecosystem) #how many colors do we need?

# copy "Ecosystem" Vector into a new column called "Color"
metadata$Color = metadata$Ecosystem

# Now, we can create a vector of colors based on the ecosystem
metadata$Color
mapvalues(metadata$Color, from = c("Marine","Terrestrial"), to = c("Blue", "Green"))
metadata$Color = as.character(mapvalues(metadata$Color, from = c("Marine","Terrestrial"), to = c("Blue", "Green")))

# which PC scores have most explanatory power?
plot(pc) # <- First three scores explain most of the variance


# plot in 3D using rgl
?plot3d
plot3d(pc$scores[,1:3], col = metadata$Color)


# other visualizations

# heatmap of OTU abundance
?heatmap
heatmap(rare_otu_table, distfun = vegdist)
heatmap(t(rare_otu_table), distfun = vegdist, col = grey.colors(100), 
        ColSideColors = metadata$Color, Rowv = NA,
        main = "Heatmap of OTU Abundance", xlab = "Sample ID")

# Save the image as a file

?png
png("~/Desktop/example_heatmap.png")
heatmap(t(rare_otu_table), distfun = vegdist, col = grey.colors(100), 
        ColSideColors = metadata$Color, Rowv = NA,
        main = "Heatmap of OTU Abundance", xlab = "Sample ID")
dev.off()

# Bar chart of OTU abundance
# using ggplot
rowSums(rare_otu_table) # should all be 2000 because of rarefaction
colSums(rare_otu_table) # OTU abundances across all samples

# make a data frame containing otu abundances and their taxonomic assignments
total_abund = data.frame(Tax_Assignment = Taxonomy, Abundance = colSums(rare_otu_table))

ggplot(total_abund, mapping = aes(x=Taxonomy, y=Abundance)) +
  geom_bar(stat = 'identity')  # draw bars based on values
# that one is awful

ggplot(total_abund, mapping = aes(x=reorder(Taxonomy, 1/Abundance), y=Abundance)) + # coerce order to descending
  geom_bar(stat = 'identity')
# need to be able to read the names

ggplot(total_abund, mapping = aes(x=reorder(Taxonomy, 1/Abundance), y=Abundance)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, face = 'bold.italic', size = 6)) + # fixed text
  labs(x = "NCBI Top Hit") + # fixed axis label
  ggtitle("Overall OTU Abundance") # added title


# try relative abundance instead of rarefaction to normalize otu_table for read depth
norm.otu_table = decostand(otu_table, method = "total") # convert OTU counts to relative abundance values for each sample

MDS.norm = metaMDS(norm.otu_table)
MDS1.norm = MDS.norm$points[,1]
MDS2.norm = MDS.norm$points[,2]

NMDS2 = data.frame(MDS1 = MDS1.norm, MDS2 = MDS2.norm, Ecosystem = metadata$Ecosystem, 
                   Lat = metadata$Latitude, Lon = metadata$Longitude, Host = metadata$Host)

NMDS.plot.2 = ggplot(NMDS2, mapping = aes(MDS1, MDS2, col = Ecosystem)) +
  geom_point() + ggtitle("Normalized by Relative Abundance")

NMDS.plot.1 = ggplot(NMDS, mapping = aes(MDS1, MDS2, col = Ecosystem)) +
  geom_point() + ggtitle("Normalized by rarefaction")

# view plots
NMDS.plot.1
NMDS.plot.2

#save your otu_table as a csv file.
?write.table()
write.csv(otu_table, file = "~/Desktop/otu_table_from_R.csv")

#save your handy ordination data frame
write.csv(NMDS2, file = "~/Desktop/NMDS_results_from_R.csv")

