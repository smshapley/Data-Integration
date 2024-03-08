##########################################################
#### Log2 Log2 dynamic plots allows for either a static image or a plotly HTML to be created when comparing Log2FCs across datasets
#### For best use: set up ANOVAout to replicate products from https://github.com/edammer/parANOVA pipeline
##########################################################
#Load necessary libraries
library(ggplot2)
library(ggpubr)
library(ggpp) #helps annotate quadrants, provides the stat_quadrant_counts() function

########## Graphing Log2 v Log2 graphics #################
#Set working directory where your ANOVAout (or comparable) files are located
rootdir = ""
setwd(rootdir)

ANOVAout1 <- read.csv(file = "", header = T) #read in your first ANOVAout file == dataset1
colnames(ANOVAout1) <- c("Gene.ID", "UniprotID","F.Value", "FDR","ExpvCTL","diff.ExpvCTL") #name column names consistently
    ## ANOVAout1.sig <- ANOVAout1[ANOVAout1$ExpvCTL<0.05 & ANOVAout1$diff.ExpvCTL>1,] #can subset this initial data frame to ONLY significant and FC > 1 to compare significantly enriched proteins between datasets

ANOVAout2 <- read.csv(file = "", header = T) #read in your second ANOVAout file == dataset2

#create merged dataframe
# Use ANOVA out 1 & 2 for, merge by GeneIDs (can merge by other identifiers as well)
df_log <- merge(ANOVAout1, ANOVAout2, by.x = "UniprotID.1", by.y = "uniprotIDs.2", all=F) 
    #tips: by=0 is by row.names (best option if you're confident in your row names matching between ANOVAouts)
    #Alternatively, change this to X if that's how your data set imports if you didn't define row.names before

#######  Use these two lines of script if your row names are not unique and you need to force uniqueness (look, it happens to the best of us)
# names <- make.unique(df_log$UniprotID.1, sep = "-") #change to whatever you have as your row.name column
# df_log <- cbind(names, df_log)

row.names(df_log) <- df_log$names #if your row names are already unique, then start here!

#define new data frame, make sure you tab complete column names to ensure they actually exist (reduce those syntax errors!)
df.log <- data.frame(x=df_log$diff.ExpvCTL, y=df_log$`ExpvCTL-2`, check.names = T) #CHANGE OUT COLUMN HEADERS FOR COMPARISONS
row.names(df.log) <- names

write.csv(df.log, file = "dataset-overlap/ANOVAout_ExpvCTL2_ExpvCTL1_byUniprotID.csv", row.names = T) #write this to a csv for easy reference later

######### 
## Visualization

# PDF graphics
# pdf("P301LvMock_WTvMock_Log2.pdf")
pdf("ExpvCTL2_ExpvCTL1_byUniprot.pdf")

# Code of the plot
############
#create scatterplot using ggplot2
# x = ExpvCTL2, y = ExpvCTL1
p<- ggplot(df.log, aes(x=df.log$y, y=df.log$x, text = df_log$names)) + #Graph so Set 1 is on the y axis, Set 2 is along the x axis. x and y columns are defined above. Really shouldn't have to change this script 
  geom_point() + 
  theme_pubr() + #minimal theme "publication ready" and easy to modify from here
  stat_quadrant_counts(
    mapping = NULL,
    data = NULL,
    geom = "text_npc",
    position = "identity",
    quadrants = NULL,
    pool.along = c("none", "x", "y", "xy"),
    xintercept = 0,
    yintercept = 0,
    label.x = NULL,
    label.y = NULL,
    digits = 2,
    na.rm = FALSE,
    show.legend = FALSE,
    inherit.aes = TRUE) #allows us to count the occurrences in each quadrant and label on the graph

p+geom_hline(yintercept=0, linetype="dashed", color = "grey") + #feed prior ggplot object into here to modify the h/vline aesthetics
  geom_vline(xintercept = 0, linetype="dashed", color = "grey") + 
  labs(x = "Log2FC(ExpvCTL2)") + #label the x axis
  labs(y = "Log2FC(ExpvCTL1)") #label the y axes

dev.off() #this saves object to PDF
############

# x = ExpvCTL2, y = ExpvCTL1
#Also, can save as variable to feed into plotly
p<- p + geom_hline(yintercept=0, linetype="dashed", color = "grey") + #feed prior ggplot object into here to modify the h/vline aesthetics
  geom_vline(xintercept = 0, linetype="dashed", color = "grey") + 
  labs(x = "Log2FC(ExpvCTL2)") + #label the x axis
  labs(y = "Log2FC(ExpvCTL1)") #label the y axes

Log2 <- ggplotly(p, tooltip = "text") #load prior ggplot object with the new labels into plotly. This tooltip="text" allows us to hover for row.names text
tempfilename="scatterplot_Log2ExpvCTL2andExpvCTL1-byUniprot.html" 
htmlwidgets::saveWidget(Log2, tempfilename) #save to html for plotly object

