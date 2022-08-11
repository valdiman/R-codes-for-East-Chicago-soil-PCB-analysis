# R codes for creating figures and maps for East Chicago PCB Soil paper

# Packages and libraries needed -------------------------------------------------------------------
# Install packages
install.packages("readxl")
install.packages("ggplot2")
install.packages("ggmap")
install.packages("ggrepel")
install.packages("scales")
install.packages("ggpubr")
install.packages("reshape2")
install.packages("egg")
install.packages("ggfortify")
install.packages('pangaear')
install.packages("tidyverse")

# Libraries
library(pangaear) # Read data from Pangaea
library(ggplot2) # make_bbox
library(ggmap) # make_bbox
library(ggrepel) #geom_label_repel
library(scales) # comma_format
library(ggpubr) # ggarrange
library(egg) # ggarrange too
library(reshape2) # melt
library(ggfortify) # PCA analysis
library(tidyverse) # data manipulation

# Read data from Pangaea and format data ------------------------------

# Read data from Pangaea repository
# Citation:
# Martinez, Andres; Hua, Jason; Haque, Ezazul; Thorne, Peter;
# Hornbuckle, Keri C (2022): Dataset of concentrations of
# individual polychlorinated biphenyl congeners and total
# organic carbon in soils from East Chicago, Indiana, USA in
# 2017/2018. PANGAEA, https://doi.pangaea.de/10.1594/PANGAEA.941894

# Set cache path to the project folder
pg_cache$cache_path_set(full_path = "/Users/andres/OneDrive - University of Iowa/work/ISRP/Project6/Soil/R/ECSS")

# Download original datasets from Pangaea
s.0 <- pg_data(doi = '10.1594/PANGAEA.941881') # soil dataset
b.0 <- pg_data(doi = '10.1594/PANGAEA.941686') # blank dataset
toc.0 <- pg_data(doi = '10.1594/PANGAEA.941884') # toc dataset

# Obtain just concentrations from Pangaea dataset
s <- data.frame(s.0[[1]]$data) # ng/g
b <- data.frame(b.0[[1]]$data)
toc <- data.frame(toc.0[[1]]$data)

# Format concentration data for plotting
# Remove metadata from soil dataset
s.1 <- subset(s, select = -c(Event:Distance..m...Distance.to.Fork.))
# Rename PCB congener in columns in s.1
names(s.1) <- substr(names(s.1), 18, 36) # Leave PCB names
names(s.1) <- sub("\\.", "+", names(s.1)) # change "." to "+"
names(s.1) <- sub("\\.", "+", names(s.1)) # change "." to "+"
names(s.1) <- sub("\\.", "+", names(s.1)) # change "." to "+"
names(s.1) <- sub("\\.", "+", names(s.1)) # change "." to "+"
names(s.1) <- gsub(".{1}$", "", names(s.1)) # Remove last "+"
# Collect metadata
meta.s <- subset(s, select = c(Sample.label, Longitude, Latitude,
                               Distance..m...Distance.to.Fork.))

# Remove metadata from blank dataset
b.1 <- subset(b, select = -c(Event:Date.time.end))
# Rename PCB congener in columns in b.1
names(b.1) <- substr(names(b.1), 16, 36) # Leave PCB names
names(b.1) <- sub("\\.", "+", names(b.1)) # change "." to "+"
names(b.1) <- sub("\\.", "+", names(b.1)) # change "." to "+"
names(b.1) <- sub("\\.", "+", names(b.1)) # change "." to "+"
names(b.1) <- gsub(".{1}$", "", names(b.1)) # Remove last "+"

# Create LOQ, i.e., upper 95 CI% (mean + 1.96*sd/(n)^0.5)
b.1[b.1 == 0] <- NA # Not the best way, but zeros need to be removed for lo10 transformation
log10b.1 <- log10(b.1) # log10 blank data
loq <- colMeans(log10b.1, na.rm = TRUE) + 1.96*sapply(log10b.1, sd, na.rm = TRUE)/sqrt(8)
loq <- data.frame(t(loq))
names(loq) <- sub("\\.", "+", names(loq)) # change "."
names(loq) <- sub("\\.", "+", names(loq)) # change "."
names(loq) <- sub("\\.", "+", names(loq)) # change "."

# Make comparison between s.1 and loq
# If s.1 > loq, then s.1, if not 0
# Create matrix to storage s.1 or loq values in s.2
s.2 <- matrix(NA, nrow = dim(s.1)[1], ncol = dim(s.1)[2])
# Do comparison
for(i in 1:dim(s.1)[1]) {
  for(j in 1:dim(s.1)[2]) {
    if (log10(s.1[i, j]) > loq[j]) {
      s.2[i, j] <- s.1[i, j]
    } else {
      s.2[i, j] <- 0
    }
  }
}

colnames(s.2) <- colnames(s.1) # add PCB names to columns.

# Remove congeners with < 75% detection frequency
s.2 <- s.2[, colMeans(s.2 > 0) >= 0.75]

# Remove the same congeners in loq
loq.2 <- loq %>% select(colnames(s.2))

# Add loq value = (loq/(2)^0.5)), where zeros are in s.2
# Create new matrix to storage final concentrations
s.3 <- data.frame(matrix(NA, nrow = dim(s.2)[1], ncol = dim(s.2)[2]))

for(i in 1:dim(s.3)[1]) {
  for(j in 1:dim(s.3)[2]) {
    if (s.2[i, j] == 0) {
      s.3[i, j] <- 10^(loq.2[j])/sqrt(2)
    } else {
      s.3[i, j] <- s.2[i, j]
    }
  }
}

# Final dataset for concentrations, s.3
colnames(s.3) <- colnames(s.2) # add PCB names to columns.

# Maps and figures --------------------------------------------------------

# Map with ggmap
# Create a square map around samples
IN.box <- make_bbox(lon = s$Longitude, lat = s$Latitude, f = 0.6)
in.map <- get_stamenmap(bbox = IN.box, zoom = 13)
# Change the location names
s$Sample.label <- gsub('EC-', 's', s$Sample.label)

# Map sampling locations (Figure S1)
ggmap(in.map) +
  geom_point(data = s, aes(x = Longitude, y = Latitude)) +
  geom_label_repel(aes(x = Longitude, y = Latitude,
                       label = Sample.label),
                   data = s, size = 3, color = "red",
                   box.padding = unit(0.5, "lines"),
                   max.overlaps = getOption("ggrepel.max.overlaps",
                                            default = 30)) +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(title = "East Chicago soil locations")

# Sum individual PCB congeners to yield total PCB (tPCB)
tPCB <- rowSums(s.3)

# Add coordinates from sample locations
tPCB <- data.frame(cbind(tPCB, s$Longitude, s$Latitude))

# Rename column names
colnames(tPCB) <- c("tPCB", "long", "lat")

# Map tPCB (Figure 1). Includes: The Fork, Lake Michigan and IHSC.
ggmap(in.map) +
  geom_point(data = tPCB, aes(x = long, y = lat,
                              size = tPCB)) +
  geom_label_repel(aes(x = long, y = lat,
                       label = formatC(signif(tPCB, digits = 2))),
                   data = tPCB, size = 3, color = "red",
                   box.padding = unit(0.5, "lines"),
                   max.overlaps = getOption("ggrepel.max.overlaps",
                                            default = 28)) +
  scale_size_area(breaks = c(250, 500, 750, 1000, 1500),
                  labels = c('0-250', '251-500', '501-750',
                             '751-1000', '1001-1500'),
                  name = expression(Sigma*"PCB (ng/g DW)")) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = c(0.133, 0.81),
        legend.title = element_text(size = 10),
        legend.background = element_rect(linetype = "solid",
                                         colour ="black"),
        legend.text = element_text(size = 8)) +
  annotate('text', x = -87.471, y = 41.646,
           label = 'The Fork', colour = 'black', size = 4,
           fontface = 2) +
  annotate('text', x = -87.415, y = 41.65,
           label = 'Lake Michigan', colour = 'black', size = 4,
           fontface = 2) +
  annotate('text', x = -87.438, y = 41.668,
           label = 'IHSC', colour = 'black', size = 4,
           fontface = 2)

# Map individual congeners. E.g., PCB 8 (Figure S2).
ggmap(in.map) +
  geom_point(data = s, aes(x = Longitude, y = Latitude,
                           size = PCB.soil..ng.g...PCB8.)) +
  geom_label_repel(aes(x = Longitude, y = Latitude,
                       label = formatC(signif(PCB.soil..ng.g...PCB8.,
                                              digits = 2))),
                   data = s, size = 3, color = "red",
                   box.padding = unit(0.5, "lines"),
                   max.overlaps = getOption("ggrepel.max.overlaps",
                                            default = 28)) +
  scale_size_area(breaks = c(0.3, 2, 4),
                  labels = c('0.10-2.0', '2.1-4.0', '4.1-7.0'),
                  name = expression("PCB 8 (ng/g DW)")) +
  theme(legend.position = c(0.133, 0.83),
        legend.title = element_text(size = 10),
        legend.background = element_rect(linetype = "solid",
                                         colour ="black"),
        legend.text = element_text(size = 9)) +
  annotate('text', x = -87.471, y = 41.646,
           label = 'The Fork', colour = 'black', size = 4,
           fontface = 2) +
  annotate('text', x = -87.415, y = 41.65,
           label = 'Lake Michigan', colour = 'black', size = 4,
           fontface = 2) +
  annotate('text', x = -87.438, y = 41.668,
           label = 'IHSC', colour = 'black', size = 4,
           fontface = 2)

# Map TOC (Figure S3)

ggmap(in.map) +
  geom_point(data = toc, aes(x = Longitude, y = Latitude,
                           size = TOC....), color = "black") +
  geom_label_repel(aes(x = Longitude, y = Latitude,
                       label = formatC(signif(TOC...., digits = 2))),
                   data = toc, size = 3, color = "red",
                   box.padding = unit(0.5, "lines"),
                   max.overlaps = getOption("ggrepel.max.overlaps",
                                            default = 28)) +
  scale_size(name = bquote(atop("TOC (%)")),
             range = c(1, 6)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = c(0.08, 0.84),
        legend.title = element_text(size = 10),
        legend.background = element_rect(linetype = "solid",
                                         colour ="black"),
        legend.text = element_text(size = 10)) +
  annotate('text', x = -87.471, y = 41.646,
           label = 'The Fork', colour = 'black', size = 4,
           fontface = 2) +
  annotate('text', x = -87.415, y = 41.65,
           label = 'Lake Michigan', colour = 'black', size = 4,
           fontface = 2) +
  annotate('text', x = -87.438, y = 41.668,
           label = 'IHSC', colour = 'black', size = 4,
           fontface = 2)

# Spatial plots (Figure 2) ------------------------------------------------

# Organize site location name for plotting
site <- c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9",
          "s10", "s11", "s12", "s13", "s14", "s15", "s16", "s17",
          "s18", "s19", "s20", "s21", "s22", "s23", "s24", "s25",
          "s26", "s27", "s28", "s29", "s30", "s31", "s32", "s33")

# Spatial plot (Figure 2)
# Individual PCB congener plots + tPCB
# Data from s
# PCB 8 plot
p8 <- ggplot(s, aes(y = PCB.soil..ng.g...PCB8., x = factor(Sample.label,
                                         levels = site))) + 
  geom_bar(stat = 'identity', width = 0.8, fill = "black") +
  scale_y_continuous(labels = comma_format(big.mark = ".",
                                           decimal.mark = ",")) +
  theme_bw() +
  theme(aspect.ratio = 6/35) +
  ylab(expression("PCB 8 (ng/g DW)")) +
  xlab(expression("")) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = c(1.5, 3.5, 4.5, 5.5, 7.5, 11.5, 12.5, 13.5,
                            14.5, 15.5, 17.5, 18.5, 19.5, 21.5, 23.5,
                            25.5, 29.5, 32.5),
             color = "#999999")

# PCB 11 plot
p11 <- ggplot(s, aes(y = PCB.soil..ng.g...PCB11., x = factor(Sample.label,
                                           levels = site))) + 
  geom_bar(stat = 'identity', width = 0.8, fill = "black") +
  scale_y_continuous(labels = comma_format(big.mark = ".",
                                           decimal.mark = ",")) +  
  theme_bw() +
  theme(aspect.ratio = 6/35) +
  ylab(expression("PCB 11 (ng/g DW)")) +
  xlab(expression("")) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = c(1.5, 3.5, 4.5, 5.5, 7.5, 11.5, 12.5, 13.5,
                            14.5, 15.5, 17.5, 18.5, 19.5, 21.5, 23.5,
                            25.5, 29.5, 32.5),
             color = "#999999")

# PCB 52 plot
p52 <- ggplot(s, aes(y = PCB.soil..ng.g...PCB52., x = factor(Sample.label,
                                           levels = site))) + 
  geom_bar(stat = 'identity', width = 0.8, fill = "black") +
  scale_y_continuous(labels = comma_format(big.mark = ".",
                                           decimal.mark = ",")) +
  theme_bw() +
  theme(aspect.ratio = 6/35) +
  ylab(expression("PCB 52 (ng/g DW)")) +
  xlab(expression("")) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = c(1.5, 3.5, 4.5, 5.5, 7.5, 11.5, 12.5, 13.5,
                            14.5, 15.5, 17.5, 18.5, 19.5, 21.5, 23.5,
                            25.5, 29.5, 32.5),
             color = "#999999")

# PCB 136 plot
p136 <- ggplot(s, aes(y = PCB.soil..ng.g...PCB136., x = factor(Sample.label,
                                             levels = site))) + 
  geom_bar(stat = 'identity', width = 0.8, fill = "black") +
  scale_y_continuous(labels = comma_format(big.mark = ".",
                                           decimal.mark = ",")) +
  theme_bw() +
  theme(aspect.ratio = 6/35) +
  ylab(expression("PCB 136 (ng/g DW)")) +
  xlab(expression("")) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = c(1.5, 3.5, 4.5, 5.5, 7.5, 11.5, 12.5, 13.5,
                            14.5, 15.5, 17.5, 18.5, 19.5, 21.5, 23.5,
                            25.5, 29.5, 32.5),
             color = "#999999") 

# PCB 206 plot
p206 <- ggplot(s, aes(y = PCB.soil..ng.g...PCB206., x = factor(Sample.label,
                                             levels = site))) + 
  geom_bar(stat = 'identity', width = 0.8, fill = "black") +
  scale_y_continuous(labels = comma_format(big.mark = ".",
                                           decimal.mark = ",")) +
  theme_bw() +
  theme(aspect.ratio = 6/35) +
  ylab(expression("PCB 206 (ng/g DW)")) +
  xlab(expression("")) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = c(1.5, 3.5, 4.5, 5.5, 7.5, 11.5, 12.5, 13.5,
                            14.5, 15.5, 17.5, 18.5, 19.5, 21.5, 23.5,
                            25.5, 29.5, 32.5),
             color = "#999999")

# PCB 209 plot
p209 <- ggplot(s, aes(y = PCB.soil..ng.g...PCB209., x = factor(Sample.label,
                                             levels = site))) + 
  geom_bar(stat = 'identity', width = 0.8, fill = "black") +
  scale_y_continuous(labels = comma_format(big.mark = ".",
                                           decimal.mark = ",")) +
  theme_bw() +
  theme(aspect.ratio = 6/35) +
  ylab(expression("PCB 209 (ng/g DW)")) +
  xlab(expression("")) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = c(1.5, 3.5, 4.5, 5.5, 7.5, 11.5, 12.5, 13.5,
                            14.5, 15.5, 17.5, 18.5, 19.5, 21.5, 23.5,
                            25.5, 29.5, 32.5),
             color = "#999999")

# tPCB plot
ptPCB <- ggplot(s, aes(y = rowSums(s.3), x = factor(Sample.label,
                                                    levels = site))) + 
  geom_bar(stat = 'identity', width = 0.8, fill = "black") +
  theme_bw() +
  theme(aspect.ratio = 6/35) +
  ylab(expression(Sigma*"PCB (ng/g DW)")) +
  xlab(expression("")) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 7,
                                   angle = 60, hjust = 1),
        axis.title.x = element_text(face = "bold", size = 7)) +
  geom_vline(xintercept = c(1.5, 3.5, 4.5, 5.5, 7.5, 11.5, 12.5, 13.5,
                            14.5, 15.5, 17.5, 18.5, 19.5, 21.5, 23.5,
                            25.5, 29.5, 32.5),
             color = "#999999")

# Plot all plots in one panel (Figure 2)
ggarrange(p8, p11, p52, p136, p206, p209, ptPCB, ncol=1)

# Total PCB concentration comparison plot (Figure 3) ----------------------

# Create matrix to storage external data
s.tpcb <- matrix(NA, nrow = 10, ncol = 73)
# Add column names for s.2
colnames(s.tpcb) <- paste0("s", 0:(ncol(s.tpcb)-1))
colnames(s.tpcb)[1] <- c('location')

# Add location name to rows to s.tpcb
s.tpcb[,1] <- c("(k) East Chicago, USA", "(a) Cedar Rapids, USA",
             "(d) USA background", "(e) Global background",
             "(c) London, UK", "(b) China", "(f) Gadsden, USA",
             "(g) New Bedford, USA",
             "(h) Wenling, China", "(j) Kalamazoo, USA")

# Fill matrix s.2 with total PCB data (ng/g) from EC and other studies
# Row 1 East Chicago data from this study
s.tpcb[1, 2:34] <- tPCB$tPCB

# Row 2 Cedar Rapids, IA data
# Ref: 
# Martinez, A., N. R. Erdman, Z. L. Rodenburg, P. M. Eastling
# and K. C. Hornbuckle (2012). "Spatial distribution of
# chlordanes and PCB congeners in soil in Cedar Rapids,
# Iowa, USA." Environmental Pollution 161: 222-228.
# https://doi.org/10.1016/j.envpol.2011.10.028 

s.tpcb[2, 2:65] <- c(6.253491724, 16.80387245, 12.63058738, 17.70738503,
                  5.059038158, 15.11833429, 15.30164564, 4.4212167,
                  8.133366499, 15.80452114, 10.77403669, 4.561438133,
                  2.966628839, 50.67778504, 10.31410547, 9.919715784,
                  15.88149424, 17.31904254, 109.8589162, 14.28944635,
                  25.20520682, 101.7270658, 112.2568732, 15.08458713,
                  16.78324533, 8.887912669, 111.1558474, 1237.792597,
                  22.32334494, 20.00868752, 48.11045486, 32.24177906,
                  14.41009595, 13.85213435, 53.55265063, 63.87656227,
                  46.99782271, 27.43837787, 19.43270661, 17.69035474,
                  30.59548161, 30.55263431, 27.20893854, 46.69509645,
                  58.98510302, 12.80152206, 18.77768408, 24.05794805,
                  8.294087523, 6.795078005, 46.45813946, 30.0188002,
                  238.1083633, 39.58545627, 28.24760348, 29.91904074,
                  41.72994497, 6.313088911, 17.26741528, 299.9327233,
                  16.38622096, 82.68447304, 22.93793716, 44.96459011)

# Row 3 USA background data
# Ref:
# USEPA (2007). Pilot Survey of Levels of Polychlorinated
# Dibenzo-P-Dioxins (PCDDs), Polychlorinated Dibenzofurans
# (PCDFs), Polychlorinated Biphenyls (PCB) and Mercury in
# Rural Soils of the U.S. Washington, DC.

s.tpcb[3, 2:25] <- c(1.366, 0.475, 2.604, 15.7, 2.037, 0.358, 1.115,
                  2.464, 4.028, 3.023, 1.543, 1.019, 0.845, 0.303,
                  4.93, 4.954, 24.57, 0.713, 0.57, 0.509, 3.274,
                  1.3, 2.419, 1.224)

# Row 4 Global background data
# Ref:
# Meijer, S. N., W. A. Ockenden, A. Sweetman, K. Breivik,
# J. O. Grimalt and K. C. Jones (2003). "Global distribution
# and budget of PCBs and HCB in background surface soils:
# Implications or sources and environmental processes."
# Environmental Science & Technology 37(4): 667-672.
# https://pubs.acs.org/doi/10.1021/es025809l

s.tpcb[4, 2:22] <- c(0.026, 96.9, 5.9661, 3.3517, 0.6442, 8.7147,
                  1.5107, 3.165, 3.165, 0.0472, 9.5003, 16.0376,
                  24.9879, 7.4148, 7.3494, 48.2592, 5.5822,
                  0.5264, 0.3557, 37.121, 6.5972)

# Row 5 London, UK, data
# Ref:
# Vane, C. H., A. W. Kim, D. J. Beriro, M. R. Cave, K. Knights,
# V. Moss-Hayes and P. C. Nathanail (2014). "Polycyclic
# aromatic hydrocarbons (PAH) and polychlorinated biphenyls
# (PCB) in urban soils of Greater London, UK." Applied Geochemistry
# 51: 303-314.
# https://doi.org/10.1016/j.apgeochem.2014.09.013

s.tpcb[5, 2:73] <- c(28.85, 26.58, 36.03, 9.35,	16.89,	43.79,
                  82.53, 20.46,	55.69,	779.69,	92.26,	12.78,
                  15.51,	12.3,	28.76,	43.96,	42.14,	19.29,
                  103.19,	10.4,	21.66,	53.14,	47.13,	19.23,
                  173.2,	95.18, 161.38,	84.82,	125.38,	31.2,
                  32.55,	41.77, 77.97,	52.15,	26.04,	33.82,
                  15.61,	95.01, 106.68, 18.28,	23.17,	27.86,
                  44.89,	22.17, 50.59,	32.12, 205.45,	19.69,
                  132.71,	34.27, 60.69,	16.62, 15.18,	24.78,
                  11.11,	33.58, 17.81,	79.75, 40.37,	18.47,
                  59.92,	56.21, 174.09, 40.01, 2645.84, 37.88,
                  35.06,	65.9,	79.39, 128.19, 45.33, 14.62)

# Row 6 China data
# Ref:
# Yu, H. Y., Y. F. Liu, X. Q. Shu, L. M. Ma and Y. W. Pan (2020).
# "Assessment of the spatial distribution of organochlorine
# pesticides (OCPs) and polychlorinated biphenyls (PCBs) in urban
# soil of China." Chemosphere 243.
# https://doi.org/10.1016/j.chemosphere.2019.125392

s.tpcb[6, 2:25] <- c(537.64,	31.35, 113.26, 22.68, 10.16,
                  27.22, 32.38,	2.42,	38.77, 1.78, 2.74,
                  0.3, 1.3, 0.84, 281.41, 1.14, 1.9,
                  1.28,	2.06, 0.44, 14.72, 5.63, 0.44, 0.72)

# Row 7 New Bedford, USA, data
# Ref:
# Vorhees, D. J., A. C. Cullen and L. M. Altshul (1999).
# "Polychlorinated Biphenyls in House Dust and Yard Soil
# near a Superfund Site." Environmental Science & Technology
# 33(13): 2151-2156.
# https://pubs.acs.org/doi/10.1021/es9812709
# Mean value only

s.tpcb[7,2] <- c(200)

# Row 8 Gadsden, USA, data
# Ref:
# Carey, A. E., P. Douglas, H. Tai, W. G. Mitchell
# and G. B. Wiersma (1979). "Pesticide-residue
# concentrations in soils of 5 United States cities,
# 1971 - Urban Soils Monitoring Program." Pesticides
# Monitoring Journal 13(1): 17-22.
# Mean value only

s.tpcb[8,2] <- c(220)

# Row 9 Wenling, China, data
# Ref:
# Tang, X., C. Shen, D. Shi, S. A. Cheema,
# M. I. Khan, C. Zhang and Y. Chen (2010).
# "Heavy metal and persistent organic compound
# contamination in soil from Wenling: An emerging e-waste
# recycling city in Taizhou area, China." Journal of
# Hazardous Materials 173(1): 653-660.
# https://doi.org/10.1016/j.jhazmat.2009.08.134
# Mean value only

s.tpcb[9,2] <- c(123467.75)

# Row 10 Kalamazoo, USA, data
# Ref:
# Blankenship, A. L., M. J. Zwiernik, K. K. Coady,
# D. P. Kay, J. L. Newsted, K. Strause, C. Park,
# P. W. Bradley, A. M. Neigh, S. D. Millsap, P. D. Jones
# and J. P. Giesy (2005). "Differential accumulation
# of polychlorinated biphenyl congeners in the terrestrial
# food web of the Kalamazoo River superfund site, Michigan.
# " Environmental Science & Technology 39(16): 5954-5963.
# https://pubs.acs.org/doi/10.1021/es048317c
# Mean value only

s.tpcb[10,2] <- c(6530)

# Modify and adjust format
# Transform matrix to data.frame
s.tpcb <- data.frame(s.tpcb)

# Organize s.tpcb in a long data.frame format
s.tpcb <- melt(s.tpcb, id.vars = "location")

# transform value into numeric value
s.tpcb$value <- as.numeric(s.tpcb$value)

# Plot Figure 3
ggplot(s.tpcb, aes(x = location, y = value)) +
  geom_boxplot(width = 0.6, outlier.colour = "black",
               outlier.shape = 1) +
  scale_y_log10() +
  theme_bw() +
  theme(aspect.ratio = 25/20) +
  ylab(expression(bold(Sigma*"PCB soil concentration (ng/g DW)"))) +
  xlab(expression(bold(""))) +
  theme(axis.text.y = element_text(face = "bold", size = 9,
                                   color = "black"),
        axis.title.y = element_text(face = "bold", size = 9,
                                    color = "black")) +
  theme(axis.text.x = element_text(face = "bold", size = 8,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.6, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm")) +
  geom_bracket(xmin = 1, xmax = 5, y.position = 3.7,
               label = " ", size = 0.8, color = "blue") +
  geom_bracket(xmin = 9.5, xmax = 10.5, y.position = 3.7,
               label = " ", size = 0.8, color = "blue") +
  geom_bracket(xmin = 3, xmax = 10, y.position = 4,
               label = " * Tukey test, p < 0.05", label.size = 3.5,
               size = 0.8, color = "blue")

# Individual PCB boxplot (Figure 4) ---------------------------------------
# need to change the column names!

ggplot(stack(s.3), aes(x = ind, y = values)) +
  scale_y_log10(limits = c(0.005, 300),
                breaks = c(0.01, 0.1, 1, 10, 100),
                labels = c(0.01, 0.1, 1, 10, 100)) +
  geom_boxplot(width = 0.6, outlier.colour = "black",
               outlier.shape = 1) +
  theme_bw() +
  theme(aspect.ratio = 25/135) +
  xlab(expression("")) +
  ylab(expression(bold("Soil concentration (ng/g DW)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10,
                                   color = "black"),
        axis.title.y = element_text(face = "bold", size = 10,
                                    color = "black")) +
  theme(axis.text.x = element_text(face = "bold", size = 5,
                                   angle = 60, hjust = 1,
                                   color = "black"),
        axis.title.x = element_text(face = "bold", size = 8)) +
  theme(axis.ticks = element_line(size = 0.6, color = "black"), 
        axis.ticks.length = unit(0.2, "cm")) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# TOC vs. PCB concentrations plots (Figure 5) -----------------------

# Include TOC data to s.3
s.toc <- cbind(toc$TOC...., s.3)
# Remove samples (row) where toc was not measured (NA)
s.toc <- na.omit(s.toc)
# Change name in s.toc to toc
colnames(s.toc)[1] <- "toc"

# Perform linear regressions, individual PCB congeners and tPCB vs. TOC
# PCB 8 vs. TOC
a.pcb8 <- log10(s.toc$PCB8)
fit1 <- lm(a.pcb8 ~ s.toc$toc)
# Plot PCB 8 vs. TOC linear regression
p8.toc <- ggplot(s.toc, aes(y = PCB8, x = toc)) +
  geom_point(shape = 21, colour = "black", fill = "white",
             size = 1.7, stroke = 0.8) +
  stat_smooth(method = "lm", col = "black", se = FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(aspect.ratio = 4/6) +
  xlab(expression(bold("TOC %"))) + 
  ylab(expression(bold("PCB 8 (ng/g DW)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)) +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  annotate("text", x = 18, y = 0.15,
           label = paste("R2 = ", signif(summary(fit1)$adj.r.squared,
                                         2),
                         " m =", signif(fit1$coef[[2]], 2),
                         " p =", signif(summary(fit1)$coef[2,4],
                                        2)),
           size = 2, fontface = 2) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# PCB 11 vs. TOC
a.pcb11 <- log10(s.toc$PCB11)
fit11 <- lm(a.pcb11 ~ s.toc$toc)
# Plot PCB 11 vs. TOC linear regression
p11.toc <- ggplot(s.toc, aes(y = PCB11, x = toc)) +
  geom_point(shape = 21, colour = "black", fill = "white",
             size = 1.7, stroke = 0.8) +
  stat_smooth(method = "lm", col = "black", se = FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(aspect.ratio = 4/6) +
  xlab(expression(bold("TOC %"))) + 
  ylab(expression(bold("PCB 11 (ng/g DW)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)) +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  annotate("text", x = 18, y = 0.1,
           label = paste("R2 = ", signif(summary(fit11)$adj.r.squared,
                                         2),
                         " m =", signif(fit11$coef[[2]], 2),
                         " p =", signif(summary(fit11)$coef[2,4],
                                        2)),
           size = 2, fontface = 2) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# PCB 52 vs. TOC
a.pcb52 <- log10(s.toc$PCB52)
fit52 <- lm(a.pcb52 ~ s.toc$toc)
# Plot PCB 52 vs. TOC linear regression
p52.toc <- ggplot(s.toc, aes(y = PCB52, x = toc)) +
  geom_point(shape = 21, colour = "black", fill = "white",
             size = 1.7, stroke = 0.8) +
  stat_smooth(method = "lm", col = "black", se = FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(aspect.ratio = 4/6) +
  xlab(expression(bold("TOC %"))) + 
  ylab(expression(bold("PCB 52 (ng/g DW)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)) +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  annotate("text", x = 18, y = 1.7,
           label = paste("R2 = ", signif(summary(fit52)$adj.r.squared,
                                         2),
                         " m =", signif(fit52$coef[[2]], 2),
                         " p =", signif(summary(fit52)$coef[2,4],
                                        2)),
           size = 2, fontface = 2) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# PCB 136 vs. TOC
a.pcb136 <- log10(s.toc$PCB136)
fit136 <- lm(a.pcb136 ~ s.toc$toc)
# Plot PCB 136 vs. TOC linear regression
p136.toc <- ggplot(s.toc, aes(y = PCB136, x = toc)) +
  geom_point(shape = 21, colour = "black", fill = "white",
             size = 1.7, stroke = 0.8) +
  stat_smooth(method = "lm", col = "black", se = FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(aspect.ratio = 4/6) +
  xlab(expression(bold("TOC %"))) + 
  ylab(expression(bold("PCB 136 (ng/g DW)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)) +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  annotate("text", x = 18, y = 0.09,
           label = paste("R2 = ", signif(summary(fit136)$adj.r.squared,
                                         2),
                         " m =", signif(fit136$coef[[2]], 2),
                         " p =", signif(summary(fit136)$coef[2,4],
                                        2)),
           size = 2, fontface = 2) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# PCB 209 vs. TOC
a.pcb209 <- log10(s.toc$PCB209)
fit209 <- lm(a.pcb209 ~ s.toc$toc)
# Plot PCB 209 vs. TOC linear regression
p209.toc <- ggplot(s.toc, aes(y = PCB209, x = toc)) +
  geom_point(shape = 21, colour = "black", fill = "white",
             size = 1.7, stroke = 0.8) +
  stat_smooth(method = "lm", col = "black", se = FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(aspect.ratio = 4/6) +
  xlab(expression(bold("TOC %"))) + 
  ylab(expression(bold("PCB 209 (ng/g DW)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)) +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  annotate("text", x = 18, y = 0.08,
           label = paste("R2 = ", signif(summary(fit209)$adj.r.squared,
                                         2),
                         " m =", signif(fit209$coef[[2]], 2),
                         " p =", signif(summary(fit209)$coef[2,4],
                                        2)),
           size = 2, fontface = 2) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# Remove metadata
s.toc.2 <- subset(s.toc, select = -c(toc))
# Perform linear regression for tPCB vs. TOC
a <- log10(rowSums(s.toc.2))
fittPCB <- lm(a ~ s.toc$toc)
# Plot tPCB vs. TOC
ptPCB.toc <- ggplot(s.toc, aes(x = toc, y = rowSums(s.toc.2))) +
  geom_point(shape = 21, colour = "black", fill = "white",
             size = 1.7, stroke = 0.8) +
  stat_smooth(method = "lm", col = "black", se = FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(aspect.ratio = 4/6) +
  xlab(expression(bold("TOC %"))) + 
  ylab(expression(bold(Sigma*"PCB (ng/g DW)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)) +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  annotate("text", x = 18, y = 22,
           label = paste("R2 = ", signif(summary(fittPCB)$adj.r.squared,
                                         2),
                         " m =", signif(fittPCB$coef[[2]], 2),
                         " p =", signif(summary(fittPCB)$coef[2,4],
                                        2)),
           size = 2, fontface = 2) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# Plot all plots in one panel (Figure 5)
ggarrange(p8.toc, p11.toc, p52.toc, p136.toc, p209.toc, ptPCB.toc, ncol=2)

# Distance IHSC (The Fork) vs PCB concentration plots (Figure 6) ---------------------------

# Change name to distance in s
colnames(s)[8] <- "distanceFork"
d <- s$distanceFork # m
# PCB 56 vs distance to IHSC (The Fork)
a.pcb56 <- log10(s.3$PCB56)
fit56 <- lm(a.pcb56 ~ d)
p56.dis <- ggplot(s.3, aes(x = d, y = PCB56)) +
  geom_point(shape = 21, colour = "black", fill = "white",
             size = 1.7, stroke = 0.8) +
  stat_smooth(method = "lm", col = "black", se = FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(aspect.ratio = 4/6) +
  xlab(expression(bold("distance (m)"))) + 
  ylab(expression(bold("PCB 56 (ng/g DW)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)) +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  annotate("text", x = 2500, y = 140,
           label = paste("R2 = ", signif(summary(fit56)$adj.r.squared,
                                         2),
                         " m =", signif(fit56$coef[[2]], 2),
                         " p =", signif(summary(fit56)$coef[2,4],
                                        2)),
           size = 2, fontface = 2) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))
  
# PCB 66 vs distance to IHSC (The Fork)
a.pcb66 <- log10(s.3$PCB66)
fit66 <- lm(a.pcb66 ~ d)
p66.dis <- ggplot(s.3, aes(y = PCB66, x = d)) +
  geom_point(shape = 21, colour = "black", fill = "white",
             size = 1.7, stroke = 0.8) +
  stat_smooth(method = "lm", col = "black", se = FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(aspect.ratio = 4/6) +
  xlab(expression(bold("distance (m)"))) + 
  ylab(expression(bold("PCB 66 (ng/g DW)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)) +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  annotate("text", x = 2500, y = 140,
           label = paste("R2 = ", signif(summary(fit66)$adj.r.squared,
                                         2),
                         " m =", signif(fit66$coef[[2]], 2),
                         " p =", signif(summary(fit66)$coef[2,4],
                                        2)),
           size = 2, fontface = 2) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))
  
# PCB 110 vs distance to IHSC (The Fork)
a.pcb110 <- log10(s.3$PCB110)
fit110 <- lm(a.pcb110 ~ d)
p110.dis <- ggplot(s.3, aes(y = PCB110, x = d)) +
  geom_point(shape = 21, colour = "black", fill = "white",
             size = 1.7, stroke = 0.8) +
  stat_smooth(method = "lm", col = "black", se = FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(aspect.ratio = 4/6) +
  xlab(expression(bold("distance (m)"))) + 
  ylab(expression(bold("PCB 110 (ng/g DW)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)) +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  annotate("text", x = 2500, y = 80,
           label = paste("R2 = ", signif(summary(fit110)$adj.r.squared,
                                         2),
                         " m =", signif(fit110$coef[[2]], 2),
                         " p =", signif(summary(fit110)$coef[2,4],
                                        2)),
           size = 2, fontface = 2) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# PCB 187 vs distance to IHSC (The Fork)
a.pcb187 <- log10(s.3$PCB187)
fit187 <- lm(a.pcb187 ~ d)
p187.dis <- ggplot(s.3, aes(y = PCB187, x = d)) +
  geom_point(shape = 21, colour = "black", fill = "white",
             size = 1.7, stroke = 0.8) +
  stat_smooth(method = "lm", col = "black", se = FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(aspect.ratio = 4/6) +
  xlab(expression(bold("distance (m)"))) + 
  ylab(expression(bold("PCB 187 (ng/g DW)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)) +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  annotate("text", x = 2500, y = 18,
           label = paste("R2 = ", signif(summary(fit187)$adj.r.squared,
                                         2),
                         " m =", signif(fit187$coef[[2]], 2),
                         " p =", signif(summary(fit187)$coef[2,4],
                                        2)),
           size = 2, fontface = 2) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# PCB 206 vs distance to IHSC (The Fork)
a.pcb206 <- log10(s.3$PCB206)
fit206 <- lm(a.pcb206 ~ d)
p206.dis <- ggplot(s.3, aes(y = PCB206, x = d)) +
  geom_point(shape = 21, colour = "black", fill = "white",
             size = 1.7, stroke = 0.8) +
  stat_smooth(method = "lm", col = "black", se = FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(aspect.ratio = 4/6) +
  xlab(expression(bold("distance (m)"))) + 
  ylab(expression(bold("PCB 206 (ng/g DW)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)) +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  annotate("text", x = 2500, y = 14,
           label = paste("R2 = ", signif(summary(fit206)$adj.r.squared,
                                         2),
                         " m =", signif(fit206$coef[[2]], 2),
                         " p =", signif(summary(fit206)$coef[2,4],
                                        2)),
           size = 2, fontface = 2) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# PCB 209 vs distance to IHSC (The Fork)
a.pcb209 <- log10(s.3$PCB209)
fit209 <- lm(a.pcb209 ~ d)
p209.dis <- ggplot(s.3, aes(y = PCB209, x = d)) +
  geom_point(shape = 21, colour = "black", fill = "white",
             size = 1.7, stroke = 0.8) +
  stat_smooth(method = "lm", col = "black", se = FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(aspect.ratio = 4/6) +
  xlab(expression(bold("distance (m)"))) + 
  ylab(expression(bold("PCB 209 (ng/g DW)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)) +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  annotate("text", x = 2500, y = 100,
           label = paste("R2 = ", signif(summary(fit209)$adj.r.squared,
                                         2),
                         " m =", signif(fit209$coef[[2]], 2),
                         " p =", signif(summary(fit209)$coef[2,4],
                                        2)),
           size = 2, fontface = 2) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# Total PCB vs distance to IHSC (The Fork)
a <- log10(rowSums(s.3))
fittPCB <- lm(a ~ d)
ptPCB.dis <- ggplot(s.3, aes(x = d, y = rowSums(s.3))) +
  geom_point(shape = 21, colour = "black", fill = "white",
             size = 1.7, stroke = 0.8) +
  stat_smooth(method = "lm", col = "black", se = FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(aspect.ratio = 4/6) +
  xlab(expression(bold("distance (m)"))) + 
  ylab(expression(bold(Sigma*"PCB (ng/g DW)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 7),
        axis.title.y = element_text(face = "bold", size = 7)) +
  theme(axis.text.x = element_text(face = "bold", size = 7),
        axis.title.x = element_text(face = "bold", size = 7)) +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  annotate("text", x = 2300, y = 2600,
           label = paste("R2 = ", signif(summary(fittPCB)$adj.r.squared,
                                         2),
                         " m =", signif(fittPCB$coef[[2]], 2),
                         " p =", signif(summary(fittPCB)$coef[2,4],
                                        2)),
           size = 2, fontface = 2) +
  annotation_logticks(sides = "l",
                      short = unit(0.5, "mm"),
                      mid = unit(1.5, "mm"),
                      long = unit(2, "mm"))

# Plot all plots in one panel (Figure 6)
ggarrange(p56.dis, p110.dis, p187.dis, p206.dis, p209.dis,
          ptPCB.dis, ncol=2)

# PCB profile plot (Figure 7) ---------------------------------------------

# Create average PCB congener profiles

tmp <- rowSums(s.3, na.rm = TRUE)
prof <- sweep(s.3, 1, tmp, FUN = "/")
prof.ave <- data.frame(colMeans(prof))
colnames(prof.ave) <- c("mean")
prof.sd <- data.frame(apply(prof,2,sd))
colnames(prof.sd) <- c("sd")
congener <- row.names(prof.ave)
prof.ave <- cbind(congener, prof.ave$mean, prof.sd$sd)
colnames(prof.ave) <- c("congener", "mean", "sd")
prof.ave <- data.frame(prof.ave)
prof.ave$mean <- as.numeric(as.character(prof.ave$mean))
prof.ave$sd <- as.numeric(as.character(prof.ave$sd))
prof.ave$congener <- as.character(prof.ave$congener)
#Then turn it back into a factor with the levels in the correct order
prof.ave$congener <- factor(prof.ave$congener,
                            levels = unique(prof.ave$congener))

# PCB Profile plot (Figure 7)
# Need to check the labels
ggplot(prof.ave, aes(x = congener, y = mean)) +
  geom_bar(position = position_dodge(), stat = "identity",
           fill = "black") +
  geom_errorbar(aes(ymin = mean, ymax = (mean+sd)), width = 0.9,
                position = position_dodge(0.9)) +
  xlab("") +
  ylim(0, 0.12) +
  theme_bw() +
  theme(aspect.ratio = 3/12) +
  ylab(expression(bold("Mass fraction "*Sigma*"PCB"))) +
  theme(axis.text.y = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 13)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  annotate("text", x = 11, y = 0.019, label = "PCB 11", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 15, y = 0.04, label = "PCBs 18+30", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 17, y = 0.052, label = "PCBs 20+28", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 23, y = 0.044, label = "PCB 31", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 37, y = 0.096, label = "PCB 52", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 42.9, y = 0.09, label = "PCBs 61+70+74+76", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 61, y = 0.092, label = "PCBs 90+101+113", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 66, y = 0.098, label = "PCB 95", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 74, y = 0.088, label = "PCB 110", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 77, y = 0.063, label = "PCB 118", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 80.9, y = 0.09, label = "PCBs 129+138+163", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 93, y = 0.069, label = "PCBs 147+149", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 95.2, y = 0.072, label = "PCBs 153+168", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 108, y = 0.074, label = "PCBs 180+193", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 111, y = 0.049, label = "PCB 187", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 119, y = 0.055, label = "PCBs 198+199", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 125, y = 0.037, label = "PCB 206", size = 4,
           fontface = 1, angle = 90) +
  annotate("text", x = 127.8, y = 0.035, label = "PCB 209", size = 4,
           fontface = 1, angle = 90)

# PCA plot (Figures 8 and S4) -----------------------------------------------------

# Prepare data from above section
t.prof <- data.frame(t(prof))
colnames(t.prof) <- s$Sample.label

# Create matrix to include other PCB profiles
s.prof <- matrix(NA, nrow = 128, ncol = 15)

# Add profile and PCBs names
colnames(s.prof) <- c("1016.Koh", "1242.Koh", "1232.Frame", "1221.Frame",
                      "1221.Koh", "1248.Koh", "1254.Koh", "1254.Frame",
                      "1260.Frame", "1262.Frame", "CR", "EPA", "emission",
                      "ECair", "LM")

rownames(s.prof) <- rownames(t.prof)

# Aroclors references
# 1016, 1242, 1248, 1254 and 1221:
# Reference
# Koh, W. X., K. C. Hornbuckle and P. S. Thorne (2015).
# "Human Serum from Urban and Rural Adolescents and Their
# Mothers Shows Exposure to Polychlorinated Biphenyls Not
# Found in Commercial Mixtures." Environmental Science &
# Technology 49(13): 8105-8112.
# https://pubs.acs.org/doi/10.1021/acs.est.5b01854

# 1221*, 1232*, 1254*, 1260*, 1262* 
# Frame, G. M., J. W. Cochran and S. S. Bowadt (1996).
# "Complete PCB congener distributions for 17 aroclor
# mixtures determined by 3 HRGC systems optimized for
#comprehensive, quantitative, congener-specific analysis."
# Hrc-Journal of High Resolution Chromatography 19(12): 657-668.
# http://dx.doi.org/10.1002/jhrc.1240191202

# Aroclor 1016
s.prof[1:128, 1] <- c(0.005862388, 0.000308547, 0.001748432,
                      0.036614214, 0.001954129, 0.01707292,
                      0.003394014, 0.083513319, 0.005965237,
                      0.001748432, 0.000205698, 0.022421063,
                      0.039391134, 0.038259796, 0.088038671,
                      0.009256402, 0.092666872, 0.064589119,
                      0.035791422, 0.007302273, 0.016867222,
                      0.005450992, 0.09102129, 0.00092564,
                      0.018204258, 0.023346704, 0.008433611,
                      0.015221639, 0.002571223, 0.044739278,
                      0.025815078, 0.004422503, 0.015633035,
                      0.030443279, 0.009770647, 0.002674072,
                      0.046590558, 0.000205698, 0.001439885,
                      0.005348144, 0.000719942, 0.026226473,
                      0.001131338, 0.022523912, 0.010182043,
                      0.001234187, 0.000102849, 0.000205698,
                      0, 0, 0, 0, 0, 0, 0.001234187, 0, 0.000102849,
                      0.000102849, 0, 0.000205698, 0.000822791,
                      0.001542734, 0.000308547, 0.000205698,
                      0.000102849, 0.006376633, 0.000308547,
                      0.000308547, 0.000719942, 0.000102849,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0)

# Arolcor 1242
s.prof[1:128, 2] <- c(0.00595115945003078,	0.000410424789657295,
                      0.0021547301457008,	0.032628770777755,
                      0.0017443053560435,	0.0152883234147342,
                      0.00297557972501539,	0.0726451877693412,
                      0.00523291606813051,	0.00153909296121486,
                      0.000102606197414324,	0.0187769341268213,
                      0.0341678637389698,	0.0335522265544839,
                      0.0769546480607428,	0.00831110199056023,
                      0.0770572542581572,	0.0533552226554484,
                      0.0298584034475682,	0.00615637184485943,
                      0.0140570490457624,	0.00471988508105889,
                      0.0758259798891853,	0.00082084957931459,
                      0.0179560845475067,	0.0204186332854504,
                      0.00687461522675969,	0.0130309870716191,
                      0.0021547301457008,	0.0379642930432998,
                      0.021855120049251,	0.00359121690950133,
                      0.0129283808742048,	0.0255489431561666,
                      0.0082084957931459,	0.00205212394828648,
                      0.038990355017443,	0.00102606197414324,
                      0.0167248101785348,	0.00451467268623025,
                      0.0120049250974759,	0.0580751077365073,
                      0.00153909296121486,	0.0193925713113072,
                      0.0332444079622409,	0.00133388056638621,
                      0.000102606197414324,	0.000205212394828648,
                      0,	0.00266776113277242,	0.000205212394828648,
                      0.000102606197414324,	0.00225733634311512,
                      0.000615637184485943,	0.00338600451467268,
                      0.00235994254052945,	0.00307818592242971,
                      0.00471988508105889,	0,	0.000513030987071619,
                      0.00707982762158834,	0.00205212394828648,
                      0.00123127436897189,	0.000205212394828648,
                      0.000102606197414324,	0.00677200902934537,
                      0.000307818592242971,	0.00471988508105889,
                      0.000718243381900266,	0.000102606197414324,
                      0.00400164169915863,	0.000513030987071619,
                      0.000307818592242971,	0.00790067720090293,
                      0.000307818592242971,	0.000205212394828648,
                      0.00636158423968807,	0.000205212394828648,
                      0.000205212394828648,	0.000718243381900266,
                      0.000102606197414324,	0,	0.000307818592242971,	0,
                      0.000102606197414324,	0.000205212394828648,
                      0.000102606197414324,	0.000102606197414324,	0,
                      0.000102606197414324,	0,	0.000102606197414324,
                      0.000513030987071619,	0,	0.000307818592242971,	0,
                      0.000102606197414324,	0,	0,	0,	0,
                      0,	0,	0,          0,	0,	0,
                      0,	0,	0,	0,	0,	0,
                      0,	0,	0,	0,          0,	0,
                      0,	0,	0,	0,	0,	0,
                      0,	0,	0)

# Aroclor 1232*
s.prof[1:128, 3] <- c(0.176721684689812,	0.0223107569721116,
                      0.117017643710871,	0.0608992601024472,
                      0.0144564598747866,	0.121969265793967,
                      0.0125782583949915,	0.00671599317017644,
                      0.0342629482071713,	0.00563460443938532,
                      0.0122367672168469,	0.0365964712578258,
                      0.0325554923164485,	0.0553215708594195,
                      0.020774046670461,	0.020375640295959,
                      0.0471257825839499,	0.0122367672168469,
                      0.000113830392714855,	0.00904951622083096,
                      0.00421172453044963,	0.00529311326124075,
                      0.00136596471257826,	0,	0,	0.0105293113261241,
                      0.000569151963574274,	0,	0.00700056915196357,
                      0.00142287990893569,	0.0076835515082527,
                      0.0261809903244166,	0.00421172453044963,
                      0.0155378486055777,	0.00404097894137735,
                      0.00119521912350598,	0.00523619806488332,	0,	0,
                      0.0210017074558907,	0,	0.0196357427433125,
                      0.00990324416619237,	5.69151963574274E-05,
                      0.000967558338076267,	0.00113830392714855,
                      0.00221969265793967,	0,	0.00216277746158224,
                      0,	0,	0,	0.00432555492316449,	0,
                      0.000569151963574274,	0,	0,	0,	0.00369948776323278,
                      0,	0.000113830392714855,	0.00113830392714855,
                      0.00341491178144565,	0.00250426863972681,
                      0.00193511667615253,	0.00244735344336938,
                      0.000512236767216847,	0.00239043824701195,	0,
                      0.00216277746158224,	0,	0,	0.000170745589072282,
                      0,	0,	0.00244735344336938,	0,	0,
                      5.69151963574274E-05,	0,	0,	0,	0,	0,
                      0,	0,	0,
                      0,	0,	0.000569151963574274,	0,	0,	0,
                      0.000853727945361412,	0,	5.69151963574274E-05,	0,          0,
                      0,	0,	0,	0,	0,	0,	0,
                      0,	0,	0,          0.000113830392714855,	0,	0,
                      0,	0,	0,	0,	0,          0,	0,
                      0,	0,	0,	0,	0,	0,	0,
                      0,	0,	0)

# Aroclor 1221*
s.prof[1:128, 4] <- c(0.365,	0.0389,	0.208,	0.0631,	0.0177,	0.126,
                      0.0173,	0.00816,	0.039,	0.00755,	0.0174,	0.0426,
                      0.00489,	0.00795,	0.00347,	0.00316,	0.00612,	0.00173,
                      0,	0.00143,	0.000918,	0.000816,	0.00051,	0,
                      0,	0.00122,	0,	0,	0.000612,	0,
                      0.000918,	0.00265,	0.000408,	0.00153,	0.000306,	0,
                      0.000408,	0,	0,	0.00224,	0,	0.00214,
                      0.00102,	0,	0,	0,	0.000102,	0,
                      0.000204,	0,	0,	0,	0.00051,	0,
                      0.000204,	0,	0,	0,	0.000714,	0,
                      0,	0,	0.00051,	0.000408,	0.000306,	0.000306,
                      0,	0.000408,	0,	0.000204,	0,	0,
                      0,	0,	0,	0.00051,	0,	0,
                      0,	0,	0,	0,	0,	0,
                      0,	0,	0,	0,	0,	0,
                      0,	0,	0,	0,	0,	0,
                      0,	0,	0,	0,	0,	0,
                      0,	0,	0,	0,	0,	0,
                      0,	0,	0,	0,	0,	0,
                      0,	0,	0,	0,	0,	0,
                      0,	0,	0,	0,	0,	0,
                      0,	0)

# Aroclor 1221
s.prof[1:128, 5] <- c(0.38319276496291,	0.0363784168275582,
                      0.191443958947262,	0.0609694136774718,
                      0.00772279239914643,	0.0338380245909969,
                      0.0156488161772178,	0.105883548419876,
                      0.0159536632456051,	0.00589370998882228,
                      0.000711309826237171,	0.0249974596077634,
                      0.00538563154151001,	0.00558886292043492,
                      0.0115841885987196,	0.00132100396301189,
                      0.0115841885987196,	0.0082308708464587,
                      0.00447109033634793,	0.00111777258408698,
                      0.00233716085763642,	0.000812925515699624,
                      0.0108728787724825,	0.000203231378924906,
                      0.00254039223656133,	0.00243877654709887,
                      0.000914541205162077,	0.00162585103139925,
                      0.000304847068387359,	0.00477593740473529,
                      0.00264200792602378,	0.000508078447312265,
                      0.00162585103139925,	0.00315008637333604,
                      0.00101615689462453,	0.000304847068387359,
                      0.00508078447312265,	0.000203231378924906,
                      0.00193069809978661,	0.000609694136774718,
                      0.00142261965247434,	0.00619855705720963,
                      0.000203231378924906,	0.00233716085763642,
                      0.00375978051011076,	0.000203231378924906,
                      0.000101615689462453,	0.000101615689462453,
                      0,	0.000304847068387359,	0,	0,
                      0.000304847068387359,	0.000101615689462453,
                      0.000508078447312265,	0.000304847068387359,
                      0.000508078447312265,	0.000609694136774718,
                      0.000101615689462453,	0.000101615689462453,
                      0.00111777258408698,	0.000203231378924906,
                      0.000203231378924906,	0.000101615689462453,
                      0,	0.00101615689462453,	0.000101615689462453,
                      0.000609694136774718,	0.000101615689462453,	0,
                      0.000508078447312265,	0.000101615689462453,
                      0.000101615689462453,	0.00111777258408698,
                      0.000101615689462453,	0.000101615689462453,
                      0.000812925515699624,	0,	0,
                      0.000203231378924906,	0,	0,
                      0.000101615689462453,	0,	0,
                      0.000101615689462453,	0,	0,	0,
                      0,	0,	0,          0.000203231378924906,
                      0,	0.000101615689462453,   0,	0,
                      0,	0,	0,	0,	0,
                      0,	0,	0,	0,	0,
                      0,	0,	0,	0,	0,
                      0,	0,	0,	0,	0,
                      0,	0,          0,	0,	0,
                      0,	0,	0,	0,	0,	0)

# Aroclor 1248
s.prof[1:128, 6] <- c(0.000506534292371593,	0,	0.000202613716948637,
                      0.00314051261270388,	0.000101306858474319,
                      0.00131698916016614,	0.000202613716948637,
                      0.00830716239489413,	0.000405227433897275,
                      0.000101306858474319,	0,	0.00202613716948637,
                      0.0115489818660723,	0.0110424475737007,
                      0.037179617060075,	0.00222875088643501,
                      0.0384966062202411,	0.0239084185999392,
                      0.0141829601864046,	0.00121568230169182,
                      0.00445750177287002,	0.00131698916016614,
                      0.0496403606524161,	0.000101306858474319,
                      0.00699017323472799,	0.0317090467024617,
                      0.00780062810252254,	0.018133927666903,
                      0.00263397832033228,	0.0635194002633978,
                      0.0261371694863742,	0.00455880863134434,
                      0.0167156316482626,	0.0402188228143045,
                      0.011447675007598,	0.00243136460338365,
                      0.0746631546955729,	0.000607841150845912,
                      0.0289737615236551,	0.00466011548981866,
                      0.018741768817749,	0.105966973964137,
                      0.00233005774490933,	0.032418194711782,
                      0.0576436024718873,	0.00111437544321751,
                      0.000101306858474319,	0.000202613716948637,
                      0,	0.00384966062202411,	0.000202613716948637,
                      0.000101306858474319,	0.00729409381015094,
                      0.00182352345253774,	0.0113463681491237,
                      0.00881369668726572,	0.0103332995643805,
                      0.0151960287711478,	0,	0.00131698916016614,
                      0.0238071117414649,	0.00628102522540776,
                      0.00405227433897275,	0.000506534292371593,
                      0.000303920575422956,	0.0218822814304528,
                      0.000709148009320231,	0.0157025630635194,
                      0.00212744402796069,	0.000202613716948637,
                      0.0160064836389423,	0.0016209097355891,
                      0.000911761726268868,	0.029885523249924,
                      0.00111437544321751,	0,	0.0249214871846824,
                      0.000607841150845912,	0.000709148009320231,
                      0.00455880863134434,	0.000303920575422956,
                      0.000101306858474319,	0.00182352345253774,
                      0,	0.000303920575422956,	0.000911761726268868,
                      0.000506534292371593,	0.000405227433897275,
                      0.000101306858474319,	0.000607841150845912,
                      0.000202613716948637,	0.000405227433897275,
                      0.00253267146185797,	0,	0.00222875088643501,
                      0,	0.000506534292371593,	0.000202613716948637,
                      0.000202613716948637,	0.000101306858474319,	0,
                      0.000101306858474319,	0,	0,	0.000101306858474319,
                      0,	0,	0.000303920575422956,	0.000101306858474319,
                      0,	0.000202613716948637,	0,	0,
                      0,	0,	0,	0,          0,
                      0.000101306858474319,	0,	0,	0,
                      0,	0,          0.000101306858474319,	0,
                      0,	0)

# Aroclor 1254
s.prof[1:128, 7] <- c(0.00031126789790413,	0,	0.000207511931936086,
                      0.000103755965968043,	0,	0.000103755965968043,
                      0,	0.00031126789790413,	0,	0,	0,
                      0.000103755965968043,	0.000207511931936086,
                      0.000103755965968043,	0.000622535795808259,	0,
                      0.000622535795808259,	0.00031126789790413,
                      0.000207511931936086,	0,	0.000103755965968043,
                      0,	0.0014525835235526,	0.000103755965968043,
                      0.000207511931936086,	0.00269765511516912,
                      0.000207511931936086,	0.00124507159161652,
                      0.000207511931936086,	0.00871550114131563,
                      0.000933803693712389,	0.000207511931936086,
                      0.000830047727744345,	0.00477277443452999,
                      0.000830047727744345,	0.000103755965968043,
                      0.0162896866569828,	0.000207511931936086,
                      0.01431832330359,	0.00031126789790413,
                      0.00747042954969911,	0.0875700352770284,
                      0.000830047727744345,	0.00425399460468977,
                      0.0314380576883171,	0.000207511931936086,	0,
                      0,	0,	0.0014525835235526,	0.000726291761776302,
                      0,	0.01431832330359,	0.00249014318323304,
                      0.0157709068271426,	0.0221000207511932,
                      0.0229300684789375,	0.0368333679186553,	0,
                      0.000622535795808259,	0.0619423116829218,
                      0.00539531023033824,	0.00715916165179498,
                      0.000207511931936086,	0.000103755965968043,
                      0.0232413363768417,	0.000207511931936086,
                      0.0466901846856194,	0.000933803693712389,
                      0.000103755965968043,	0.0796845818634572,
                      0.00830047727744345,	0.00498028636646607,
                      0.0873625233450924,	0.00446150653662586,	0,
                      0.144324548661548,	0.00217887528532891,
                      0.00280141108113717,	0.0723179082797261,
                      0.00373521477484955,	0.000933803693712389,
                      0.016912222452791,	0.000415023863872173,
                      0.00217887528532891,	0.0052915542643702,
                      0.00259389914920108,	0.00612160199211455,
                      0.00124507159161652,	0.00684789375389085,
                      0.00124507159161652,	0.0045652625025939,
                      0.0209587051255447,	0,	0.0378709275783358,
                      0.000207511931936086,	0.00850798920937954,
                      0.00280141108113717,	0.00466901846856194,
                      0.00124507159161652,	0.000415023863872173,
                      0.00166009545548869,	0.000103755965968043,
                      0.000207511931936086,	0.00103755965968043,
                      0.000207511931936086,	0.00031126789790413,
                      0.00498028636646607,	0.00124507159161652,
                      0.000103755965968043,	0.00103755965968043,
                      0.00031126789790413,	0.000726291761776302,
                      0.000207511931936086,	0.000207511931936086,
                      0.000103755965968043,	0.000103755965968043,	0,
                      0.000103755965968043,	0,	0,	0,
                      0.000103755965968043,	0,	0.000103755965968043,
                      0,	0,	0)

# Aroclor 1254*                  
s.prof[1:128, 8] <- c(0,	0,	0,	0.000594369495614195,	0,	0.000208027604820148,
                      0, 0.00136809534786731, 0, 0, 0, 0.000287941677446609,
                      0.000911940602049987, 0.000872593197713001,
                      0.00260191465781369, 0, 0.00198918354245831,
                      0.00166862915030772, 0.000392709496013762, 0,
                      0.000354387775151916, 0, 0.00289613393320768, 0,
                      0.000717473099934463, 0.00280218082867534,
                      0.000138044759409004,	0.00151912144856345, 0,
                      0.0252219166794649, 0.000503623968330042, 0,
                      0.00123083316446924, 0.0112542011649302,
                      0.00120142796937609, 0, 0.055246137716281, 0,
                      0.00563601169278585, 0.000238776347213613,
                      0.00184436561170285, 0.0446233530545292,
                      0.000205402680853391, 0.00603125372209212,
                      0.0103412243675552, 0, 0, 0, 0, 0.000290445621195225,
                      0, 0, 0.0113595720524287, 0.00492966434048139,
                      0.023834190338366, 0.0131132075648303, 0.032154025385486,
                      0.0412355766380525, 0, 0.000924312063840261,
                      0.0824108341127801, 0.00952750988406525,
                      0.0132823498927887, 0, 0.000246418581422322,
                      0.06422757904171, 0.000457144407564762,
                      0.031050918744746, 0.00156625125735295,
                      0.000258161741838105, 0.030745347080332, 0,
                      0.00295206342791111, 0.0954006495783798,
                      0.0018738946099829, 0.00237594434774808,
                      0.0755192954357059, 0.000992046427009297,
                      0.00154052010640044, 0.0740419234600759,
                      0.00614708382596712, 0.00194463369165087,
                      0.0235083368236707, 0.00111427851233479,
                      0.00383721481635986, 0.0132945107622253,
                      0.00718289050276573, 0.00435156477580248,
                      0.00152408245038403, 0.0101093765937139,
                      0.00244163405989748, 0.0068920803064644,
                      0.0385830901970306, 0, 0.038728675474907,
                      0.000382943954148511, 0.00834354490684839,
                      0.00415066491819365, 0.00534023264978126,
                      0.0013963620275577, 0.000712503113277994,
                      0.00349992097724663, 0, 0.000432299471893389,
                      0.00207927195179908, 0.000345100279361987,
                      0.00105576995124697, 0.00720826691712498,
                      0.00180870636882629, 0, 0.0025925505781069,
                      0.000109005844561664, 0.000767710670561565, 0,
                      0.000141070657354089, 0, 0, 0, 0.000140371398223474,
                      0, 0, 0, 0.00017327594035687, 0, 0.000316450629917424,
                      0, 0.000133403101312959, 0)

# Aroclor 1260*
s.prof[1:128, 9] <- c(0.00023523316498943,	0,	0,	0.000191582782409939,	0,	
                      7.40992179141525E-05,	0,	0.000429500611295873,	0,	0,	
                      0,	7.61257063477238E-05,	0.000142851656494913,	0.000154332837109262,	
                      0.000470543307929344,	0,	0.000339866910875251,	0.000301291981962105,	
                      0.000147205397190366,	0,	0,	0,	0.000377187974937975,	0,	
                      7.11275621173537E-05,	0.000051167602277935,	0,	5.46778873259085E-05,	0,	
                      0.000348939967936172,	0,	0,	0,	0.000139952587054788,	0,	
                      0,	0.00244183048162357,	0,	0.000181784457337379,	0,	
                      0.000408349559802985,	0.000921043379014284,	0,	0.000106796994620973,	
                      0.000233112255603598,	0,	0,	0,	0,	0,	0,	
                      0,	0,	5.14662569146743E-05,	0.00110128699121323,	0.000102225190982825,	
                      0.000986165401492073,	0.00412767324018134,	0,	0,	0.0316588019904055,	
                      7.40427546675558E-05,	0.00303521056551929,	0,	0,	0.0247888125752524,	
                      0,	0.000438278256594998,	0,	0,	0.00222431239894061,	0,	
                      0.000070901851233696,	0.0134219608188498,	0,	0,	0.00489253045797529,	
                      0,	0,	0.091901177777766,	0.00221509681016245,	0.000656636887151579,	
                      0.0293297809594372,	0.000656636887151579,	0.00346588037065508,	0.0416394067968281,	
                      0.0147162944740267,	0.000234361001531615,	0,	0.0264852939214048,	
                      0.00617176965102022,	0.0115887044760929,	0.088420740460793,	0,	
                      0.09482173387334,	0,	0.00588683387753021,	0.00701502950494056,	
                      0.0415547793537635,	0.0122237128996653,	0.00703413659915134,	0.0500820437271584,	
                      0.00176050181851026,	0.00594372472548385,	0.0259222190312305,	0.00836445396189673,	
                      0.020489341204431,	0.120405160961562,	0.0243890126592658,	0.00553933847331223,	
                      0.0545213470857056,	0.00103104070694597,	0.00831909253365024,	0.00170028036851201,	
                      0.0208971876525863,	0.00846531348503495,	0.0109844843853129,	0.000708704422257327,	
                      0.0189602766772492,	0.00254189519341989,	0.0024546875728075,	0.00338005114148336,	
                      0.0141184921393349,	0.00102154556386432,	0.00533404655680236,	0.000499650885165893,	
                      0.00127582745020947,	0)
                      
# Aroclor 1262*
s.prof[1:128, 10] <- c(0.000349088878028346,	0,	6.98177756056692E-05,	0.000767995531662361,	0,	
                       0.00160580883893039,	0,	0,	0.000349088878028346,	0,	0,	
                       0.000349088878028346,	0.00139635551211338,	0.00202471549256441,	0.000698177756056692,	
                       0.000698177756056692,	0.00167562661453606,	0.000488724429239685,	0,	
                       0.000279271102422677,	0,	0.000139635551211338,	0,	0,	0,	
                       0.000418906653634015,	0,	0,	6.98177756056692E-05,	0,	
                       0.000279271102422677,	0.00111708440969071,	0,	0.000767995531662361,	
                       0,	0,	0,	0,	0,	0.00195489771695874,	0,	
                       0.0009076310828737,	0.000418906653634015,	0,	0,	0,	0,	
                       0,	0,	0,	0,	0,	0.0054457864972422,	0,	
                       0.00111708440969071,	0,	0,	0,	0.0157788172868812,	0,	
                       0,	6.98177756056692E-05,	0.0129861062626545,	0.00153599106332472,	
                       0.000279271102422677,	0.000628359980451023,	0,	0.000628359980451023,	0,	
                       0.000558542204845354,	0,	0,	0,	0,	0,	
                       0.00125671996090205,	0,	0,	0,	0.0027228932486211,	0,	
                       0,	0,	0.00823849752146897,	0,	0.00174544439014173,	0,	
                       0,	0,	0.0991412413600503,	0.0140333728967395,	0,	
                       6.98177756056692E-05,	0.0599734692452699,	0,	0.0507575228653215,	0,	
                       0.00223416881938142,	0.0105424841164561,	0.0121482929553864,	0,	0,	
                       0.0387488654611464,	0.0883893039167772,	0.0401452209732598,	0.00970467080918801,	
                       0.0127068351602318,	0,	0.206451162465964,	0,	0.0168260839209663,	
                       0,	0.0018152621657474,	0.000488724429239685,	0,	0.0592054737136075,	
                       0.00865740417510298,	0.00188507994135307,	0.00900649305313133,	0.0693988689520352,	
                       0.0198980660476157,	0.0566222160161977,	0.0150806395308246,	0.00237380437059275,	
                       0.00383997765831181,	0.00244362214619842,	0.0175940794526286,	0)

# Cedar Rapids soil
# Reference:
# Martinez, A., N. R. Erdman, Z. L. Rodenburg,
# P. M. Eastling and K. C. Hornbuckle (2012).
# "Spatial distribution of chlordanes and PCB congeners
# in soil in Cedar Rapids, Iowa, USA." Environmental
# Pollution 161: 222-228.
# https://doi.org/10.1016/j.envpol.2011.10.028
s.prof[1:128, 11] <- c(0,	0,	0.000695726446668982,	0.000543899755431695,	0,	
                       3.80925064577676E-05,	5.58000754472744E-06,	0.00134346806868051,	9.71752629817895E-06,	
                       0,	0.00019702894297409,	0.000734166346003567,	0.000547369910401551,	
                       0.000727112127538938,	0.00388753700900152,	0,	0.0097654551842655,	
                       0.00541464162202491,	0.00202108152668011,	0.000737981481547752,	0.00162501145262242,	
                       0.000049765227736859,	0.00692661515030913,	5.41466011686214E-05,	0.00356347470214116,	
                       0.00146177170217015,	0.000235526832022201,	0.00160079509497262,	0.00529254974806944,	
                       0.00396031785147721,	0.000964248566598812,	0.000054434488078561,	0.000461458793821235,	
                       0.00373814046348401,	0.000103868612892997,	0.000341079415589403,	0.0192311492123222,	
                       0,	0.00508236196176955,	0.000160094067277145,	0.000501853942627791,	
                       0.0147793505839514,	0.00206072647413032,	0.00279567636625839,	0.00424248329226847,	
                       4.22455670371034E-05,	0,	0.000108032871789382,	0.000469234924720485,	
                       0.000222569076578602,	0,	0,	0.00101248115834472,	0.00149321729853721,	
                       0.00554740506034819,	0.00715664143114321,	0.00497336493419445,	0.00783046622878104,	
                       0.000588893453038876,	0,	0.0582522070681208,	0.00166133164996668,	
                       0.00448923415948433,	1.54923174870925E-05,	0,	0.0411639703670326,	0,	
                       0.0212485213533608,	1.82020003935801E-05,	0,	0.0362364015496706,	
                       0.000757982264014871,	0.000534705692268866,	0.0938207244349853,	6.81232550475565E-05,	
                       0.00454358979143425,	0.081111338310827,	0,	0.000206356577992057,	
                       0.150895526488588,	0.000461387498788824,	0,	0.0167265317977068,	0,	
                       0.000503141432375496,	0.00575810394369525,	0.000823235911669785,	0.00453207508500385,	
                       2.80094698643753E-05,	0.00432092071935958,	0.000105910631956616,	0.00595254334808654,	
                       0.0581271676348037,	0,	0.0977937778166942,	0.00231048905914966,	
                       0.00308788477993534,	0,	0.0186435673957838,	0.00115981588877121,	
                       0.0010381578238995,	0,	6.25643201808999E-06,	0.000149716896940696,	
                       0.00774751020749529,	0.002276875691155,	0.004170566749132,	0.0484050656650062,	
                       0.00835160377685219,	0.00026911898394651,	0.0274025627604381,	0,	
                       0.0013330710188098,	5.12468185989664E-06,	0.00207922541723306,	0.000280362927298974,	
                       0.000199109827636975,	1.14853262370639E-05,	0.00404670670071086,	8.24767252368197E-05,	
                       7.59586203052473E-05,	0.000313154594167487,	0.00111513166119744,	0,	0.013851659057069,	
                       0.000069277371266649,	0.00393695488982279,	0.0220242854301517)

# EPA background soil
# Reference:
# USEPA (2007). Pilot Survey of Levels of Polychlorinated
# Dibenzo-P-Dioxins (PCDDs), Polychlorinated Dibenzofurans
# (PCDFs), Polychlorinated Biphenyls (PCB) and Mercury in
# Rural Soils of the U.S. Washington, DC.
s.prof[1:128, 12] <- c(0.00346023260988123,	0.00164062753054713,	0.00346023260988123,	0.00584659992704071,	
                       0,	0.00325142546962978,	0,	0.00885938866495453,	0,	0,	
                       0.00483239381724793,	0.00435512035381604,	0.00328156261621622,	0.00441477953674501,	
                       0.00876989989056105,	0.00187926426226308,	0.0103623680966912,	0.00574561129413521,	
                       0.00432529076235155,	0.00101420610979277,	0.00232670813423048,	0.000298349813947501,	
                       0.00975427640888933,	0,	0.0045042683111385,	0.00439303567392396,	
                       0.000707824466504409,	0.00193892344519207,	0.00439647141352963,	0.0129758722870546,	
                       0.00234440645633471,	0,	0.00175994589640511,	0.0076960345978393,	
                       0.00184943467079859,	0.00111753730443138,	0.0159751351968186,	0,	
                       0.00387784689038414,	0.000894887743934802,	0.00250568568301745,	0.0210000323910034,	
                       0,	0.00444460912820951,	0.00942615090277991,	0,	0,	0,	
                       0.000389789048939191,	0.00372869893306167,	0,	0,	0.00295312955498485,	
                       0.0014425361523963,	0.00677131726244,	0.0053095181327225,	0.00771670545897373,	
                       0.0121498024563788,	0.00103046203901015,	0,	0.0350795995622442,	
                       0.00290704403430299,	0.0067414876709755,	0,	0,	0.0241321394947752,	
                       0,	0.0205273273134559,	0,	0,	0.0187329834397018,	
                       0.00340653692106054,	0.00202841221958555,	0.0323888212641896,	0,	
                       0.00337089298344507,	0.036153464854966,	0,	0.000927411280193042,	
                       0.0755285255880972,	0.00378835811599066,	0,	0.00874007029909656,	0,	
                       0.00116724488339322,	0.014624658498266,	0.00519034891482185,	0.00638353257340158,	
                       0,	0.00614489584168563,	0.00158096834761815,	0.0113949039394365,	
                       0.0278608384278368,	0,	0.0748126153929494,	0.00586827083784095,	
                       0.00536932646360881,	0,	0.0204034405617134,	0.00417614280502907,	
                       0.00441477953674501,	0,	0,	0.00125284284150872,	0.011186096799185,	
                       0.00745739786612334,	0.00659233971365304,	0.0423580198795806,	0.0178002549201667,	
                       0.000573588815525642,	0.0346321556902768,	0.00149147957322467,	0.00566762237825374,	
                       0,	0.0196577007751012,	0.00626421420754361,	0.0085610927503096,	
                       0.000306279152759819,	0.0301278873791383,	0.00219940653025763,	0.00214773058544353,	
                       0.00584659992704071,	0.0170326967262257,	0.00140199079883119,	0.0295312955498485,	
                       0.00519034891482185,	0.0120213253601909,	0.0356761913915341)

# IHSC gas emissions 
# Reference (Data):
# Martinez, A., Awad, A. M., Herkert, N. J., Hornbuckle, K. C. (2018).
# PCB congener data of gas-phase, freely-dissolved water, air-water
# fugacity ratios and air-water fluxes in Indiana Harbor and Ship Canal,
# IN, USA. PANGAEA, https://doi.org/10.1594/PANGAEA.894908,
# Reference (Paper):
# Martinez, A., Awad, A. M., Herkert, N. J., Hornbuckle, K. C. (2019)
# Determination of PCB fluxes from Indiana Harbor and Ship Canal using
# dual-deployed air and water passive samplers. Environmental Pollution,
# 244, 469-476, https://doi.org/10.1016/j.envpol.2018.10.048
s.prof[1:128, 13] <- c(0.000204, 0.000241, 0.000132, 0.0197, 0.0000802, 0.00103, 0.000148, 0.00198, 0.000214,
                       0.0013, 0.000402, 0.0134, 0.0227, 0.0399, 0.0985, 0.0161, 0.0595, 0.00649, 0.00924,
                       0.00723, 0.0113, 0.00777, 0.036, 0.00035, 0.00666, 0.0285, 0.00666, 0.0306, 0.00538,
                       0.105, 0.0389, 0.0133, 0.0168, 0.0623, 0.0286, 0.00813, 0.0993, 0.000676, 0.0105,
                       0.00668, 0.00512, 0.0306, 0.0013, 0.0319, 0.0179, 0.000766, 0.000121, 0.000199, 0.000204,
                       0.00116, 0.0000568, 0.0000571, 0.00275, 0.000947, 0.00834, 0.00261, 0.00436, 0.0049,
                       0.0000327, 0.000784, 0.00972, 0.00395, 0.00192, 0.000495, 0.000241, 0.0155, 0.00117, 0.00503,
                       0.00118, 0.000151, 0.00183, 0.00014, 0.000293, 0.00937, 0.000132, 0.000287, 0.00326, 0.0000961,
                       0.0001, 0.0014, 0.0000828, 0.0000286, 0.000689, 0.0000204, 0.000138, 0.000763, 0.000669,
                       0.0000532, 0.0000309, 0.000275, 0.0000976, 0.000186, 0.00167, 0.00000684, 0.00109, 0.0000186,
                       0.000117, 0.0000831, 0.000136, 0.0000591, 0.0000279, 0.000222, 0.0000101, 0.0000372, 0.000125,
                       0.0000455, 0.000142, 0.000331, 0.000119, 0.0000203, 0.000258, 0.00000699, 0.0000308, 0.00000769,
                       0.0000385, 0.0000197, 0.0000265, 0.0000056, 0.0000228, 0.0000113, 0.00000792, 0.0000171,
                       0.0000327, 0.00000282, 0.000014, 0.00000315, 0.00000536, 0.00000733)

# East Chicago air PCB  
# Reference:
# Marek, R. F., P. S. Thome, N. J. Herkert, A. M. Awad and
# K. C. Hornbuckle (2017). "Airborne PCBs and OH-PCBs
# Inside and Outside Urban and Rural US Schools."
# Environmental Science & Technology 51(14): 7853-7860.
# https://pubs.acs.org/doi/10.1021/acs.est.7b01910
s.prof[1:128, 14] <- c(0.0623551465438699,	0.00864998886996662,	0.0368693088893442,	0.0428708083072982,	
                       0.000242567288084444,	0.0118102160592177,	0.00328222881997639,	0.0475783478593331,	
                       0.00436255263401141,	0.00179077585686455,	0.023440617481507,	0.0108191436856317,	
                       0.0187002886293166,	0.0185287263601378,	0.0450427957464993,	0.00823074551710456,	
                       0.0319285038115658,	0.0188287649627368,	0.0115134465473421,	0.00155717692885623,	
                       0.00625024093713563,	0.00275600770102887,	0.031158751717337,	0.000347762621538083,	
                       0.00351474669983469,	0.0110148734348883,	0.00275470294838588,	0.00736710796324867,	
                       0.00150385658540303,	0.0325083600354721,	0.00713239833404103,	0.00249983038581995,	
                       0.00671074923978415,	0.021090777838774,	0.00725958881710749,	0.00182561203203396,	
                       0.0625200301468843,	0.000274111966440388,	0.00523637232715869,	0.00249050307686171,	
                       0.00296759528635737,	0.0328771851672126,	0.000606927381653134,	0.0126475688568827,	
                       0.011587841419723,	0.000461159340326436,	2.54172004686119E-05,	2.32775118856395E-05,	
                       0.000209155999432289,	0.000443080934013614,	6.13979906441399E-05,	0,	
                       0.00213574128480548,	0.00118041332269262,	0.0093945020337131,	0.00317833846978302,	
                       0.00755761794158186,	0.0119036531519991,	9.20689634081032E-06,	0.000250559891982451,	
                       0.0326860936510236,	0.00435440545820426,	0.00566185353829537,	0.000170152178117742,	
                       0.000111124009515486,	0.0324073949144107,	0.000283460314373287,	0.011782635987312,	
                       0.000968443828234564,	0.000190120053249235,	0.00538952141208785,	0.00106467068624822,	
                       0.000617177710090503,	0.0296908095382985,	0.000139862998453105,	0.000451173219549266,	
                       0.0183801536627304,	4.96713323207171E-05,	8.47268830640934E-05,	0.0232377350302627,	
                       0.00120779315635536,	0.000351880247574193,	0.00928766524005977,	0.000149859772002686,	
                       0.00180484583855641,	0.00842074573437242,	0.00429374960175442,	0.00115074993887306,	
                       0.000402831789182228,	0.00385058943684151,	0.00116927979229585,	0.00262855157692209,	
                       0.0211241177304563,	6.47728410410996E-06,	0.0178696300482147,	0.000136589299373177,	
                       0.00220100100300501,	0.00116279470223509,	0.00178155349606701,	0.000856003989013376,	
                       0.000255589669265129,	0.0023940241087499,	0.000067421639294609,	0.000428192579156261,	
                       0.00151364967702005,	0.000545437777732682,	0.00205509434413492,	0.00491304458805774,	
                       0.00211580910621333,	0.000398388238961998,	0.00528685024827212,	0,	
                       0.000350559294280105,	0.000079764566809351,	0.000934684515567536,	9.52799003879923E-05,	
                       0.000715780150371951,	8.98927245198691E-06,	0.00233373278434139,	0.00018026295185288,	
                       0.000478240990029441,	0.00107790598343666,	0.00162438628227095,	0,	0.00143953268003797,	
                       0.000222130856985418,	0.000601368297860797,	0.000166807724016178)

# Lake Michigan emission
# Reference:
# Boesen, A. C., A. Martinez and K. C. Hornbuckle (2020).
# "Air-water PCB fluxes from southwestern Lake Michigan
# revisited." Environmental Science and Pollution Research
# 27(9): 8826-8834.
# https://link.springer.com/article/10.1007/s11356-019-05159-1
s.prof[1:128, 15] <- c(0.00699336152201969,	0.00213197469418328,	0.00411529745234851,	0.0217750572433371,	
                       0.00103034494300888,	0.00865877768485397,	0.00356210938371506,	0.0406726245561402,	
                       0.00268989615522224,	0.00137321584396853,	0.0582857610192436,	0.0184317570775284,	
                       0.0234316732391625,	0.0232624592077576,	0.0447675692252999,	0.00834080240446592,	
                       0.0427211000587293,	0.0217607929860258,	0.0157804073263246,	0.00490100347200848,	
                       0.0099767394495305,	0.00538078778489643,	0.0408657863294315,	0.00184332438529258,	
                       0.0100969818068506,	0.0186164984133849,	0.00298729858261293,	0.0105149785839879,	
                       0.00218019738963307,	0.0456910986540881,	0.0104478762188781,	0.00401368340397422,	
                       0.00730663870911504,	0.023298823154929,	0.0107015634053377,	0.00103330643922971,	
                       0.0581718052613227,	0.000654558383600541,	0.00948093633599099,	0.0044619730170353,	
                       0.00443185209965115,	0.0438329966536452,	0.00128510732841969,	0.0186637033702838,	
                       0.0189128843105683,	0.000782830924866229,	0.000606330262480037,	0.000295001489949703,	
                       0.000215623917656018,	0.00158018968433101,	0.000147407813416395,	0.000237149684423799,	
                       0.00355225884350056,	0.00106852072890002,	0.0133842018703051,	0.00202939868214987,	
                       0.0100022851157315,	0.0156445997964005,	0.00141095871180609,	0.000458028063134958,	
                       0.0384318660871084,	0.00401580556437118,	0.00723052019450804,	0.000550816058070531,	
                       0.000172170281582251,	0.0341039154859239,	0.000539327580624991,	0.0152880658134926,	
                       0.00138236606690595,	0.000149644139047909,	0.00493513913903684,	0.00150414535765117,	
                       0.00105054602402996,	0.0380506352157902,	0.000580781814885316,	0.00129748440334172,	
                       0.0174911446260907,	0.00033920199119707,	0.00039897076520757,	0.0101722381301663,	
                       0.000598992521060009,	0.000282084099884211,	0.00638227700659961,	0.0002157345899337,	
                       0.00131762453051479,	0.00747180138080876,	0.00412074175755651,	0.000515966504344526,	
                       0.000285274885875181,	0.0021738913063059,	0.00115193387456638,	0.00162861682827093,	
                       0.0145438270323471,	0,	0.00823918187419527,	0.000186200009552022,	
                       0.000953999591325476,	0.000515966504344526,	0,	0,	0,	0,	
                       0,	0,	0,	0,	0,	0,	0,	0,	
                       0,	0,	0,	0,	9.01322257496156E-05,	4.46230840090351E-05,	
                       6.13662529511216E-05,	0,	0.000281620558518321,	0,	0,	
                       0.000118792923403944,	0.000164074276500117,	0,	1.83210862646205E-05,	0,	
                       0,	0)

# Combine East Chicago soil and other PCB profiles
prof.f <- data.frame(cbind(t.prof, s.prof))
t.prof.f <- t(prof.f)

# Perform PCA all samples
PCA <- prcomp(t.prof.f)

# Plot PCA results (Figure 8)
ggplot2::autoplot(PCA, data = t.prof.f, frame = TRUE,
                  frame.type = "t", size = 1.5, shape = 21) +
  stat_ellipse(color = "red") +
  theme_bw() +
  xlim(-0.4, 0.61) +
  ylim(-0.4, 0.61) +
  theme(aspect.ratio = 10/10) +
  annotate("text", x = 0.54, y = -0.135, label = "Aroclor 1221*",
           size = 3.5, fontface = 1) +
  annotate("text", x = 0.54, y = -0.09, label = "Aroclor 1221",
           size = 3.5, fontface = 1) +
  annotate("text", x = 0.35, y = 0.024, label = "Aroclor 1232*",
           size = 3.5, fontface = 1) +
  annotate("text", x = 0.14, y = 0.3, label = "Aroclor 1016",
           size = 3.5, fontface = 1) +
  annotate("text", x = 0.13, y = 0.27, label = "Aroclor 1242",
           size = 3.5, fontface = 1) +
  annotate("text", x = -0.08, y = 0.26, label = "Aroclor 1248",
           size = 3.5, fontface = 1) +
  annotate("text", x = -0.02, y = -0.29, label = "Aroclor 1260*",
           size = 3.5, fontface = 1) +
  annotate("text", x = 0.1, y = -0.06, label = "Aroclor 1262*",
           size = 3.5, fontface = 1) +
  annotate("text", x = -0.15, y = -0.08, label = "Aroclor 1254*",
           size = 3.5, fontface = 1) +
  annotate("text", x = -0.145, y = -0.1, label = "Aroclor 1254",
           size = 3.5, fontface = 1) +
  annotate("text", x = 0.02, y = -0.134, label = "EPA soil",
           size = 3.5, fontface = 1) +
  annotate("text", x = -0.07, y = 0.3, label = "Emission IHSC",
           size = 3.5, fontface = 1) +
  annotate("text", x = -0.135, y = -0.25, label = "CR soil",
           size = 3.5, fontface = 1) +
  annotate("text", x = 0.1, y = 0.07, label = "EC air",
           size = 3.5, fontface = 1) +
  annotate("text", x = 0.11, y = 0.17, label = "Volatilization LM",
           size = 3.5, fontface = 1)

# Perform PCA removing Aroclors 1221, 1232 and 1262 (Figure S4)
# Remove Aroclors 1221, 1232
prof.f.2 <- subset(prof.f, select = -c(X1221.Koh, X1221.Frame, X1232.Frame,
                                       X1262.Frame))
t.prof.2 <- t(prof.f.2)

# Perform 
# PCA 2
PCA.2 <- prcomp(t.prof.2)

# PCA plot (Figure S4)
ggplot2::autoplot(PCA.2, data = t.prof.2,
                  frame = TRUE, frame.type = "t", size = 1.5,
                  shape = 21) +
  stat_ellipse(color = "red") +
  theme_bw() +
  theme(aspect.ratio = 10/10) +
  annotate("text", x = 0.3, y = -0.32, label = "Aroclor 1016",
           size = 3, fontface = 1) +
  annotate("text", x = 0.24, y = -0.25, label = "Aroclor 1242",
           size = 3, fontface = 1) +
  annotate("text", x = 0.27, y = 0.03, label = "Aroclor 1248",
           size = 3, fontface = 1) +
  annotate("text", x = -0.22, y = -0.4, label = "Aroclor 1260*",
           size = 3, fontface = 1) +
  annotate("text", x = -0.057, y = 0.255, label = "Aroclor 1254*",
           size = 3, fontface = 1) +
  annotate("text", x = -0.08, y = 0.285, label = "Aroclor 1254",
           size = 3, fontface = 1) +
  annotate("text", x = -0.095, y = -0.17, label = "EPA soil",
           size = 3, fontface = 1) +
  annotate("text", x = 0.3, y = -0.15, label = "Emission IHSC",
           size = 3, fontface = 1) +
  annotate("text", x = -0.25, y = 0.0, label = "CR soil",
           size = 3, fontface = 1) +
  annotate("text", x = 0.1, y = -0.105, label = "EC air",
           size = 3, fontface = 1) +
  annotate("text", x = 0.23, y = -0.055, label = "Volatilization LM",
           size = 3, fontface = 1)

