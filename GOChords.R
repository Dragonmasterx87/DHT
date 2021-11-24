# Installation of the latest released version
install.packages('GOplot')
library(GOplot)
packageVersion("GOplot")

gene.data <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\Coding Scripts\GOChord\Data\Massspec_data\UP\beta.DHT.de.data.csv)")
up.go <- read.csv(r"(C:\Users\mqadir\Box\Lab 2301\Coding Scripts\GOChord\Data\Massspec_data\UP\DHTGOup.csv)")
# down.go <- read.csv("D:/R-Projects/DHT/Data output/Down/DHTGOdown.csv")

head(gene.data)
head(up.go)
circ <- circle_dat(up.go, gene.data)
circ
process <- List("dephosphorylation of RNA polymerase II C-terminal domain", "actin filament capping", 
                "adenylate cyclase-modulating G protein-coupled receptor signaling pathway")
             
process
chord <- chord_dat(data = circ, genes = gene.data)
chord <- chord_dat(data = circ, process = process)
#chord <- chord_dat(data = circ, genes = gene.data, process = process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)


GOBubble(circ, labels = 3)
GOBubble(circ, title = 'Bubble plot', colour = c('darkgreen', 'red', 'darkblue'), display = 'multiple', labels = 3) 

# Reduce redundant terms with a gene overlap >= 0.75...
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# ...and plot it
GOBubble(reduced_circ, labels = 2.8)
GOBubble(reduced_circ, title = 'Bubble plot', colour = c('darkgreen', 'red', 'darkblue'), display = 'multiple', labels = 3) 
