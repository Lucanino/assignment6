# You must install the package HMM require(HMM) 
library("HMM", lib.loc="~/R/win-library/3.6")

require(HMM)
#HMM Paket wird geladen
nSim = 2000 
#Variable für die Länge der Simulation
States = c("Fair", "Unfair") 
#Variable für die beiden Emissionszustände
Symbols = 1:6 
#Festlegung der Ergebnismöglichkeiten für Würfe
transProbs = matrix(c(0.99, 0.01, 0.02, 0.98), c(length(States), length(States)), byrow = TRUE) 
#Festlegung der Transitionswahrscheinlichkeiten als Matrix
emissionProbs = matrix(c(rep(1/6, 6), c(rep(0.1, 5), 0.5)), c(length(States), length(Symbols)), byrow = TRUE) 
#Festlegung der Emissionswahrscheinlichkeiten als Matrix 
hmm = initHMM(States, Symbols, transProbs = transProbs, emissionProbs = emissionProbs) 
#die oben definierten Parameter werden zu einem Modell verknüpft
sim = simHMM(hmm, nSim) 
#mit der Länge nSim wird das Modell hmm ausgeführt
vit = viterbi(hmm, sim$observation) 
#vit ist ein Vektor, der die versteckten Emissionszustände beinhaltet
f = forward(hmm, sim$observation) 
#f ist ein Vektor mit den akkumulierten Wahrscheinlichkeiten dafür den in sim$observation beobachteten Status jeweils im fairen und unfairen modus zu erhalten
i <- f[1, nSim] 
#ein Vektor aus den Wahrscheinlichkeiten den letzten zustand durch die kette der vorherigen zustände zu erreichen für den fairen wurf
j <- f[2, nSim] 
#ein Vektor mit den Wahrscheinlichkeiten den letzten zustand durch die kette der vorherigen zustände zu erreichen für den unfairen Wurf
probObservations = (i + log(1 + exp(j - i))) 
#
######################################### 
## NO MORE DOCUMENTATION BELOW THIS LINE 
######################################### 
x = list(hmm = hmm, sim = sim, vit = vit) 
# PLOT simulated throws at top ##################### 
mn = "Fair and unfair die" 
xlb = "Throw nr." 
ylb = "" 
plot(x$sim$observation, ylim = c(-7.5, 6), pch = 3, main = mn, xlab = xlb, ylab = ylb, bty = "n", yaxt = "n") 
axis(2, at = 1:6) 
# PLOT Simulated, which die was used (ground truth) ########### 
text(0, -1.2, adj = 0, cex = 0.8, col = "black", "True: green = fair die") 
for (i in 1:nSim) { 
  if (x$sim$states[i] == "Fair") 
    rect(i, -1, i + 1, 0, col = "green", border = NA) 
  else rect(i, -1, i + 1, 0, col = "red", border = NA) 
} 
# PLOT Most probable path (viterbi) ####################### text(0, -3.2, adj = 0, cex = 0.8, col = "black", "Most probable path") 
for (i in 1:nSim) { 
  if (x$vit[i] == "Fair") 
    rect(i, -3, i + 1, -2, col = "green", border = NA) 
  else rect(i, -3, i + 1, -2, col = "red", border = NA) 
} 
# PLOT Differences #################### 
text(0, -5.2, adj = 0, cex = 0.8, col = "black", "Difference")
differing = !(x$sim$states == x$vit) 
for (i in 1:nSim) {
  if (differing[i]) 
    rect(i, -5, i + 1, -4, col = rgb(0.3, 0.3, 0.3), border = NA) 
  else rect(i, -5, i + 1, -4, col = rgb(0.9, 0.9, 0.9), border = NA) 
}


#neues Programm Aufgabe 3

nSim = 2000 
#Variable für die Länge der simulierten Sequenz
States = c("codierend", "nichtcodierend") 
#Variable für die beiden Emissionszustände
Symbols = c(1:4) 
#Festlegung der Codes für die Nucleotide
transProbs = matrix(c(0.99, 0.01, 0.02, 0.98), c(length(States), length(States)), byrow = TRUE) 
#Festlegung der Transitionswahrscheinlichkeiten als Matrix (hier jetzt Wahrscheinlichkeiten von Intron auf Exon zu wechseln)
emissionProbs = matrix(c(rep(0.2,0.3,0.3,0.2), c(rep(0.25, 4))), c(length(States), length(Symbols)), byrow = TRUE) 
#Festlegung der Emissionswahrscheinlichkeiten für die Nucleotide als Matrix
hmm = initHMM(States, Symbols, transProbs = transProbs, emissionProbs = emissionProbs) 
#die oben definierten Parameter werden zu einem Modell verknüpft
sim = simHMM(hmm, nSim) 
#mit der Länge nSim wird das Modell hmm ausgeführt
vit = viterbi(hmm, sim$observation) 

x = list(hmm = hmm, sim = sim, vit = vit) 
# PLOT simulated throws at top ##################### 
mn = "Nukleotid" 
xlb = "Stelle in Sequenz" 
ylb = "" 
plot(x$sim$observation, ylim = c(-7.5, 6), pch = 3, main = mn, xlab = xlb, ylab = ylb, bty = "n", yaxt = "n") 
axis(2, at = 1:4) 
# 
text(0, -1.2, adj = 0, cex = 0.8, col = "black", "True: green = codierend") 
for (i in 1:nSim) { 
  if (x$sim$states[i] == "codierend") 
    rect(i, -1, i + 1, 0, col = "green", border = NA) 
  else rect(i, -1, i + 1, 0, col = "red", border = NA) 
} 
# PLOT Most probable path (viterbi) ####################### text(0, -3.2, adj = 0, cex = 0.8, col = "black", "Most probable path") 
for (i in 1:nSim) { 
  if (x$vit[i] == "codierend") 
    rect(i, -3, i + 1, -2, col = "green", border = NA) 
  else rect(i, -3, i + 1, -2, col = "red", border = NA) 
} 
# PLOT Differences #################### 
text(0, -5.2, adj = 0, cex = 0.8, col = "black", "Difference")
differing = !(x$sim$states == x$vit) 
for (i in 1:nSim) {
  if (differing[i]) 
    rect(i, -5, i + 1, -4, col = rgb(0.3, 0.3, 0.3), border = NA) 
  else rect(i, -5, i + 1, -4, col = rgb(0.9, 0.9, 0.9), border = NA) 
}

# Erstellung des plots

#andere Parameter

transProbs = matrix(c(0.999, 0.001, 0.01, 0.99), c(length(States), length(States)), byrow = TRUE) 
#Festlegung der Transitionswahrscheinlichkeiten als Matrix (hier jetzt Wahrscheinlichkeiten von Intron auf Exon zu wechseln)
emissionProbs = matrix(c(rep(0.1,0.4,0.4,0.1), c(rep(0.25, 4))), c(length(States), length(Symbols)), byrow = TRUE) 
#Festlegung der Emissionswahrscheinlichkeiten für die Nucleotide als Matrix
hmm = initHMM(States, Symbols, transProbs = transProbs, emissionProbs = emissionProbs) 
#die oben definierten Parameter werden zu einem Modell verknüpft
sim = simHMM(hmm, nSim) 
#mit der Länge nSim wird das Modell hmm ausgeführt
vit = viterbi(hmm, sim$observation) 

x = list(hmm = hmm, sim = sim, vit = vit) 
# PLOT simulated throws at top ##################### 
mn = "Nukleotid" 
xlb = "Stelle in Sequenz" 
ylb = "" 
plot(x$sim$observation, ylim = c(-7.5, 6), pch = 3, main = mn, xlab = xlb, ylab = ylb, bty = "n", yaxt = "n") 
axis(2, at = rep(1:4,1)) 
# 
text(0, -1.2, adj = 0, cex = 0.8, col = "black", "True: green = codierend") 
for (i in 1:nSim) { 
  if (x$sim$states[i] == "codierend") 
    rect(i, -1, i + 1, 0, col = "green", border = NA) 
  else rect(i, -1, i + 1, 0, col = "red", border = NA) 
} 
# PLOT Most probable path (viterbi) ####################### text(0, -3.2, adj = 0, cex = 0.8, col = "black", "Most probable path") 
for (i in 1:nSim) { 
  if (x$vit[i] == "codierend") 
    rect(i, -3, i + 1, -2, col = "green", border = NA) 
  else rect(i, -3, i + 1, -2, col = "red", border = NA) 
} 
# PLOT Differences #################### 
text(0, -5.2, adj = 0, cex = 0.8, col = "black", "Difference")
differing = !(x$sim$states == x$vit) 
for (i in 1:nSim) {
  if (differing[i]) 
    rect(i, -5, i + 1, -4, col = rgb(0.3, 0.3, 0.3), border = NA) 
  else rect(i, -5, i + 1, -4, col = rgb(0.9, 0.9, 0.9), border = NA) 
}

