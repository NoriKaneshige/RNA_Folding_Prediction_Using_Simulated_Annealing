# RNA_Folding_Prediction_Using_Simulated_Annealing

## [Comments about the result]
### 1) I generated the average base-pair matrix over 200 iterations for both UUCG-tetraloop and Human Telomerase RNA.
### 2) When I compare each average base-pair matrix and NMR structures, I found discrepancies. My simulated_annealing model with 200 iterations and given hyperparameters could not really predict UUCG tetra loop in both RNAs and also failed in predicting GCUCC internal loop in Human Telomerase RNA.
### 3) A possible improvement would be to increase the iteration number. I used 200 iteration this time. Another improvement would be to tune hyperparameters such as distribution parameter.

## [Result]
## UUCG-tetraloop 2D structure generated in [RNApdbee](http://rnapdbee.cs.put.poznan.pl/) using PDB data.
![Image 1](UUCG_tetraloop_NMR.png)

## Vusualize UUCG-tetraloop 2D structure as a Line graph
![Image 2](UUCG_tetraloop_from_NMR.png)

## Averaged UUCG-tetraloop 2D structure predicted by Simulated Annealing
![Image 3](predicted_averaged_UUCG-tetraloop_2D_structure.png)

## Human telomerase RNA 2D structure generated in [RNApdbee](http://rnapdbee.cs.put.poznan.pl/) using PDB data.
![Image 1](HUMAN_TELOMERASE_RNA_NMR.png)

## Vusualize Human telomerase RNA 2D structure as a Line graph
![Image 2](UUCG_tetraloop_from_NMR.png)

## Averaged Human telomerase RNA 2D structure predicted by Simulated Annealing
![Image 3](predicted_averaged_Human_Telomerase_RNA_2D_structure.png)



