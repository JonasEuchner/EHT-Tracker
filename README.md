# Contractility analysis of engineered heart tissues
This code was developed as a quick solution for a collaboration partner to extract meaningful parameters for EHT. As a result, standard coding practices were not strictly followed, and the code may lack documentation and optimization. Use at your own discretion. 

## Getting started
### ImageJ script
Open ImageJ and select `Plugins` -> `Marcos` -> `Run...`. Then navigate to the downloaded git repository and select `EHT-pillarTracker.ijm`. 
Select the folder where all tiff-stacks are located.
Now follow the instruction of the GUI (Use the Multi-Point Tool to select a feature of the pillar that should be tracked).
Add constraints to the per frame displacement in the X&Y-direction as well as the pixel size of the pattern to be followed

### R script
Open RStudio and run the script `ParameterExtraction_v2.0.R`. 
Select the folder where all "*.txt" files are located that were produced from the ImageJ script
Now provide the needed information into the pop-up windows.

## Additional info
Written  to measure contraction parameters of EHT<br>
- "MuscleMotion" by Sala et al. is unable to measure force during contraction (Circ. Res.: 10.1161/CIRCRESAHA.117.312067)
- "EHT analysis" by Rivera-Arbeláez et al. was not working as painting pillars with carbon black was avoided (PloS one: 10.1371/journal.pone.0266834)