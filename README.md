# Identify buried residues in protein subunits and calculate its relative moment 

## Background knowledge

There are many protein complexes in nature with various structures and functions. They can be classified into two categories: 
- Homomeric complex: consists of multiple copies of the same protein chains
- Heteromeric complex: consists of multiple copies of different protein chains

Understanding the protein complex assembly mechanism is fundamental, especially for the complicated protein complex machines that play important role in various biological process. Any assembly must start from two protein subunits, so here a two-subunit complex is the study object for the tools to run on.

To visulise protein complexes, the figure below shows three example that represent all conditions of heteromeric protein complex (from *Saccharomyces cerevisiae*):

1. **5W3F** --- Simple heteromeric complex: only contain two subunits
2. **3O8O** --- Complex heteromeric complex: contain multiple copy (4 at here) of the same two-subunit complex
3. **7QO4** --- Complex heteromeric complex: contain multiple copy of different subunits

A single two-subunit complex are highlighed in each model, with one subunit coloured in teal and another one coloured in pink. The buried residues in each subunit that contributed to form the protein-protein interface is highlighted in orange and red respectively. The background subunits are in white.

![PDB models](figures/PDB_models.svg)


## Description
This repository contain two Python scripts that provide tools to work on protein complex that contains at least two subunits.

1. **Tool CBR: cbr.py / cbr.ipynb** (cbr stands for Calculate Buried Residues)  
   It runs an external program "naccess" to calculate the Solvent Surface Accessability (ASA) of each residues in the amino acid sequence (chains) of the two subunits in a two-subunit complex, then identify the buried residues within each subunit, if there is any.  
   (This calculation is recommend to run the cbr.py script in commend line. The cbr.ipynb is only to shows steps)
   
1. **Tool CRM: crm.ipynb** (crm stands for Calculate Relative Moments)  
   It uses the output results from the first step to calculate a value called Relative Moment (M-rel) for each subunit in the a two-subunit complex.  
   (You need one cal.ipynb file for each subunit in one combination of two-subunit complex)


### 1. Tool CBR: cbr.py / cbr.ipynb

**Principle:**

This tool identify the buried residues in each subunit within a two-subunit complex, if there is any. The rational is illustrated by this schematic diagram.

For a two-subunit complex contain subunit 1 (blue) and subunit 2 (red), the (protein-protein) interface (yellow) comes from the overlap of the (buried) surface 1 (light blue) and the (buried) surface 2 (pink) from the corresponding subunits. 

The program "naccess" written by Simon Hubbard (1992) calculates the Absolute Solvent Accessability (ASA) (water, by default) of each residues in one or multiple chains environment (naccess produce more data but we only care about ASA here). 

This tool utilises the "naccess", for each two-subunit complex, calculate the ASA for each residues in the two 'single chain only' status and the one 'within a complex' status. If a residues contribute to the formation of the buried residues, they will be less accessible in the 'within a complex' status and therefore have smaller values. The residues do not contribute to the formation will be equivalently accessibly in both statuses and therefore have the same values. 

By subtracting the ASA values in 'within a complex' status using the two 'single chain only' status for each residues in the two chains, the buried residues will have a negative value, and therefore being determined.

![principle_1](figures/principle_1.png)

**Calculation process:**


![calculation_1](figures/calculation_1.png)

**Example result:**

 Using PDB model 5W3F as an example:

 - The output result will look like the data in the two tables. 
 - The bar plots highlight the location of the buried residues on the scale of the respective subunit protein sequence. 
 - The histograms show the size of the ASA values of these buried residues.


![result_1](figures/result_1.svg)



### 2. Tool CRM: crm.ipynb

**Principle:**

To determine the relative location of the buried surface of a subunit, a concept of Relative Moment (M-rel) of a subunit was raised as the overall contribution of two values:  
- (1) `the relative absolute ASA (∆a)`  --- reflects the contribution from each residue to the total ASA of the whole subunit  
- (2) `the relative distance (∆d)` --- applies a position weight to each residue when considering their contribution/effect to the total ASA of the whole subunit. 

The **N-terminal sequence** and the **C-terminal sequence** are defined as the two halves of a protein sequence that towards the different two ends when separated from the middle point. The effect of residues on the N-terminal of the sequence counteracts the effect of residues on the C-terminal. 

Since `∆d` is in the range from `-0.5 to 0.5` and `∆a` is in the range from `0 to 1`. The value of M-rel would be in the range of `-1 to 1`, which enables comparison between subunits from different complexes. The direction of the relative Moment (M-rel) (smaller or greater than 0) indicates whether the buried surface of a particular subunit is located on the N terminal or the C terminal. The more extreme the number, the more the buried surface localises towards the two ends in the protein sequence of a subunit. 

![principle_2](figures/principle_2.png)



**Calculation process:**


![calculation_2](figures/calculation_2.png)



**Example result:**

 Using PDB model 5W3F as an example:
 - The two extremes of each protein subunit have the M-rel values of -0.5 and 0.5. The location of M-rel = 0.0 is drown by dash lines in the middle. The calculated M-rel values (-0.01 and 0.09) of the two subunits are visualised as the red dots in the bar plots.

![result_2](figures/result_2.png)



## Requirements
- Python 3
- Naccess (download by link in references)
- Pandas library
- Biopython package

## Workflow
The figure shows the workflow of the two scripts.  
### CBR  
>**Section 1:** Download PDB files  
>**Section 2:** Calculate ASA values  
### CRM  
>**Section 3:** Mapping data  
>**Section 4:** Calculate (average) M-rel

![flowchart_simple](figures/flowchart_simple.svg)

In details, this flowchart shows the function of each step and the relevant input and output files.

![flowchart](figures/flowchart.svg)

## Usage Instructions
Plese see the repository wiki for detailed usage instructions
1. [Tool CBR](https://github.com/Luna120120/Calculate-ProteinSubunit-RelativeMoment/wiki/Usage-Instructions-%E2%80%90%E2%80%90%E2%80%90-1.-Tool-CBR:-cbr.py---cbr.ipynb)
2. [Tool CRM](https://github.com/Luna120120/Calculate-ProteinSubunit-RelativeMoment/wiki/Usage-Instructions-%E2%80%90%E2%80%90%E2%80%90-2.-Tool-CRM:-crm.ipynb)

---

## Reference
Naccess --- Solvent accessible area calculations (S. Hubbard and J. Thornton, 1992)
http://www.bioinf.manchester.ac.uk/naccess/nac_readme.html  
http://www.bioinf.manchester.ac.uk/naccess/nacdownload.html
