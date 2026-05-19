# Original model codes
These files are from the study by Chandrasekaran and Jiang (2019).  
Original GitHub repository: https://github.com/JiangLabUCI/AbResistantDoseResponse 


# Modified model code
This is a simplified version of the original model code that produces outputs similar to the original implementation. The code is more compact
and combines the analysis and plotting steps into a single script. However, sensitivity analyses and model comparison components from the original
implementation are not included. In addition, the code directly uses the reported parameter values for the exponential and Beta-Poisson models,
making the workflow simpler and easier to understand.


# Final code for Vibrio
This code follows a similar structure to the modified model code but is adapted for _Vibrio_. The model parameters were adjusted accordingly, 
and two antibiotics were incorporated instead of one. Additionally, only the Beta-Poisson model was used, as it was identified as the 
best-fitting model based on the literature. The implementation was verified to produce results comparable to those reported in 
Chandrasekaran and Jiang (2019).
