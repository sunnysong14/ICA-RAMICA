[2018/10/15] I have tried to clean the codes and comments as much as possible though I did not go through all of them. 

The below describes how to use the codes in the 'matlab' folder.

o To configurate the experimental environment, please run MAKEUP() in the first place.

o Script script_Exp1.m is for the first experiments, i.e. BSS on the synthetic data. 
    * mainSyn_fun() is the main function called in this script to compute and save our experimental results. 

o Script script_Exp2.m is for the second experiments, i.e. blind image separation
    * mainIS_fun() is the main function called in this script to compute and save our experimental results. 

o Folder 'ICAs' contains the ICA methods including FastICA, JADE, Informax, and our proposed RAMICA. 
    * The subfolder 'basic' contains the common functions used among ICA methods.

o Folder 'basic' contains the basic function for conducting our experiments, such as the generation matrix, how to inject Gaussian noise to the data, how to compute Amari error and so on. 

>>>> Best wishes for the acceptance! 
