# IndependentComponentAnalysis
 This repository contains the codes used for the following work: 
 ```
 Liyan Song, Shuo Zhou, and Haiping Lu. "Direct ICA on Data Tensor via Random Matrix Modeling".
 ```
 

The code only requires basic Matlab libraries and a tensor toolbox, which can all be downloaded with a purchased Matlab licenses.


There is one coding folder "matlab/" that contains the codes implementing the proposed method, namely RAMICA.
To start the code, one needs to configure the files and paths by running the script "config.m" or type in the command window the following line:
```
>> config()
```
This function configures the directories between code scripts and the data set and among all scripts, so that the scripts can call each other and load the data set as if they are in a one-layer directory.


A quick start of running the implementation of the proposed RAMICA can refer to the script "example_run_SynB_RVM.m" in the "matlab/" folder, by typing in the command window the following line

    >> example_ramica()


For example, if the experimental setup is as below:

```
ica_name_cell = {'infomax', 'FastICA', 'jade', 'ramica'};
syn_data_name = 'tdistr'; 
nb_source = 4; 
nb_dimension = 32; 
nb_sample = 64;
seed_array = 1:10;
noise_level = 0;
```

The output should be expected as:

```
p_noise is 0.000, the amari errors of ICA methods are as below

amari_errors_ave =

    8.9652
    6.4198
    6.1490
    4.6827

Suceed!
```

Of course, people are free to change the experimental setup freely to play with this toy. 
For example, if one wants to see the recovery results when suffering from the label noise level 0.1, one can set the experiment as

```
noise_level = 0.1;
```

Or if one want to conduct more runs than 10, the below setup can be set:
```
seed_array = 1:100;
```


<p align="right">Enjoy~</p>

<p align="right">Liyan Song, June 2021</p>
<p align="right">Email: songly@sustech.edu.cn</p>
<p align="right">Southern University of Science and Technology, China</p>
