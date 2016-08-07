Aggregated Wasserstein Metric and State Registration for Hidden Markov Models
===========
Aggregated Wasserstein is an alternative to traditional KL based distance for measuring the dissimilarity between hidden Markov Models.

Please see our ECCV 2016 (spotlight) paper for technical details:

[Yukun Chen](http://www.personal.psu.edu/yzc147), [Jianbo Ye](http://www.personal.psu.edu/jxy198), [Jia Li](http://stat.psu.edu/~jiali/), "A Distance for HMMs based on Aggregated Wasserstein Metric and State Registration".

The list of motions we used for our motion retrieval and motion classification experiments is provided in *index_to_use.txt*

## Usage (tested on MATLAB 2015a for Mac OSX)
Edit line 26 in initialization_script.m to set the path to the [mosek](https://www.mosek.com). (You can refer to [mosek installation guide](http://docs.mosek.com/7.0/toolsinstall.pdf), [mosek MATLAB mannual](http://docs.mosek.com/7.0/toolbox/Installation.html) for more detailed instruction on mosek installation) Then run the following at the folder initialization_script.m resides,
```
clear;
initialization_script;
```
To reproduce Figure 1 in the paper, run:
```
syntheticdata_exp_mu;
syntheticdata_exp_sigma;
```
To reproduce Figure 3 in the paper (and Figure 6 in the supplement material), run:
```
perturbation_exps_mu;
perturbation_exps_sigma;
perturbation_exps_transmat;
```

Or to reproduce all above, run:
```
run_all;
```
All resulting figures are also saved in the *imgs* folder.
