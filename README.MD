Aggregated Wasserstein Metric and State Registration for Hidden Markov Models
===========
Aggregated Wasserstein is an alternative to traditional KL based distance for measuring the dissimilarity between hidden Markov Models.

Please see our ECCV 2016 (spotlight) paper for technical details:

[Yukun Chen](http://www.personal.psu.edu/yzc147), [Jianbo Ye](http://www.personal.psu.edu/jxy198), [Jia Li](http://stat.psu.edu/~jiali/), "A Distance for HMMs based on Aggregated Wasserstein Metric and State Registration".

```
@inproceedings{chen2016distance,
  title={A Distance for HMMs based on Aggregated Wasserstein Metric and State Registration},
  author={Chen, Yukun and Ye, Jianbo and Li, Jia},
  booktitle={Proceedings of The 14th European Conference on Computer Vision (ECCV)},
  page={},
  year={2016}
}
```

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



----
Copyright [2016] [Yukun Chen, Jianbo Ye, Jia Li]

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
