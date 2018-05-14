# Robust and automatic MoCap data recovery (Matlab)

Matlab implementation of an original method for automatic and robust MoCap data recovery, based on a probabilistic averaging of different individual MoCap data recovery models.

This code requires the MoCap Toolbox for Matlab available at https://www.jyu.fi/hytk/fi/laitokset/mutku/en/research/materials/mocaptoolbox :

Burger, B. & Toiviainen, P. (2013). MoCap Toolbox – A Matlab toolbox for computational analysis of movement data. In R. Bresin (Ed.), Proceedings of the 10th Sound and Music Computing Conference, (SMC). Stockholm, Sweden: KTH Royal Institute of Technology.

This repository also provides an extension to the MoCap Toolbox.

## Installation

1. Download the repository (MocapRecovery)
2. Download the [MoCap Toolbox](https://www.jyu.fi/hytk/fi/laitokset/mutku/en/research/materials/mocaptoolbox)
3. Move the folder "mocaptoolbox" from the MoCap Toolbox to the MocapRecovery folder
4. Follow the script Example.m

## Method summary

![Method block diagram](diagram.jpg)

The above diagram shows the global approach of our data recovery method, which can be divided in different steps. First, parameters are extracted from each marker trajectory of the motion sequence. These parameters mainly represent relations between markers, and allow identification of neighboring markers and their distance distributions. These parameters are needed in the following steps. Then, different automatic individual recovery methods, based either on interpolation or machine learning (regression) techniques, are applied on the incomplete MoCap sequence, resulting in several temporary recovered sequences. For each recovered sequence, a correction is applied on each recovered trajectory to respect movement continuity (simple linear ramp for each coordinate). It allows to remove the trajectory disconinuity eventually occurring at each gap border. From all individually recovered sequences, a weighted average is then applied, as inspired from ensemble learning systems. The weight of each individual recovered sequence is computed based on neighbor markers distances distributions. This allows to give more credit to an individual method giving a more coherent recovered trajectory, providing a more robust final reconstruction, independently of any factor, including gap length, the number of available neighbor markers, or the length of the sequence. Finally, a space constraint is applied on the recovered trajectory, enforcing plausibility of its distance with neighboring markers. Basically, if the distance to a neighbor marker is outside its confidence interval, the recovered marker is translated to the limit of that confidence interval.

## Videos

Examples of MoCap data recovered with this method are displayed in the videos below (generated with the [MotionMachine](https://github.com/numediart/ofxMotionMachine)):

* Rolling motion capture sequence recovery: https://youtu.be/bgRZOWVPcWU
* Dancing and falling motion capture sequence recovery: https://youtu.be/Mq0NAZilIYY
* Karate motion capture sequence recovery: https://youtu.be/BwDAmHU8MKA

## Reference

An article describing and validating this method is under submission (PLoS ONE). If accepted, the reference will be provided here, and should be used to cite this work.

This repository can also be cited as: Mickael Tits. (2018, February 22). titsitits/MocapRecovery: First release (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.1182937
