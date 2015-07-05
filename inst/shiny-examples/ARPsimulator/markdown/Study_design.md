### Study design features

In addition to entering information about the characteristics of the behavior and how it responds to treatment, you will need to specify some basic information about the design of the single-case study you plan to conduct. Currently, two different study designs are supported: multiple baseline designs and treatment reversal designs. Select the design that you want to simulate from the options in the first input box in column 2 of the simulator. The required further information depends on which option you select. 

#### Treatment reversal designs

Treatment reversal designs begin with a baseline phase, followed by a treatment phase. In designs with more than one AB pair, this process is repeated---sometimes several times---where each repetition entails removing the treatment that was in place (a __return to baseline__) and then re-introducing it. This __within-case replication__ provides multiple opportunities to test the functional relationship between the treatment and the outcome. The most common and well-known type of treatment reversal is the ABAB design, which provides three opportunities to demonstrate the functional relationship. 

To simulate a treatment reversal design, you will need to specify three further pieces of information:

1. Input the desired __number of cases__ in the second input box of column two. The default value for treatment reversal designs is 1 case. Further cases are essentially exact replicates of the first case, so you can keep things simple by just looking at one case at a time. 
2. Input the __number of (AB) phases__ in the third input box of column two. The default value of 2 corresponds to an ABAB design; a value of 3 corresponds to an ABABAB design, etc.
3. Input the number of __sessions per phase__ in the fourth input box of column two. For simplicity, the simulator assumes that each phase of the design has an equal number of phases. The default value of 5 sessions per phase is the minimum length recommended by the What Works Clearinghouse standards for single-case designs. 

#### Multiple baseline designs

In contrast to treatment reversal designs, multiple baseline designs use __cross-case replication__ in order to test the relationship between the treatment and the outcome, with each case providing one opportunity to demonstrate a functional relationship. The study begins by measuring each case in a baseline phase. The treatment is then introduced to each case, but the point of introduction is staggered so that not every case changes phases at the same time. Staggering the introduction of treatment allows the researcher to control for other common influences on the set of cases, which might be confounded with the introduction of the treatment if all cases entered treatment simultaneously. 

To simulate a multiple baseline design, you will need to specify three further pieces of information:

1. Input the desired __number of cases__ in the second input box of column two. The default value for multiple baseline designs is 3 cases, which is the minimum number recommended by the What Works Clearinghouse standards for single-case designs. 
2. Input the __total number of sessions__ per case in the third input box of column two. Note that the total number of sessions includes both baseline- and treatment-phase sessions. 
3. Input the __phase change times__ for each case in the fourth input box of column two. These values are equivalent to the lengths of the baseline phase for each case, and will usually be different for each case in the study. By default, the phase change times will be spread as evenly as possible across the total number of sessions in the study. For example, in a study with 3 cases and 20 sessions, the first case would enter treatment after 5 sessions, the second case would enter after 10 sessions, and the third cases would enter after 15 sessions.

#### Other designs

Currently, the ARPsimulator only provides options for the treatment reversal and multiple baseline designs. However, an alternating treatment design could be simulated by choosing the treatment reversal design and setting the number of sessions per phase to be one. (Of course, this approach does not allow you to randomize the order of treatments, as is sometimes done in alternating treatment designs.) Other types of designs, such as changing criterion designs or hybrid/complex designs, are not currently supported because they would introduce a number of additional complexities into the model. Feel free to contact me if you have suggestions for additional design features.  
