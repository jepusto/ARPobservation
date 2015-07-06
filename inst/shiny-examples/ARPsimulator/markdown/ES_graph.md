### Effect size graph

Another way to think about the results of a single-case study is by considering an estimate of __effect size__. Effect sizes are quantitative measures of the magnitude of treatment effects, which reflect the degree of change in an outcome for a given case as a result of intervention. A wide variety of effect size metrics has been proposed for use with single-case designs, as described in detail below. For any particular study, an effect size can only be estimated (as oppose to observed with 100% certainty) because the outcome measurements are influenced by chance fluctuations. However, since the simulator can produce hypothetical data from many identical repetitions of a study, it is possible consider the sampling distribution of the effect sizes. The __sampling distribution__ is a way of summarizing the range of possible effect size estimates that could be observed in the study, given the specified behavioral parameters, study design, and measurement procedures.

The __Effect sizes__ tab in the lower pane of the simulator displays the sampling distribution of an effect size, given a set of assumptions as specified in the input boxes of the upper pane. Initially, no graph will be displayed. You should begin by examining and modifying the options in the left-hand panel, which are described further below. Then hit the __Simulate!__ button to produce a graph of the estimated sampling distribution of the effect size estimate. This graph is a __density plot__, or smoothed histogram, which is a common way of representing a sampling distribution. Separate density plots will be displayed for each case in the study. For a given case, the horizontal axis of the graph corresponds to the range of possible values for the effect size estimate. The vertical height of the density corresponds to the relative frequency with which a given value of the effect size estimate is obtained. For example, suppose that the height of the density at an effect size of 80 is twice the height at an effect size of 45. This means that you are twice as likely to observe an effect size estimate around 80 as you are to observe an effect size around 45. Also, the area of the density plot is proportional to the probability of obtaining an effect size estimate in a given range. 

#### Options

The effect size sampling distribution tab has four further options:

1. 
2. 
3. The number of __samples per case__ controls how many hypothetical studies will be simulated in order to estimate the sampling distribution of the effect size for each case. 
4. 

#### Effect size measures

Currently, seven different effect size measures are available, including many of the non-overlap measures as well as the within-case standardized mean difference statistic. These seven measures were selected because they are commonly used as effect sizes for single-case designs. Other effect size measures that account for time trends and auto-correlation are excluded because the basic model embedded in the ARPsimulator assumes stability of baseline trends and independence of measurements across sessions. The seven included effect sizes are:

* __Percentage of non-overlapping data (PND)__
* __Percentage exceeding the median (PEM)__
* __Percentage of all non-overlapping data (PAND)__
* __Improvement rate difference (IRD)__
* __Non-overlap of all pairs (NAP)__
* __Tau__
* __Within-case standardized mean difference (SMD)__