ARPobservation provides a set of tools for simulating the sort of data generated during direct observation of a behavior over time. It works by first simulating a behavior stream based on an alternating renewal process, based on specified distributions of event durations and interim times. The following distributions are currently supported:

* Exponential
* Gamma with fixed shape parameter
* Mixture of two gammas
* Uniform(0,x)
* Constant (degenerate)

Different methods of recording data can then be applied to the simulated behavior stream. Functions are provided for the following recording methods: 

* continuous duration recording 
* event counting
* momentary time sampling
* partial interval recording
* whole interval recording

