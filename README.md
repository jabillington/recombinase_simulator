# Recombinate
Project to design an in-silico DNA recombinase simulator. Given two DNA sequences as genbanks (.gb) or FASTA (.fa) files, this program should identify any recombinase att sites within them, perform recombination in silico and output a new DNA sequence file of the recombined sequence.

Two methods are provided for doing so,

recombinase_site_detector- for two just sequences of interest

Usage: 

python **recombinase_site_detector.py** -i  inputfile1 inputfile2 outputfile

batch_recombine - for recombining multiple sequences with a single user specified template

Usage:

python **batch_recombine.py** -i path/to/batch common_template
  


![alt text](https://github.com/jambomber/recombinase_simulator/blob/main/recombinase_icon_low.png)
