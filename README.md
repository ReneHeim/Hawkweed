## Abstract

This code analysis spectral profiles of plant species from Mount Kosciuszko
National Park, NSW, Australia. Initially, outlying spectra from each species 
will be removed.Then, a Random Forest classifier will be trained and validated. Follwoing questions will be
addressed:

- Can spectral profiles of all species be accurately classified? 
- Is the classification still accurate if spectral profiles are resampled according to the specifications of following sensors:
  + Micasense Parrot Sequoia 
  + Sentinel-2

## Instructions

1. Please install R and R Studio.
2. Download and open the provided file _Analysis_scriptonly.R_
3. Uncomment and run line 18-28 to install required packages.
4. Download _R/_ folder where custom made functions are contained. Save this folder in you R project working directory.
5. Run entire code to reproduce analysis.
6. All results can be found in the _output/_ folder.
7. For further explanations, the _Analysis.Rmd_ can be rendered.

The use of this code must be cited:

Heim, 2018, Spectral classification of invasive Hawkweek (_Hieracium aurantiacum_) at Mt. Kosciuszko NP. https://github.com/ReneHeim/Hawkweed

## Data

Data will be provided on request.

## Licences

Data: CC-0 attribution requested in reuse
Manuscript: CC-BY-4.0
Code: MIT year: 2017, copyright holder: René Heim

Copyright © 2018 René Hans-Jürgen Heim

__contact: rene.heim@hdr.mq.edu.au__ 
