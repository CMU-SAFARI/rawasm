# Overview
**Rawasm** is the first software tool that enables the construction of genome assembly from raw nanopore signals.
It mostly reuses the **miniasm** (https://github.com/lh3/miniasm) features, but adds support to FAST5, POD5 and SLOW5 formats.
Rawasm can be used in pipelining with **RawHash2** (https://github.com/CMU-SAFARI/RawHash), using **Rawsamble** overlapping feature.

# Installation

To install **Rawasm**, do the following:
* Clone the repository
  ```bash
  git clone --recursive https://github.com/CMU-SAFARI/rawasm.git rawasm
  ```
* install Rawasm
  ```bash
  cd rawasm && make
  ```

Rawasm makefile downloads a new **miniasm** and patches it. 
Default installation directory is _rawasm/miniasm_, but you can change it in the Makefile.
If a local copy of **miniasm** exists, you can specify it in the Makefile.
If this is the case, to skip the miniasm download and install use:
```bash
make install
```




# Usage

# Cite Rawsamble
