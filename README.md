# Overview
**Rawasm** is the first software tool that enables the construction of genome assembly from raw nanopore signals.
It mostly reuses the **[miniasm](https://github.com/lh3/miniasm)** features, but adds support to _FAST5_, _POD5_ and _SLOW5_ formats.
Rawasm can be used in pipelining with **[RawHash2](https://github.com/CMU-SAFARI/RawHash)**, using **Rawsamble** overlapping feature.

# Installation

To install **Rawasm**, do the following:
* Clone the repository
  ```bash
  git clone https://github.com/CMU-SAFARI/rawasm.git rawasm
  ```
* install Rawasm
  ```bash
  cd rawasm && make
  ```

**Rawasm** makefile downloads a new **miniasm** and patches it. 
Default installation directory is _rawasm/miniasm_, but you can change it in the Makefile.
If a local copy of **miniasm** exists, you can specify it in the Makefile.
If this is the case, to skip the miniasm download and install use:
```bash
make install
```

**Rawasm** is a self-contained implementation that can be downloaded and run. It uses a set of pre-compiled static libraries.
However, recompiling might be required depending on your system. If it is the case, you can choose two options:
- Refer to **[RawHash2](https://github.com/CMU-SAFARI/RawHash)** repo for _SLOW5_, _POD5_ and _FAST5_ libraries compilation.
- Compile the libraries using _gcc ar_ from the original repos:
  - _libhdf5.a_ : compile from **[HD5 Group](https://github.com/HDFGroup/hdf5)** repo.
  - _libpod5_format.a, libarrow.a, libjemalloc_pic.a, libzstd.a_ : compile from **[POD5 format](https://github.com/nanoporetech/pod5-file-format)** repo.
  - _libslow5.a_ : compile from **[Slow5 Tools](https://github.com/hasindu2008/slow5tools)** repo.

Finally, **Rawasm** requires **[libuuid](https://github.com/cloudbase/libuuid/tree/master)**. you can do the following:
- Download libuuid
```bash
git clone https://github.com/cloudbase/libuuid/tree/master libuuid
```
- Build .o files (run the following for each .c or just write a Makefile/script to automate)
```bash
gcc -c -g -Wall -O2 -Wno-all  -Wno-write-strings -Wno-deprecated-declarations -Wcpp -I. file.c -o file.o
```
- Generate the library
```bash
ar rcs libuuid.a *.o
```
All static libraries must be in the _lib_ directory before running the _make_ or _make install_ command.
# Usage

**Rawasm** supports by default all of **miniasm** features. Moreover, it introduces two new features:
- Processing of single or multiple _FAST5, POD5, S/BLOW5_ input files
- Output the assembly as _FAST5, POD5, S/BLOW5_ unitigs files. For more info about unitigs, check **[miniasm](https://github.com/lh3/miniasm)**.

Using **Rawasm** is straightforward. The following com
```bash
./miniasm -f input_data[.fast5/pod5/slow5/blow5] overlaps.paf -H outdir > assembly.gfa
```
_input_data_ can be either a directory containing multiple files, or a single file (fast5, pod5, slow5, blow5). In case of a directory, do not mix different-type files.
_overlaps.paf_ is the all vs all overlaps file produced by **[RawHash2](https://github.com/CMU-SAFARI/RawHash)**.
_outdir_ specifies the unitig files output directory. **Rawasm** creates a unitig file for each distinct unitig that makes the assembly. The format type is the same of the input.
_assembly.gfa_ assembly text output.

# Cite Rawsamble
If you use Rawasm in your work, please consider citing the following papers:

```bibtex
TBD
```
