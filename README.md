# Introduction to genomics
This repository contains the material for the genomics classes delivered 01-09 April 2026 at the University of Messina.
The material is partially derivative from the classes delivered by Dr. Andrea Talenti and Prof. James Prendergast for the [GCRF-STAR masterclass in Genomics](https://github.com/evotools/GCRF_tutorial) delivered in May 2022 at the Nelson Mandela Institute (Arusha, TZ).

This repository contains the practicals for the following modules:
1. [Introduction to BASH](BASH/BASH_tutorial.md)
1. [Variant calling](VariantCalling/WGS_tutorial.md)
1. [Assembly and pangenomes](Pangenomes/Assembly_tutorial.md)
1. [RNA-seq analysis](Epigenetic/RNA_seq.md)

> **Note 1 for MacOS users**: some of the tools used in this tutorial will not be available to run on the Apple Silicon chips (M1 to M5) due to the ARM64 architecture and lack of compiled binaries. For these cases, users will have to install each dependency manually from the sources or attempt compiling.
> **Note 2 for MacOS users**: a possible solution is the use of containers, such as [docker](), [apptainer]() or [podman](); while these are good options, the limited time available makes them inmprossible to adopt for this course.

## Credits
This repository uses material created by:
1. [Andrea Talenti](https://www.gla.ac.uk/schools/bohvm/staff/andreatalenti/)
1. [James Prendergast](https://edwebprofiles.ed.ac.uk/profile/james-prendergast)