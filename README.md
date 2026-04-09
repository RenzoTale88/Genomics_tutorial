# Introduction to genomics
This repository contains the material for the genomics classes delivered 01-09 April 2026 at the University of Messina.
The material is partially derivative from the classes delivered by Dr. Andrea Talenti and Prof. James Prendergast for the [GCRF-STAR masterclass in Genomics](https://github.com/evotools/GCRF_tutorial) delivered in May 2022 at the Nelson Mandela Institute (Arusha, TZ).

The slides of the lectures are available at the following links:
1. [01/04/2026 - History of genomics](https://www.dropbox.com/scl/fi/7qvccgf0jfz87l9ctf0lz/day1_history_genomics.pdf?rlkey=09b1s8f6k85xpgmim87cktgvn&st=ad5zza9a&dl=1)
1. [02/04/2026 - Next Generation Sequencing and Variant calling](https://www.dropbox.com/scl/fi/165er33su7zwkbzbjdize/day2_shortread_analysis.pdf?rlkey=2wug8k6o9bh585lanim185wl5&st=01s8megh&dl=1)
1. [07/04/2026 - Genome assembly and Pangenomics](https://www.dropbox.com/scl/fi/q1w91xs9auimthvadksa5/day3_longread_pangenomics.pdf?rlkey=3tapnus8lsod7g9oq331d1xjt&st=817r9xue&dl=1)
1. [08/04/2026 - Epigenetics](https://www.dropbox.com/scl/fi/1l39vyn5c788kr2qjaz06/day4_epigenetics.pdf?rlkey=nqa3ch8xq00rppax88s63xiuy&st=j44lrtd3&dl=1)

Additionally, the repository contains the practicals for the following modules:
1. [Introduction to BASH](BASH/BASH_tutorial.md)
1. [Variant calling](VariantCalling/WGS_tutorial.md)
1. [Assembly and pangenomes](Pangenomes/Assembly_tutorial.md)
1. [RNA-seq analysis](Epigenetic/RNA_seq.md)

> **Note 1 for MacOS users**: some of the tools used in this tutorial will not be available to run on the Apple Silicon chips (M1 to M5) due to the ARM64 architecture and lack of compiled binaries. For these cases, users will have to install each dependency manually from the sources or attempt compiling.

> **Note 2 for MacOS users**: a possible solution is the use of containers, such as [docker](https://www.docker.com/), [apptainer](https://apptainer.org/) or [podman](https://podman.io/); while these are good options, the limited time available makes them inmprossible to adopt for this course.

## Credits
This repository uses material created by:
1. [Andrea Talenti](https://www.gla.ac.uk/schools/bohvm/staff/andreatalenti/)
1. [James Prendergast](https://edwebprofiles.ed.ac.uk/profile/james-prendergast)