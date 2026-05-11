#!/bin/bash

set -x
set -e

j=$(sbatch --parsable make_bowtie2_db.sbatch)
j=$(sbatch --parsable --dependency afterok:${j} bowtie2.sbatch)
j=$(sbatch --parsable --dependency afterok:${j} woltka.sbatch)
j=$(sbatch --parsable --dependency afterok:${j} make_q2_import.sbatch)
