#!/bin/bash

set -x
set -e

j=$(sbatch --parsable make_q2_import.sbatch)
j=$(sbatch --parsable --dependency afterok:${j} deblur.sbatch)
