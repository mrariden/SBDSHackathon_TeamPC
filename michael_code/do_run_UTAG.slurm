#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --output=result_%j.out
#SBATCH --partition=instructional
#SBATCH --mem=64Gb
#SBATCH --cpus-per-task=8
#SBATCH -A BDShackathon


/home/mr2ep/hack_env/bin/python run_UTAG.py
