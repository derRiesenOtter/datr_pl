# README.md DATR PL

The data analysis was performed before some changes to the exercise were made. To have comparable
results, the result data was taken from another student, where the data analysis was done after
all changes were made to the exercises: `/media/BioNAS/ag_hallab/DATR/Beck/`. The original
results for my analysis can be found here: `/media/BioNAS/ag_hallab/DATR/ender/`

The pipeline of exercise 6 was carried out on the bio cluster. The reports of exercise 6
and exercise 10 and 11 were carried out on the private computer.

## Exercise 6

Running following commands in this order from within the folder `/media/BioNAS/ag_hallab/DATR/Beck/methods`
will yield the results for this analysis:

```sh
sbatch 1_gffread.sh 
sbatch 2_fastq_cn1.sh
sbatch 2_fastq_cn2.sh
sbatch 2_fastq_ctl.sh
sbatch 2_trimmomatic_cnd2_SRR9929272.sh
sbatch 2_trimmomatic_cnd2_SRR9929279.sh
sbatch 2_trimmomatic_cnd2_SRR9929281.sh
sbatch 2_trimmomatic_cnd_SRR9929265.sh
sbatch 2_trimmomatic_cnd_SRR9929271.sh
sbatch 2_trimmomatic_cnd_SRR9929274.sh
sbatch 2_trimmomatic_ctl_SRR9929263.sh
sbatch 2_trimmomatic_ctl_SRR9929264.sh
sbatch 2_trimmomatic_ctl_SRR9929273.sh
sbatch 3_kallisto_index.sh
sbatch 3_kallisto_quant_cn1.sh
sbatch 3_kallisto_quant_cn2.sh
sbatch 3_kallisto_quant_ctl.sh
sbatch 3_bowtie2_indexing.sh
sbatch 3_bowtie2_align_cn1.sh
sbatch 3_bowtie2_align_cn2.sh
sbatch 3_bowtie2_align_ctl.sh
sbatch 3_convert.sh
sbatch 3_report_bowtie.sh
sbatch 3_stringtie.sh
```

### Exporting the results

All files needed for the further analysis were copied to a folder export and
exported to the private computer via:

```sh
scp -r studdatr@143.93.91.124:/media/BioNAS/ag_hallab/DATR/Beck/export ../material
mv ../material/export/* ../material
```

## Reports for exercise 6

To generate all report results of exercise 6 run the following command(s) in the methods directory:

```sh
Rscript 6_1.R > ../results/6_1.R
Rscript 6_2.R
Rscript 6_3_4.R > ../results/6_3_4.R
```

## Exercise 10

To generate all report results of exercise 10 run the following command(s) in the methods directory:

```sh
./download_go_annotations.sh
Rscript 10_tx2gene.R
Rscript 10.R > ../results/10.R
```

## Exercise 11

To generate all report results of exercise 11 run the following command(s) in the methods directory:

```sh
Rscript 11.R > ../results/11.R
```
