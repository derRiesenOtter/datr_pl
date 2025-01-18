# README.md DATR PL

The data analysis was performed before some changes to the exercise were made. To have comparable
results the results were taken from another directory, where the data analysis was done after
all changes were made to the exercises: `/media/BioNAS/ag_hallab/DATR/Beck/`. The original
results for my analysis can be found here: `/media/BioNAS/ag_hallab/DATR/ender/`

The pipeline of exercise 6 was carried out on the bio cluster. The reports of exercise 6
and exercise 10 and 11 were carried out on the private computer.

## Exercise 6

### Exporting the results

All files needed for the further analysis were copied to a folder export and
exported to the private computer via:

```sh
scp -r studdatr@143.93.91.124:/media/BioNAS/ag_hallab/DATR/Beck/export ../material
mv ../material/export/* ../material
```

## Reports for exercise 6

To generate all report results of exercise 6 run:

```sh
Rscript 6_1.R > ../results/6_1.R
Rscript 6_2.R
Rscript 6_3_4.R > ../results/6_3_4.R
```

## Exercise 10

All commands are called from within the `methods` directory if not stated otherwise.

```sh
./download_go_annotations.sh
Rscript 10_tx2gene.R
Rscript 10.R > ../results/10.R
Rscript 11.R > ../results/11.R
```

```

## Exercise 11
