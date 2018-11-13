# DAMACY operator

[Discriminant Analysis of Multi-Aspect CYtometry](https://www.ru.nl/science/analyticalchemistry/education/master-internships/vm/increasing-power-damacy/)

Integrating DAMACY sofware into Tercen.

## Overview

In current version of Tercen operators are coded using R, so matlab code requires to be wrapped into an R operator.

```
Tercen ==> DAMACY R operator ==> DAMACY matlab
```

## DAMACY R operator

- extract data from tercen
- transform data to fit DAMACY format
- extract and transform parameters to fit DAMACY format
- call DAMACY matlab function
- load and transform results to fit tercen format
- save result into tercen

## DAMACY matlab

The entry point is a function with the following signature :

```
DAMACY(String dataFile, String parametersFile, String outputFolder)
```

2 files, one containing the data, and one containing the parameters.

1 folder, contains DAMACY results.

## File format

### Data

format : csv

3 types of columns :

- proteins : a set of columns
- label : control = 1 not = 0
- patientID : a string identifying a patient

Each line represent a cell.

| prot1 | prot2 | prot3 | prot4 | label | patientID |
| ----- | ----- | ----- | ----- | ----- | --------- |
| | | | | | |
| | | | | | |

### Parameters

format : csv

2 columns :

- name : parameter name
- value: parameter value

| name | value |
| ---- | ----- |
|||
|||

### Output data

DAMACY result is a set of 3 tables, stored in csv files.

- Prediction

filename : prediction.csv

| patientID | predictionValue |
| --------- | --------------- |
|||
|||

- Loadings

filename : loadings.csv

| protein | loadingX | loadingY |
| ------- | -------- | -------- |
||||
||||

- Cell map

filename : cellmap.csv

| cellmapX | cellmapY | cellmapValue |
| -------- | -------- | ------------ |
||||
||||


## Packaging and deployment

Tercen R operators are packaged and deployed using git repository.

https://github.com/tercen/damacy_operator
 

## Links

https://www.ru.nl/science/analyticalchemistry/education/master-internships/vm/increasing-power-damacy/

https://www.ru.nl/science/analyticalchemistry/research/software/
