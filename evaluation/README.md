# Evaluation

This folder contains scripts for evaluating results submitted to the 11th SPE CSP.

## check_format.py

Given a folder and case letter, the script checks if the result files present in the folder are conforming to the requirements as given in the [description](https://doi.org/10.2118/218015-PA). For example,
```bash
python3 ./check_format.py -f path_to_folder -c A
```
will check the files in `path_to_folder` for the case CSP 11A. If an archive `spe11a.zip` is found, the files in that archive will be checked.
