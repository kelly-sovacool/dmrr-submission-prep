# dmrr-submission-prep
Prepare data for submission to the exRNA Data Coordination Center.
This project was created for the Tewari Lab at the University of Michigan.

## Usage

Create a configuration file in YAML format for your submission. Example:

```
group: exrna-mtewa1
user_login: sovacool
study_name: U01 Feeding Study July 2018
samples_filename: sample_sheet.csv
sample_study_names:
  - Feeding Study-30 min
  - Feeding Study-1 hr
study_id: EXR-MTEWA1FeedingStudy
working_dir: feedingStudy
templates_dir: templates
is_time_series: True
database: hg19_exrna
md5sum: daec25d670e3bb6b3ab3bbf5733df68c
fastq_filenames: HealthyControl_and_FeedingStudy_fastq_file.names.txt
tar_archive: None
```

Run the program with:
```
$ dmrr_prep.py config.yaml
```

The program will create donors metadata, biosamples metadata, and a manifest and write them to `working_dir`.
