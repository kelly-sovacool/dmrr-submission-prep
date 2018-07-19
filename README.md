# DMRR Submission Preparation

Prepare data for submission to the [exRNA Data Coordination Center](http://genboree.org/theCommons/projects/exrna-mads/wiki).
This project was created for the [Tewari Lab](http://www.tewarilab.org/) at the University of Michigan.

## Setup

You can download this repository with:
```
$ git clone https://github.com/kelly-sovacool/latex-templates.git
```

Python 3 packages to install:
 - pandas
 - docopt
 - numpy
 - yaml

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

The program will create donors metadata, biosamples metadata, and a manifest file and write them to `working_dir`.
