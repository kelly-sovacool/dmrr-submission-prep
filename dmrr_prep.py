#!/usr/local/bin/python3
""" Preparing metadata for submission to the DMRR
Author: Kelly Sovacool
Email: sovacool@umich.edu
Date: July 2018

Usage:
    dmrr_prep.py <config_filename>

Options:
    -h --help   print this help message

"""
import datetime
import docopt
import json
import numpy as np
import os
import pandas as pd
import pprint
import yaml
import time

MISSING = "#MISSING#"


class Submission:
    def __init__(self, group, user_login, study_name, samples_filename, sample_study_names, study_id, working_dir, templates_dir, md5sum, fastq_filenames, is_time_series, database):
        self.group = group
        self.user_login = user_login
        self.study_name = study_name
        self.study_id = study_id
        self.templates_dir = templates_dir
        self.md5sum = md5sum
        self.is_time_series = is_time_series
        self.database = database

        # load all the data
        self.samples = Samples(samples_filename, sample_study_names, is_time_series)
        self.participants = Participants(samples_filename, sample_study_names, is_time_series)
        self.donors = Donors(templates_dir, self.participants.df, study_id, is_time_series)
        self.biosamples = Biosamples(study_id, templates_dir, working_dir, self.participants.df, self.samples.df, is_time_series)
        self.manifest = Manifest(self.biosamples, working_dir, templates_dir, fastq_filenames, self.samples, study_id, study_name, user_login, md5sum, group, database)

        # Save all the files
        self.donors.write(os.path.join(working_dir, self.manifest['donorMetadataFileName']))
        self.biosamples.write()
        self.manifest.write()

    @classmethod
    def from_config_dict(cls, config):
        return cls(config['group'], config['user_login'], config['study_name'], config['samples_filename'], config['sample_study_names'], config['study_id'], config['working_dir'], config['templates_dir'], config['md5sum'], config['fastq_filenames'], config['is_time_series'], config['database'])


class Manifest(dict):
    def __init__(self, biosamples, working_dir, templates_dir, fastq_filenames, samples, study_id, study_name, user_login, md5sum, group, database):
        with open(os.path.join(templates_dir, 'manifest_template.manifest.json'), 'r') as template_file:
            super().__init__(json.load(template_file))
        self.study_id = study_id
        self.working_dir = working_dir

        # Fill in the manifest metadata
        # TODO: option to supply list of names, or actual tar archive
        # TODO: if tar archive exists, compute md5sum (instead of config param) and traverse to get filenames
        self['settings']['analysisName'] = 'MTEWA1_{}_{}'.format(study_name.strip('EXR-MTEWA1'), datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d'))
        self['studyName'] = study_name
        self['userLogin'] = user_login
        self['md5CheckSum'] = md5sum
        self['group'] = group
        self['db'] = database
        self['runMetadataFileName'] = study_id + '-RU.metadata.tsv'
        self['submissionMetadataFileName'] = study_id + '-SU.metadata.tsv'
        self['studyMetadataFileName'] = study_id + '-ST.metadata.tsv'
        if len(biosamples.dfs) == 1:
            self['experimentMetadataFileName'] = study_id + '-EX.metadata.tsv'
            self['biosampleMetadataFileName'] = study_id + '-BS.metadata.tsv'
        else:
            self.pop('experimentMetadataFileName')
            self.pop("biosampleMetadataFileName")
        self['donorMetadataFileName'] = study_id + '-DO.metadata.tsv'
        self['manifest'] = list()

        # Fill in the samples to filenames map
        with open(fastq_filenames, 'r') as file:
            sample_filenames = {int(f.split('_')[0]): f.strip() for f in file}  # { MT.Unique.ID: fastq_filename }
        for mt_unique_id in sorted(samples.df.index):  # only include fastq filenames in this study
            fastq_filename = sample_filenames[mt_unique_id]
            sample_name = 'MT.Unique.ID_' + str(mt_unique_id)
            sample_source = samples.df.loc[mt_unique_id, 'Source']
            if len(biosamples.dfs) == 1:
                sample_dict = {'sampleName': sample_name, 'dataFileName': fastq_filename}
            else:
                sample_dict = {'sampleName': sample_name, 'dataFileName': fastq_filename, 'biosampleMetadataFileName': study_id + sample_source + '-BS.metadata.tsv', 'experimentMetadataFileName': study_id + sample_source + '-EX.metadata.tsv'}
            self['manifest'].append(sample_dict)

    def write(self):
        with open(os.path.join(self.working_dir, self.study_id + '.manifest.json'), 'w') as file:
            json.dump(self, file, indent=4, sort_keys=True)

    @property
    def filenames(self):
        automatic_filetypes = {'biosampleMetadataFileName', 'donorMetadataFileName'}
        filenames = set()
        for key, value in self.items():
            if 'filename' in key.lower() and 'metadata' in value.lower() and key not in automatic_filetypes:
                filenames.add(value)
            elif type(value) == list:
                for thing in value:
                    if type(thing) == dict:
                        for key2, value2 in thing.items():
                            if 'filename' in key2.lower() and 'metadata' in value2.lower() and key2 not in automatic_filetypes:
                                filenames.add(value2)
        return sorted(list(filenames))


class MetaDataFrame:
    def __init__(self, csv_filename=None, dataframe=None):
        if csv_filename:
            self.df = pd.read_csv(csv_filename, sep='\t' if csv_filename.endswith('.tsv') else ',')
        elif type(dataframe) == pd.DataFrame:
            self.df = dataframe.copy()
        else:
            self.df = pd.DataFrame()

    def write(self, output_filename):
        self.df.to_csv(output_filename, sep='\t', index=True)


class Samples(MetaDataFrame):
    def __init__(self, samples_filename, sample_study_names, is_time_series):
        super().__init__(csv_filename=samples_filename)
        # Load the samples spreadsheet
        # only use healthy control study and samples that passed quality control
        self.df = self.df.loc[(self.df['Study'].isin(sample_study_names)) & (self.df['MISEQ.QC.PASS'] == 'PASS')]
        # keep only the necessary columns
        self.df = self.df[['Participant.ID', 'Sample.ID', 'MT.Unique.ID', 'Source', 'Study']]
        self.df = self.df.set_index(['MT.Unique.ID']).sort_values(by='Participant.ID')
        if is_time_series:
            self.df['timepoint'] = pd.Series([sample_id.split('-')[1] if len(sample_id.split('-')) > 1 else 'T0' for sample_id in self.df['Sample.ID']], index=self.df.index)
            self.df['timestep'] = pd.Series([study.split('-')[1] if len(study.split('-')) > 1 else 'NA' for study in self.df['Study']], index=self.df.index)


class Participants(MetaDataFrame):
    def __init__(self, samples_filename, sample_study_names, is_time_series):
        super().__init__(csv_filename=samples_filename)
        # only use healthy control study and samples that passed quality control
        self.df = self.df.loc[(self.df['Study'].isin(sample_study_names)) & (self.df['MISEQ.QC.PASS'] == 'PASS')]
        self.df["Source"] = self.df["Source"].str.capitalize()
        # remove duplicate participants -- all have plasma but some additionally have serum  TODO: make this more general to work with other study designs
        self.df = self.df.loc[(self.df['Source'] == 'Plasma')]
        if is_time_series:
            self.df['timepoint'] = pd.Series([sample_id.split('-')[1] if len(sample_id.split('-')) > 1 else 'T0' for sample_id in self.df['Sample.ID']], index=self.df.index)
            self.df['timestep'] = pd.Series([study.split('-')[1] if len(study.split('-')) > 1 else 'NA' for study in self.df['Study']], index=self.df.index)
            self.df = self.df.loc[(self.df['timepoint'] == 'T0')]
        self.df = self.df[['Participant.ID', 'Age', 'Race', 'Gender']].set_index('Participant.ID')

        # Use correct ontology terms
        # TODO: make this more generic. change given terms to lowercase and check if individual words (split by spaces and other symbols) match any keys.
        race_ontology = {'Asian': 'Asian', 'asian': 'Asian',
                         'Black or African American': 'African American',
                         'mixed/asian & white': 'Multiracial',
                         'mixed/Asian &Black': 'Multiracial',
                         'mixed/black, white, asian': 'Multiracial',
                         'Native Hawiian or other Pacific Islander': 'Native Hawaiian or Other Pacific Islander',
                         'Pacific Islander': 'Native Hawaiian or Other Pacific Islander',
                         'White': 'White', 'white': 'White',
                         np.nan: MISSING}
        for index in self.df.index:
            race = self.df.at[index, 'Race']
            self.df.loc[index, 'Race'] = race_ontology[race] if race in race_ontology else 'Multiracial'
            age = self.df.at[index, 'Age']
            if not (0 <= age <= 130):
                self.df.loc[index, 'Age'] = MISSING
            gender = self.df.loc[index, 'Gender']
            if type(gender) == float:
                self.df.loc[index, 'Gender'] = MISSING
            elif gender != MISSING:
                self.df.loc[index, 'Gender'] = gender.capitalize()


class Donors(MetaDataFrame):
    def __init__(self, templates_dir, participants, study_id, is_time_series):
        # Load the donors template
        super().__init__(csv_filename=os.path.join(templates_dir, 'Donors.template.tsv'))
        self.df = self.df.set_index('#property')
        self.df.drop(['- Ethnic Group', '-- Current Health Status', '-- Medical History', '-- Smoking History', '-- Medications', '-- Treatment History',
                      '-- Family History', '-- Treatment History', '-- Family History', '-- Developmental Stage', '- Has Expired?', '-- Estimated Date',
                      '-- Post-mortem Interval', '- Notes', '* Family Members', '*- Family Member', '*-- Relationship', '*-- DocURL', '* Aliases', '*-  Accession',
                      '*-- dbName', '*-- URL', '- Health Status', '*-- Notes'], inplace=True)
        genders = {'Male', 'Female'}
        # Fill in the donors dataframe
        for index, part_id in enumerate(participants.index):
            participant_column = 'value' + part_id
            donor_id = study_id + str(index + 1) + '-DO'
            participants.loc[part_id, 'donor.id'] = donor_id
            self.df.insert(index + 1, participant_column, self.df['value'])
            self.df.loc['Donor', participant_column] = donor_id  # for matching biosamples to donor ids
            self.df.loc['- Status', participant_column] = 'Protect' if is_time_series else 'Add'
            gender = participants.at[part_id, 'Gender']
            self.df.loc['- Sex', participant_column] = gender if gender in genders else '#MISSING#'
            self.df.loc['- Racial Category', participant_column] = participants.at[part_id, 'Race']
            self.df.loc['- Donor Type', participant_column] = 'Healthy Subject' if not is_time_series else 'Experimental'
            age = str(participants.at[part_id, 'Age'])
            self.df.loc['- Age', participant_column] = age + ' years' if age != '#MISSING#' else age
            self.df.loc['* Custom Metadata', participant_column] = 1
            self.df.loc['*- Property Name', participant_column] = 'Participant.ID'
            self.df.loc['*-- Value', participant_column] = part_id
        self.df = self.df.drop('value', axis=1)
        self.df.columns = [''.join(l for l in col if not l.isdigit() and l != '*') for col in list(self.df.columns)]
        self.df.index = pd.Index([''.join(l for l in col if not l.isdigit()) for col in list(self.df.index)], name=self.df.index.name)


class Biosamples:
    def __init__(self, study_id, templates_dir, working_dir, participants, samples, is_time_series):
        self.filename_base = os.path.join(working_dir, study_id)
        self.filename_ext = '-BS.metadata.tsv'

        timepoints = {"1 hr": {"T0": 'Fasting blood draw',
                               "T1": "1 hr post meal 1",
                               "T2": "2 hrs post meal 1",
                               "T3": "3 hrs post meal 1",
                               "T4": "4 hrs post meal 1",
                               "T5": "1 hr post meal 2",
                               "T6": "2 hrs post meal 2",
                               "T7": "3 hrs post meal 2",
                               "T8": "4 hrs post meal 2",
                               "T9": "24 hrs post meal 2",
                               "T10": "48 hrs post meal 2"},
                      "30 min": {"T0": 'Fasting blood draw',
                                 "T1": "0.5 hr post meal 1",
                                 "T2": "1 hr post meal 1",
                                 "T3": "1.5 hrs post meal 1",
                                 "T4": "2 hrs post meal 1",
                                 "T5": "2.5 hrs post meal 1",
                                 "T6": "3 hrs post meal 1",
                                 "T7": "3.5 hrs post meal 1",
                                 "T8": "4 hrs post meal 1",
                                 "T9": "24 hrs post meal 2",
                                 "T10": "48 hrs post meal 2"}}

        # load the template and clean it up
        template = MetaDataFrame(csv_filename=os.path.join(templates_dir, 'Biosamples.template.tsv'))
        template.df = template.df.set_index('#property')
        template.df = template.df.drop(['-- Age at Sampling', '-- Notes', '- Description', '--- Symptoms', '--- Pathology', '--- Disease Duration',
                                        '--- Collection Details',
                                        '---- Sample Collection Method', '---- Geographic Location',
                                        '---- Collection Date', '---- Time of Collection',
                                        '---- Collection Tube Type', '----- Other Collection Tube Type',
                                        '---- Holding Time', '---- Holding Temperature',
                                        '---- Preservatives Used', '---- Freezing Method',
                                        '---- Number of Times Freeze Thawed',
                                        '---- Contamination Removal Method', '--- Notes',
                                        '-- Cell Culture Supernatant', '--- Source', '---- Type',
                                        '---- Cell Line', '---- Start Date', '---- Harvest Date', '--- Tissue',
                                        '---- Date Obtained', '---- Tissue Type', '--- Notes',
                                        '-- Starting Amount', '-- Replicate Information',
                                        '--- Biological Replicate Number', '--- Technical Replicate Number', '-- Provider', '--- Company Name', '--- Lab Name', '--- Person Name',
                                        '* Pooled Biosamples', '*- Pooled Biosample', '*-- DocURL', '* Aliases', '*-  Accession', '*-- dbName', '*-- URL',
                                        '*-- Date Submitted to External Database', '*-- Notes'])
        template.df = template.df.T
        template.df.insert(18, '*-- DocURL', [np.nan, 'URL', np.nan, np.nan, 'Relative ID (accession) of doc, provide Document URL'])
        if is_time_series:
            template.df.insert(23, '*- Property Name2', template.df['*- Property Name'])
            template.df.insert(24, '*-- Value2', template.df['*-- Value'])
            template.df.insert(25, '*- Property Name3', template.df['*- Property Name'])
            template.df.insert(26, '*-- Value3', template.df['*-- Value'])
            template.df.insert(27, '*- Property Name4', template.df['*- Property Name'])
            template.df.insert(28, '*-- Value4', template.df['*-- Value'])
        template.df = template.df.T

        sources = set(samples['Source'])
        self.dfs = dict()  # one dataframe per source
        for source in sources:
            self.dfs[source] = MetaDataFrame(dataframe=template.df)
        start_column_length = len(template.df.columns)
        for mt_unique_id in sorted(samples.index):  # Fill in the biosamples dataframes
            donor_id = participants.loc[samples.loc[mt_unique_id, 'Participant.ID'], 'donor.id']
            sample_column = 'value' + str(mt_unique_id)
            sample_source = samples.loc[mt_unique_id, 'Source']
            index = len(self.dfs[sample_source].df.columns) - start_column_length + 1
            self.dfs[sample_source].df.insert(index, sample_column, self.dfs[sample_source].df['value'])
            self.dfs[sample_source].df.loc['Biosample', sample_column] = study_id + str(index) + '-BS'
            self.dfs[sample_source].df.loc['- Status', sample_column] = 'Protect' if is_time_series else 'Add'
            self.dfs[sample_source].df.loc['- Name', sample_column] = 'MT.Unique.ID_' + str(mt_unique_id)
            self.dfs[sample_source].df.loc['- Donor ID', sample_column] = donor_id
            self.dfs[sample_source].df.loc['-- DocURL', sample_column] = 'coll/Donors/doc/' + donor_id
            self.dfs[sample_source].df.loc['--- Scientific Name', sample_column] = 'Homo sapiens'
            self.dfs[sample_source].df.loc['--- Common Name', sample_column] = 'Human'
            self.dfs[sample_source].df.loc['--- Taxon ID', sample_column] = 9606
            self.dfs[sample_source].df.loc['-- Disease Type', sample_column] = 'Healthy Subject'
            self.dfs[sample_source].df.loc['-- Anatomical Location', sample_column] = 'Plasma cell'
            self.dfs[sample_source].df.loc['--- Biofluid Name', sample_column] = sample_source
            self.dfs[sample_source].df.loc['-- exRNA Source', sample_column] = ' total cell-free biofluid RNA'
            self.dfs[sample_source].df.loc['-- Fractionation', sample_column] = 'Yes'
            self.dfs[sample_source].df.loc['* Related Experiments', sample_column] = 1
            self.dfs[sample_source].df.loc['*- Related Experiment', sample_column] = study_id + '1-EX'
            self.dfs[sample_source].df.loc['*-- DocURL', sample_column] = 'coll/Experiments/doc/' + study_id + '1-EX'
            self.dfs[sample_source].df.loc['* Custom Metadata', sample_column] = 4 if is_time_series else 1
            self.dfs[sample_source].df.loc['*- Property Name', sample_column] = 'Participant.ID'
            self.dfs[sample_source].df.loc['*-- Value', sample_column] = samples.loc[mt_unique_id, 'Participant.ID']
            if is_time_series:
                self.dfs[sample_source].df.loc['*- Property Name2', sample_column] = 'timepoint_id'
                timepoint_id = samples.loc[mt_unique_id, 'timepoint']
                self.dfs[sample_source].df.loc['*-- Value2', sample_column] = timepoint_id
                self.dfs[sample_source].df.loc['*- Property Name3', sample_column] = 'timepoint_description'
                time_between_collections = samples.loc[mt_unique_id, 'timestep']
                self.dfs[sample_source].df.loc['*-- Value3', sample_column] = timepoints[time_between_collections][timepoint_id]
                self.dfs[sample_source].df.loc['*- Property Name4', sample_column] = 'time_between_sample_collections_0-8'
                self.dfs[sample_source].df.loc['*-- Value4', sample_column] = time_between_collections

        for key in self.dfs:  # drop the empty "value" column and rename all numbered value columns to just "value"
            self.dfs[key].df = self.dfs[key].df.drop('value', axis=1)
            self.dfs[key].df.columns = [''.join(l for l in col if not l.isdigit() and l != '*') for col in list(self.dfs[key].df.columns)]
            self.dfs[key].df.index = pd.Index([''.join(l for l in col if not l.isdigit()) for col in list(self.dfs[key].df.index)], name=self.dfs[key].df.index.name)

        if len(self.dfs) == 1:  # reassign the lone dataframe to an empty string dict key
            self.dfs[''] = self.dfs.pop(set(self.dfs.keys()).pop())

    def write(self):
        for source, df in self.dfs.items():
            df.write(self.filename_base + source + self.filename_ext)


def main(args):
    start = time.time()
    with open(args['<config_filename>'], 'r') as file:
        config = yaml.load(file)
    print('Configuration:')
    pprint.pprint(config)
    if not os.path.exists(config['working_dir']):
        os.mkdir(config['working_dir'])
    submission = Submission.from_config_dict(config)
    print('Prepared donors, biosamples, and manifest files for {} in {}'.format(config['study_id'], datetime.timedelta(seconds=time.time() - start)))
    print("Don't forget to create and fill in the following files:")
    for filename in submission.manifest.filenames:
        print('\t', filename)
    # Don't forget to manually fill in experiment, run, study, and submission metadata files!


if __name__ == "__main__":
    arguments = docopt.docopt(__doc__)
    main(arguments)
