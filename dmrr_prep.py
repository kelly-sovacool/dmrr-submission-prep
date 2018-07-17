#!/usr/local/bin/python3
""" Preparing metadata for submission to the DMRR

Usage:
    dmrr_prep.py <config_filename>

Options:
    -h --help   print this help message

"""
import abc
import datetime
import docopt
import json
import numpy as np
import os
import pandas as pd
import pprint
import yaml
import time

# TODO: handle multiple sources -- e.g. separate experiments & biosamples file for plasma & serum, or single file when there's only one source.
# TODO: time-series custom metadata


class Submission:
    def __init__(self, group, user_login, study_name, samples_filename, sample_study_names, study_id, working_dir, templates_dir, md5sum, is_time_series, database):
        self.group = group
        self.user_login = user_login
        self.study_name = study_name
        self.study_id = study_id
        self.templates_dir = templates_dir
        self.md5sum = md5sum
        self.is_time_series = is_time_series
        self.database = database

        self.samples = Samples(samples_filename, sample_study_names, is_time_series)
        self.participants = Participants(samples_filename, sample_study_names, is_time_series)
        self.donors = Donors(templates_dir, self.participants.df, study_id, is_time_series)
        self.biosamples = Biosamples(study_id, templates_dir, os.path.join(working_dir, study_id), self.participants.df, self.samples.df, is_time_series)

        with open(os.path.join(templates_dir, 'manifest_template.manifest.json'), 'r') as file:
            self.manifest = json.load(file)
        self.load_manifest()

        # Strip numbers from the value columns and property indices
        for df in [self.donors, self.biosamples]:
            df.columns = [''.join(l for l in col if not l.isdigit() and l != '*') for col in list(df.columns)]
            df.index = pd.Index([''.join(l for l in col if not l.isdigit()) for col in list(df.index)], name=df.index.name)

        # Save all the files
        self.donors.write(os.path.join(working_dir, self.manifest['donorMetadataFileName']))
        self.biosamples.write()
        with open(working_dir + study_id + '.manifest.json', 'w') as file:
            json.dump(self.manifest, file, indent=4)

    def load_manifest(self):
        # Load the list of fastq filenames
        with open('HealthyControl_and_FeedingStudy_fastq_file.names.txt', 'r') as file:
            sample_filenames = {int(f.split('_')[0]): f.strip() for f in file}  # { MT.Unique.ID: fastq_filename }

        # Fill in the manifest
        # TODO: option to supply list of names, or actual tar archive
        # TODO: if tar archive exists, compute md5sum (instead of config param) and traverse to get filenames
        self.manifest['settings']['analysisName'] = 'MTEWA1_Healthy_Controls_' + datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d')
        self.manifest['studyName'] = self.study_name
        self.manifest['userLogin'] = self.user_login
        self.manifest['md5CheckSum'] = self.md5sum
        self.manifest['group'] = self.group
        self.manifest['db'] = self.database
        self.manifest['runMetadataFileName'] = self.study_id + '-RU.metadata.tsv'
        self.manifest['submissionMetadataFileName'] = self.study_id + '-SU.metadata.tsv'
        self.manifest['studyMetadataFileName'] = self.study_id + '-ST.metadata.tsv'
        self.manifest['experimentMetadataFileName'] = self.study_id + '-EX.metadata.tsv'
        self.manifest['biosampleMetadataFileName'] = self.study_id + '-BS.metadata.tsv'
        self.manifest['donorMetadataFileName'] = self.study_id + '-DO.metadata.tsv'
        self.manifest['manifest'] = list()
        for mt_unique_id in sorted(self.samples.index):  # only include fastq filenames in this study
            fastq_filename = sample_filenames[mt_unique_id]
            sample_name = 'MT.Unique.ID_' + str(mt_unique_id)
            self.manifest['manifest'].append({'sampleName': sample_name, 'dataFileName': fastq_filename})
        len(self.manifest['manifest'])


class MetaDataFrame:

    def __init__(self, csv_filename=None, other_df=None):
        if csv_filename:
            self.df = pd.read_csv(csv_filename, sep='\t' if csv_filename.endswith('.tsv') else ',')
        elif other_df:
            self.df = other_df.copy()
        else:
            self.df = pd.DataFrame()

    def write(self, output_filename):
        self.df.to_csv(output_filename, sep='\t', index=True)


class Samples(MetaDataFrame):
    def __init__(self, samples_filename, sample_study_names, is_time_series):
        super().__init__(csv_filename=samples_filename)
        # Load the samples spreadsheet
        # only use healthy control study and samples that passed quality control
        self.samples = self.samples.loc[(self.samples['Study'].isin(sample_study_names)) & (self.samples['MISEQ.QC.PASS'] == 'PASS')]
        # keep only the necessary columns
        self.samples = self.samples[['Participant.ID', 'Sample.ID', 'MT.Unique.ID', 'Source', 'Study']]
        self.samples = self.samples.set_index(['MT.Unique.ID']).sort_values(by='Participant.ID')
        if is_time_series:
            self.samples['timepoint'] = pd.Series([sample_id.split('-')[1] for sample_id in self.samples['Sample.ID']], index=self.samples.index)
            self.samples['timestep'] = pd.Series([study.split('-')[1] for study in self.samples['Study']], index=self.samples.index)


class Participants(MetaDataFrame):
    def __init__(self, samples_filename, sample_study_names, is_time_series):
        super().__init__(csv_filename=samples_filename)
        # only use healthy control study and samples that passed quality control
        self.df = self.df.loc[(self.df['Study'].isin(sample_study_names)) & (self.df['MISEQ.QC.PASS'] == 'PASS')]
        for column in ['Gender', 'Race', 'Source']:  # strip whitespace for consistency
            self.df[column] = self.df[column].str.strip()
            if column == 'Gender':  # capitalize first letter
                self.df[column] = self.df[column].str.capitalize()
        # remove duplicate participants -- all have plasma but some additionally have serum  TODO: make this more general to work with other study designs
        self.df = self.df.loc[(self.df['Source'] == 'Plasma')]
        if is_time_series:
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
                         '#MISSING#': '#MISSING#'}
        for part_id in self.df.index:
            race = self.df.at[part_id, 'Race']
            self.df.loc[part_id, 'Race'] = race_ontology[race] if race in race_ontology else 'Multiracial'
        pprint.pprint({'Age': set(self.df['Age']), 'Race': set(self.df['Race']), 'Gender': set(self.df['Gender'])})


class Donors(MetaDataFrame):
    def __init__(self, templates_dir, participants, study_id, is_time_series):
        # Load the donors template
        super().__init__(csv_filename=os.path.join(templates_dir, 'Donors.template.tsv'))
        self.df = self.df.set_index('#property')
        self.df.drop(['- Ethnic Group', '-- Current Health Status', '-- Medical History', '-- Smoking History', '-- Medications', '-- Treatment History',
                     '-- Family History', '-- Treatment History', '-- Family History', '-- Developmental Stage', '- Has Expired?', '-- Estimated Date',
                     '-- Post-mortem Interval', '- Notes', '* Family Members', '*- Family Member', '*-- Relationship', '*-- DocURL', '* Aliases', '*-  Accession',
                     '*-- dbName', '*-- URL', '- Health Status', '*-- Notes'], inplace=True)
        # Fill in the donors dataframe
        i = 1
        for part_id in participants.index:
            print(part_id)
            participant_column = 'value' + part_id
            donor_id = study_id + str(i) + '-DO'
            participants.loc[part_id, 'donor.id'] = donor_id
            self.df.insert(i, participant_column, self.df['value'])
            self.df.loc['Donor', participant_column] = donor_id  # for matching biosamples to donor ids
            self.df.loc['- Status', participant_column] = 'Add'
            self.df.loc['- Sex', participant_column] = participants.at[part_id, 'Gender']
            self.df.loc['- Racial Category', participant_column] = participants.at[part_id, 'Race']
            self.df.loc['- Donor Type', participant_column] = 'Healthy Subject' if not is_time_series else 'Experimental'
            self.df.loc['- Age', participant_column] = str(participants.at[part_id, 'Age']) + ' years'
            self.df.loc['* Custom Metadata', participant_column] = 1
            self.df.loc['*- Property Name', participant_column] = 'Participant.ID'
            self.df.loc['*-- Value', participant_column] = part_id
            i += 1
        self.df = self.df.drop('value', axis=1)
        self.df = self.df


class Biosamples:
    def __init__(self, study_id, templates_dir, filename_base, participants, samples, is_time_series):
        self.master_biosamples = MetaDataFrame(csv_filename=os.path.join(templates_dir, 'Biosamples.template.tsv'))
        self.filename_base = filename_base
        self.filename_ext = '-EX.metadata.tsv'

        # TODO: need multiple biosamples templates for multiple sources (plasma, serum)

        self.master_biosamples = self.master_biosamples.set_index('#property')
        self.master_biosamples = self.master_biosamples.drop(['-- Age at Sampling', '-- Notes', '- Description', '--- Symptoms', '--- Pathology', '--- Disease Duration',
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
        self.master_biosamples = self.master_biosamples.T
        self.master_biosamples.insert(18, '*-- DocURL', [np.nan, 'URL', np.nan, np.nan, 'Relative ID (accession) of doc, provide Document URL'])
        if is_time_series:
            self.master_biosamples.insert(23, '*- Property Name2', self.master_biosamples['*- Property Name'])
            self.master_biosamples.insert(24, '*-- Value2', self.master_biosamples['*-- Value'])
            self.master_biosamples.insert(25, '*- Property Name3', self.master_biosamples['*- Property Name'])
            self.master_biosamples.insert(26, '*-- Value3', self.master_biosamples['*-- Value'])
        self.master_biosamples = self.master_biosamples.T

        # Fill in the biosamples dataframe
        i = 1
        for mt_unique_id in samples.index:
            donor_id = participants.loc[samples.loc[mt_unique_id, 'Participant.ID'], 'donor.id']
            sample_column = 'value' + str(mt_unique_id)
            self.master_biosamples.insert(i, sample_column, self.master_biosamples['value'])
            self.master_biosamples.loc['Biosample', sample_column] = study_id + str(i) + '-BS'
            self.master_biosamples.loc['- Status', sample_column] = 'Add'
            self.master_biosamples.loc['- Name', sample_column] = 'MT.Unique.ID_' + str(mt_unique_id)
            self.master_biosamples.loc['- Donor ID', sample_column] = donor_id
            self.master_biosamples.loc['-- DocURL', sample_column] = 'coll/Donors/doc/' + donor_id
            self.master_biosamples.loc['--- Scientific Name', sample_column] = 'Homo sapiens'
            self.master_biosamples.loc['--- Common Name', sample_column] = 'Human'
            self.master_biosamples.loc['--- Taxon ID', sample_column] = 9606
            self.master_biosamples.loc['-- Disease Type', sample_column] = 'Healthy Subject'
            self.master_biosamples.loc['-- Anatomical Location', sample_column] = 'Plasma cell'
            self.master_biosamples.loc['--- Biofluid Name', sample_column] = samples.loc[mt_unique_id, 'Source']
            self.master_biosamples.loc['-- exRNA Source', sample_column] = ' total cell-free biofluid RNA'
            self.master_biosamples.loc['-- Fractionation', sample_column] = 'Yes'
            self.master_biosamples.loc['* Related Experiments', sample_column] = 1
            self.master_biosamples.loc['*- Related Experiment', sample_column] = study_id + '1-EX'
            self.master_biosamples.loc['*-- DocURL', sample_column] = 'coll/Experiments/doc/' + study_id + '1-EX'
            self.master_biosamples.loc['* Custom Metadata', sample_column] = 2 if is_time_series else 1
            self.master_biosamples.loc['*- Property Name', sample_column] = 'Participant.ID'
            self.master_biosamples.loc['*-- Value', sample_column] = samples.loc[mt_unique_id, 'Participant.ID']
            if is_time_series:
                self.master_biosamples.loc['*- Property Name2', sample_column] = 'timepoint'
                self.master_biosamples.loc['*-- Value2', sample_column] = samples.loc[mt_unique_id, 'timepoint']
                self.master_biosamples.loc['*- Property Name3', sample_column] = 'timestep'
                self.master_biosamples.loc['*-- Value3', sample_column] = samples.loc[mt_unique_id, 'timestep']
            i += 1
        self.master_biosamples = self.master_biosamples.drop('value', axis=1)

        # TODO: split master dataframe into one per source
        self.dfs = {'Plasma': MetaDataFrame(), 'Serum': MetaDataFrame()}  # TODO: multiple biosamples files if multiple sources
        # self.dfs = {'': MetaDataFrame()}  # TODO: single biosamples file for one source

    def write(self):
        for source, df in self.dfs.values():
            df.write(self.filename_base + source + self.filename_ext)


def main(args):
    with open(args['<config_filename>'], 'r') as file:
        config = yaml.load(file)
    pprint.pprint(config)
    templates_dir = 'templates'
    submission = Submission(config['group'], config['user_login'], config['study_name'], config['samples_filename'], config['sample_study_names'], config['study_id'], config['dir'], templates_dir, config['md5sum'], config['is_time_series'], config['database'])
    # Don't forget to manually fill in experiment, run, study, and submission metadata files!


if __name__ == "__main__":
    arguments = docopt.docopt(__doc__)
    main(arguments)
