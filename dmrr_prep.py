#!/usr/local/bin/python3
""" Preparing metadata for submission to the DMRR

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

# TODO: redesign w/ OOP for reusability
# TODO: handle multiple sources -- e.g. separate experiments & biosamples file for plasma & serum, or single file when there's only one source.
# TODO: time-series custom metadata
# TODO: use YAML config file
# TODO: call snakemake from within code for file generation steps


def main(args):

    with open(args['<config_filename>'], 'r') as file:
        config = yaml.load(file)
    pprint.pprint(config)

    study_id = config['study_id']
    study_names = config['study_names']
    directory = config['dir']
    md5sum = config['md5sum']
    is_time_series = config['is_time_series']

    # Load the samples spreadsheet

    samples_df = pd.read_csv('sample_sheet.csv')
    for column in ['Gender', 'Race', 'Source']:  # strip whitespace for consistency
        samples_df[column] = samples_df[column].str.strip()
        if column == 'Gender':  # capitalize first letter
            samples_df[column] = samples_df[column].str.capitalize()
    # only use healthy control study and samples that passed quality control
    samples_df = samples_df.loc[(samples_df['Study'].isin(study_names)) & (samples_df['MISEQ.QC.PASS'] == 'PASS')]
    # keep only the necessary columns
    samples_df = samples_df[['Participant.ID', 'Sample.ID', 'MT.Unique.ID', 'Age', 'Gender', 'Race', 'Source', 'Study']]
    samples_df = samples_df.set_index(['MT.Unique.ID']).sort_values(by='Participant.ID')
    if config['is_time_series']:
        samples_df['timepoint'] = pd.Series([sample_id.split('-')[1] for sample_id in samples_df['Sample.ID']], index=samples_df.index)
        samples_df['timestep'] = pd.Series([study.split('-')[1] for study in samples_df['Study']], index=samples_df.index)

    # Load the participant info

    # remove duplicate participants -- all have plasma but some additionally have serum
    participants = samples_df.loc[(samples_df['Source'] == 'Plasma')]
    if is_time_series:
        participants = participants.loc[(participants['timepoint'] == 'T0')]
    participants = participants[['Participant.ID', 'Age', 'Race', 'Gender']].set_index('Participant.ID')

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
    for part_id in participants.index:
        race = participants.at[part_id, 'Race']
        participants.loc[part_id, 'Race'] = race_ontology[race] if race in race_ontology else 'Multiracial'

    pprint.pprint({'Age': set(participants['Age']), 'Race': set(participants['Race']), 'Gender': set(participants['Gender'])})

    # Load the donors template

    donors = pd.read_csv('templates/Donors.template.tsv', sep='\t')
    donors = donors.set_index('#property')
    donors.drop(['- Ethnic Group', '-- Current Health Status', '-- Medical History', '-- Smoking History', '-- Medications', '-- Treatment History',
                 '-- Family History', '-- Treatment History', '-- Family History', '-- Developmental Stage', '- Has Expired?', '-- Estimated Date',
                 '-- Post-mortem Interval', '- Notes', '* Family Members', '*- Family Member', '*-- Relationship', '*-- DocURL', '* Aliases', '*-  Accession',
                 '*-- dbName', '*-- URL', '- Health Status', '*-- Notes'], inplace=True)
    donors.head()

    # Fill in the donors dataframe

    i = 1
    for part_id in participants.index:
        print(part_id)
        participant_column = 'value' + part_id
        donor_id = study_id + str(i) + '-DO'
        participants.loc[part_id, 'donor.id'] = donor_id
        donors.insert(i, participant_column, donors['value'])
        donors.loc['Donor', participant_column] = donor_id  # for matching biosamples to donor ids
        donors.loc['- Status', participant_column] = 'Add'
        donors.loc['- Sex', participant_column] = participants.at[part_id, 'Gender']
        donors.loc['- Racial Category', participant_column] = participants.at[part_id, 'Race']
        donors.loc['- Donor Type', participant_column] = 'Healthy Subject' if not is_time_series else 'Experimental'
        donors.loc['- Age', participant_column] = str(participants.at[part_id, 'Age']) + ' years'
        donors.loc['* Custom Metadata', participant_column] = 1
        donors.loc['*- Property Name', participant_column] = 'Participant.ID'
        donors.loc['*-- Value', participant_column] = part_id
        i += 1
    donors = donors.drop('value', axis=1)

    # Load the biosamples template

    biosamples = pd.read_csv('templates/Biosamples.template.tsv', sep='\t')
    biosamples = biosamples.set_index('#property')
    biosamples = biosamples.drop(['-- Age at Sampling', '-- Notes', '- Description', '--- Symptoms', '--- Pathology', '--- Disease Duration',
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

    biosamples = biosamples.T
    biosamples.insert(18, '*-- DocURL', [np.nan, 'URL', np.nan, np.nan, 'Relative ID (accession) of doc, provide Document URL'])
    if is_time_series:
        biosamples.insert(23, '*- Property Name2', biosamples['*- Property Name'])
        biosamples.insert(24, '*-- Value2', biosamples['*-- Value'])
        biosamples.insert(25, '*- Property Name3', biosamples['*- Property Name'])
        biosamples.insert(26, '*-- Value3', biosamples['*-- Value'])
    biosamples = biosamples.T

    # Fill in the biosamples dataframe

    i = 1
    for mt_unique_id in samples_df.index:
        donor_id = participants.loc[samples_df.loc[mt_unique_id, 'Participant.ID'], 'donor.id']
        sample_column = 'value' + str(mt_unique_id)
        biosamples.insert(i, sample_column, biosamples['value'])
        biosamples.loc['Biosample', sample_column] = study_id + str(i) + '-BS'
        biosamples.loc['- Status', sample_column] = 'Add'
        biosamples.loc['- Name', sample_column] = 'MT.Unique.ID_' + str(mt_unique_id)
        biosamples.loc['- Donor ID', sample_column] = donor_id
        biosamples.loc['-- DocURL', sample_column] = 'coll/Donors/doc/' + donor_id
        biosamples.loc['--- Scientific Name', sample_column] = 'Homo sapiens'
        biosamples.loc['--- Common Name', sample_column] = 'Human'
        biosamples.loc['--- Taxon ID', sample_column] = 9606
        biosamples.loc['-- Disease Type', sample_column] = 'Healthy Subject'
        biosamples.loc['-- Anatomical Location', sample_column] = 'Plasma cell'
        biosamples.loc['--- Biofluid Name', sample_column] = samples_df.loc[mt_unique_id, 'Source']
        biosamples.loc['-- exRNA Source', sample_column] = ' total cell-free biofluid RNA'
        biosamples.loc['-- Fractionation', sample_column] = 'Yes'
        biosamples.loc['* Related Experiments', sample_column] = 1
        biosamples.loc['*- Related Experiment', sample_column] = study_id + '1-EX'
        biosamples.loc['*-- DocURL', sample_column] = 'coll/Experiments/doc/' + study_id + '1-EX'
        biosamples.loc['* Custom Metadata', sample_column] = 2 if is_time_series else 1
        biosamples.loc['*- Property Name', sample_column] = 'Participant.ID'
        biosamples.loc['*-- Value', sample_column] = samples_df.loc[mt_unique_id, 'Participant.ID']
        if is_time_series:
            biosamples.loc['*- Property Name2', sample_column] = 'timepoint'
            biosamples.loc['*-- Value2', sample_column] = samples_df.loc[mt_unique_id, 'timepoint']
            biosamples.loc['*- Property Name3', sample_column] = 'timestep'
            biosamples.loc['*-- Value3', sample_column] = samples_df.loc[mt_unique_id, 'timestep']
        i += 1

    biosamples = biosamples.drop('value', axis=1)

    # Strip numbers from the value columns and property indices

    for df in [donors, biosamples]:
        df.columns = [''.join(l for l in col if not l.isdigit() and l != '*') for col in list(df.columns)]
        df.index = pd.Index([''.join(l for l in col if not l.isdigit()) for col in list(df.index)], name=df.index.name)

    # Load the manifest template file

    with open('templates/manifest_template.manifest.json', 'r') as file:
        manifest = json.load(file)

    # Load the list of fastq filenames

    with open('HealthyControl_and_FeedingStudy_fastq_file.names.txt', 'r') as file:
        sample_filenames = {int(f.split('_')[0]): f.strip() for f in file}  # { MT.Unique.ID: fastq_filename }

    # Fill in the manifest
    # TODO: option to supply list of names, or actual tar archive
    # TODO: if tar archive exists, compute md5sum (instead of config param) and traverse to get filenames

    manifest['settings']['analysisName'] = 'MTEWA1_Healthy_Controls_' + datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d')

    manifest['studyName'] = "U01 Healthy Controls July 2018"
    manifest['userLogin'] = 'sovacool'
    manifest['md5CheckSum'] = md5sum
    manifest['group'] = 'exrna-mtewa1'
    manifest['db'] = 'hg19_exRNA'
    manifest['runMetadataFileName'] = study_id + '-RU.metadata.tsv'
    manifest['submissionMetadataFileName'] = study_id + '-SU.metadata.tsv'
    manifest['studyMetadataFileName'] = study_id + '-ST.metadata.tsv'
    manifest['experimentMetadataFileName'] = study_id + '-EX.metadata.tsv'
    manifest['biosampleMetadataFileName'] = study_id + '-BS.metadata.tsv'
    manifest['donorMetadataFileName'] = study_id + '-DO.metadata.tsv'

    manifest['manifest'] = list()
    for mt_unique_id in sorted(samples_df.index):  # only include fastq filenames in this study
        fastq_filename = sample_filenames[mt_unique_id]
        sample_name = 'MT.Unique.ID_' + str(mt_unique_id)
        manifest['manifest'].append({'sampleName': sample_name, 'dataFileName': fastq_filename})
    len(manifest['manifest'])

    # Save all the files

    donors.to_csv(directory + manifest['donorMetadataFileName'], sep='\t', index=True)
    biosamples.to_csv(directory + manifest['biosampleMetadataFileName'], sep='\t')

    with open(directory + study_id + '.manifest.json', 'w') as file:
        json.dump(manifest, file, indent=4)

    # Validate that the metadata files listed in the manifest exist

    for key in manifest:
        if 'FileName' in key:
            assert os.path.isfile(directory + manifest[key])

    # Don't forget to manually fill in experiment, run, study, and submission metadata files!


if __name__ == "__main__":
    arguments = docopt.docopt(__doc__)
    main(arguments)
