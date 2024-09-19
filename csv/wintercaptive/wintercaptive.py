import os, os.path

# CONFIG
TARGET_DIR = '/home/cbaldoni/data/cue_learning/csv/wintercaptive' # Directory containing csv files to be combined.
STRING_OF_COLUMN_LABELS = 'season, ID, trial, time, frame, x, y, ANGLE, speed, acceleration\n'
DESTINATION_FILE_NAME = 'wintercaptive.csv' #change name of .csv file
VERBOSE = True

# CODE

DESTINATION_FILE_NAME = os.path.join(TARGET_DIR, DESTINATION_FILE_NAME)

with open(DESTINATION_FILE_NAME, 'w') as dest:
    if VERBOSE:
        print('Added header')
    dest.write(STRING_OF_COLUMN_LABELS)

LIST_OF_TARGETS = os.listdir(TARGET_DIR)
LIST_OF_TARGETS.sort()

for File in LIST_OF_TARGETS:
    if File == os.path.basename(DESTINATION_FILE_NAME) or File[-4:] != '.csv':
        continue
    if VERBOSE:
        print('Now working on', File)
    lines = [line for line in open(os.path.join(TARGET_DIR, File))][1:]
    with open(DESTINATION_FILE_NAME, 'a') as dest:
        for line in lines:
            dest.write(line)
