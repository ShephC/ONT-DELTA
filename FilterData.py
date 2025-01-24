#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 08:55:30 2025

@author: shephc3
"""

import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import sys

#Here we'll define the directories and files for analysis.
#Change according to the dataset
dataDir = '/Users/shephc3/Desktop/untitled folder/'
stressFileRight = dataDir + 'rightForks_DNAscent_forkSense_stressSignatures.bed'
stressFileLeft = dataDir + 'leftForks_DNAscent_forkSense_stressSignatures.bed'
origins = dataDir + 'origins_DNAscent_forkSense.bed'
terminations = dataDir + 'terminations_DNAscent_forkSense.bed'
outputDir = '/Users/shephc3/Desktop/untitled folder/'


#This will process a stressfile and return the useful columns for tracks that are not cutoff
#This will also keep track of the number of cutoff tracks for calculating the percent cutoff
###THIS FILTERS OUT CUTOFF TRACKS###
def processStressFile(file):
   global cutoffCount
   f = open(file, 'r')
   chromosome = []
   trackStart = []
   trackEnd = []
   readID = []
   readStart = []
   readEnd = []
   direction = []
   firstAnalogueLength = []
   secondAnalogueLength = []
   detectionLength = []
   for line in f:
      splitLine = line.rstrip().split()
      if line[0] == '#':
         continue
         #Ignore header lines
      elif splitLine[14] != '-3.000000':
         #Only collect reads that aren't identified as cutoff (-3) in the last column
         chromosome.append(splitLine[0])
         trackStart.append(splitLine[1])
         trackEnd.append(splitLine[2])
         readID.append(splitLine[3])
         readStart.append(splitLine[4])
         readEnd.append(splitLine[5])
         direction.append(splitLine[6])
         detectionLength.append(float(splitLine[7]))
         firstAnalogueLength.append(float(splitLine[8]))
         secondAnalogueLength.append(float(splitLine[9]))
   filterDF = pd.DataFrame(
      {"Chromosome" : chromosome,
       "Track Start" : trackStart,
       "Track End" : trackEnd,
       "Read ID" : readID,
       "Read Start" : readStart,
       "Read End" : readEnd,
       "Direction" : direction,
       "First Analogue Length" : firstAnalogueLength,
       "Second Analogue Length" : secondAnalogueLength,
       "Detect Length" : detectionLength}, index=readID)
##Add in 'keys' at a later date to differentiate within new DF
   f.close()
   return filterDF

#This function returns an array of the read IDs in an origin or termination file.
#This is used for filtering the origin and termination sites from the forksense data 
#						to only give lengths for unimpeded forks
def processOriTerm(file):
   f = open(file, 'r')
   read = []
   for line in f:
      splitLine = line.rstrip().split()
      if line[0] != '#':
          read.append(splitLine[3])
   return read

rightDF = processStressFile(stressFileRight)
leftDF = processStressFile(stressFileLeft)
rightDF.to_csv(f'{outputDir}rightForksFiltered.bed', sep='\t', header=False, index=False, mode='a')
leftDF.to_csv(f'{outputDir}leftForksFiltered.bed', sep='\t', header=False, index=False, mode='a')

#This function is for use in creating bedgraph files that do not include cutoff reads
#It takes arrays of the required data and outputs it to a bedgraph file of the specified analogue.
def makeDF(analogue, chromosome, start, stop, probability, read):
	global outputDir
	dataframe = pd.DataFrame(
		{"Chromosome" : chromosome,
		 "Start" : start,
		 "Stop" : stop,
		 "Probability" :probability})
	with open(f"{outputDir}filtered{analogue}.bedgraph", 'a') as fd:
		fd.write(f"""track type=bedGraph name="{read}" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20 viewLimits=0.0:1.0\n""")
	dataframe.to_csv(f"{outputDir}filtered{analogue}.bedgraph", sep='\t', header=False, index=False, mode='a')

