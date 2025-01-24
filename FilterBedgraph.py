#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 16:38:20 2025

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
dataDir = '/Users/shephc3/Desktop/AEbedgraphFiles/'
detectFile = dataDir + '4_output.detect'
stressFileRight = dataDir + 'rightForks_DNAscent_forkSense_stressSignatures.bed'
stressFileLeft = dataDir + 'leftForks_DNAscent_forkSense_stressSignatures.bed'
origins = dataDir + 'origins_DNAscent_forkSense.bed'
terminations = dataDir + 'terminations_DNAscent_forkSense.bed'
bedgraphBrdU = dataDir + 'allBrdU.bedgraph'
bedgraphEdU = dataDir + 'allEdU.bedgraph'
outputDir = '/Users/shephc3/Desktop/AEbedgraphFiles/'

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
	cutChromosome = []
	cutTrackStart = []
	cutTrackEnd = []
	cutReadID = []
	cutReadStart = []
	cutReadEnd = []
	cutDirection = []
	cutFirstAnalogueLength = []
	cutSecondAnalogueLength = []
	cutDetectionLength = []
	for line in f:
		splitLine = line.rstrip().split()
		if line[0] == '#':
			#Ignore header lines
			continue
		elif splitLine[14] != '-3.000000':
			#Only collect reads that aren't identified as cutoff (-3) in the last column
			detectReads.append(splitLine[3])
			detectLength.append(float(splitLine[7]))
			firstAnalogueLength.append(float(splitLine[8]))
			secondAnalogueLength.append(float(splitLine[9]))
			chromosome.append(splitLine[0])
			trackStart.append(splitLine[1])
			trackEnd.append(splitLine[2])
			readID.append(splitLine[3])
			readStart.append(splitLine[4])
			readEnd.append(splitLine[5])
			direction.append(splitLine[6])
			detectionLength.append(float(splitLine[7]))
		if splitLine[14] == '-3.000000':
			cutoffCount += 1
			cutChromosome.append(splitLine[0])
			cutTrackStart.append(splitLine[1])
			cutTrackEnd.append(splitLine[2])
			cutReadID.append(splitLine[3])
			cutReadStart.append(splitLine[4])
			cutReadEnd.append(splitLine[5])
			cutDirection.append(splitLine[6])
			cutDetectLength.append(float(splitLine[7]))
			cutFirstAnalogueLength.append(float(splitLine[8]))
			cutSecondAnalogueLength.append(float(splitLine[9]))
			cutDetectionLength.append(float(splitLine[7]))
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

	unfilterDF = pd.DataFrame(
		{"Chromosome" : cutChromosome,
		 "Track Start" : cutTrackStart,
		 "Track End" : cutTrackEnd,
		 "Read ID" : cutReadID,
		 "Read Start" : cutReadStart,
		 "Read End" : cutReadEnd,
		 "Direction" : cutDirection,
		 "First Analogue Length" : cutFirstAnalogueLength,
		 "Second Analogue Length" : cutSecondAnalogueLength,
		 "Detect Length" : cutDetectionLength}, index=cutReadID)
##Add in 'keys' at a later date to differentiate within new DF
	f.close()
	return filterDF, unfilterDF

#This function will calculate the N50/Median given an inputted array
def calculateN50(values):
        values = np.sort(values)
        sum = np.sum(values)
        cumulative = 0
        for x in values:
                cumulative += x
                if cumulative >= sum / 2:
                        return x

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

#This function does the leg work of processing a forksense file.
#It process the file using processStressFile and turns that into a dataframe
#It then identifies reads in origin and termination files and drops those from the dataframe.
#It then returns the resulting dataframe for further processing or saving to a csv.
def filterFile(file):
	global origins
	global terminations
	filterDF, unfilterDF = processStressFile(file)
	reads = np.concatenate((processOriTerm(origins), processOriTerm(terminations)))
	for i in reads:
		if i in filterDF['Read ID'].to_numpy():
			filterDF.drop(i)
	return filterDF, unfilterDF

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

#Some variables used in file processing.
detectReads = []
cutDetectReads = []
detectLength = []
cutDetectLength = []
cutoffCount = 0

#Process the files and populate the previously defined arrays
rightDF, cutRightDF = filterFile(stressFileRight)
unfilterRightDF = pd.concat([rightDF, cutRightDF])
rightDF.to_csv(f'{outputDir}rightForksFiltered.bed', sep='\t', header=False, index=False, mode='a')

leftDF, cutLeftDF = filterFile(stressFileLeft)
unfilterLeftDF = pd.concat([leftDF, cutLeftDF])
leftDF.to_csv(f'{outputDir}leftForksFiltered.bed', sep='\t', header=False, index=False, mode='a')

frames = [rightDF, leftDF]
forkDF = pd.concat(frames)
cutForkDF = pd.concat([cutRightDF, cutLeftDF])
unfilterForkDF = pd.concat([unfilterRightDF, unfilterLeftDF])

#Now we parse the detect file which contains each read and the analog probability at each base.
#We'll ignore the probabilty of each analog and only focus on the read length and readID
#Importantly, this includes reads with EdU/BrdU tracks as well as those without.
file = open(detectFile, 'r')
readID = []
readLength = []
for line in file:
        if line[0] == '#':
                continue
        elif line[0] == '>':
                #ReadID and start/stop positions are located in lines marked with '>'
                splitLine = line.rstrip().split()
                readID.append(splitLine[0][1:])
                readLength.append(int(splitLine[3]) - int(splitLine[2]))
file.close()

#Some variables that we'll populate with reads/lengths that either contain EdU/BrdU positive tracks or not
positiveReads = []
negativeReads = []
positiveLength = []
negativeLength = []
i = 0

while i < len(readID):
        if readID[i] in detectReads:
                #If the readID shows up in either of the fork files, we'll collect those in the positive read array.
                positiveReads.append(readID[i])
                positiveLength.append(readLength[i])
        else:
	        #Otherwise they go to the negative array
                negativeReads.append(readID[i])
                negativeLength.append(readLength[i])
        i += 1

firstAnalogueLength = forkDF['First Analogue Length'].to_numpy()
secondAnalogueLength = forkDF['Second Analogue Length'].to_numpy()
unfilterfirstAnalogueLength = unfilterForkDF['First Analogue Length'].to_numpy()
unfiltersecondAnalogueLength = unfilterForkDF['Second Analogue Length'].to_numpy()
cutFirstAnalogueLength = cutForkDF['First Analogue Length'].to_numpy()
cutSecondAnalogueLength = cutForkDF['Second Analogue Length'].to_numpy()
cutDetectionLength = cutForkDF['Detect Length'].to_numpy()

#Calculates the N50/Median of each array given the function defined earlier
positiveN50 = calculateN50(positiveLength)
negativeN50 = calculateN50(negativeLength)
detectN50 = calculateN50(detectLength)
firstAnalogueN50 = calculateN50(firstAnalogueLength)
secondAnalogueN50 = calculateN50(secondAnalogueLength)
unfilterFirstN50 = calculateN50(unfilterfirstAnalogueLength)
unfilterSecondN50 = calculateN50(unfiltersecondAnalogueLength)
percentCutOff = float(len(forkDF)/len(unfilterForkDF))*100
avgPos = float(np.average(positiveLength))
avgNeg = float(np.average(negativeLength))
lengthRatio = np.average(secondAnalogueLength/firstAnalogueLength)

print("Unfiltered Data:")
print(f"Number of Reads: " + str(len(unfilterForkDF)))
print(f"Median First Analogue Length: {unfilterFirstN50}")
print(f"Median Second Analuge Length: {unfilterSecondN50}")
print(f"Positive Read Median: {positiveN50}")
print(f"Negative Read Median: {negativeN50}")
print(f"Percent Whole Reads: {percentCutOff}")

print("\nFiltered Data:")
print(f"Track N50: {detectN50}")
print(f"First Analogue Median: {firstAnalogueN50}")
print(f"Second Analogue Median: {secondAnalogueN50}")
print(f"Number of Filtered Reads: " + str(len(forkDF)))
print(f"Length Ratio: {lengthRatio}")

f = open(f"{outputDir}results.txt","w")
f.write(f"Unfiltered Data:\n"
	f" Number of reads: {len(unfilterForkDF)}\n"
	f"Median First Analogue Length: {unfilterFirstN50}\n"
	f"Median Second Analuge Length: {unfilterSecondN50}\n"
	f"Positive Read Median: {positiveN50}\n"
	f"Negative Read Median: {negativeN50}\n"
	f"Percent Reads Cut-Off: {percentCutOff}\n"
	f"\nFiltered Data:\n"
	f"Track N50: {detectN50}\n"
	f"First Analogue Median: {firstAnalogueN50}\n"
	f"Second Analogue Median: {secondAnalogueN50}\n"
	f"Number of Filtered Reads: {len(forkDF)}\n"
	f"Length Ratio: {lengthRatio}\n"
	)

plt.hist(firstAnalogueLength, bins=20, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Distribution of first analogue track lengths')
plt.grid(False)
totalCount = len(firstAnalogueLength)
plt.text(0.8, 0.8, f'Total count: {totalCount}', transform=plt.gca().transAxes, fontsize=12, ha='center')
plt.text(0.8, 0.7, f'Median size: {firstAnalogueN50}', transform=plt.gca().transAxes, fontsize=12, ha='center')
filename = f"{outputDir}First_Analogue_Distribution.png"
plt.savefig(filename)


plt.clf()

plt.hist(secondAnalogueLength, bins=20, color='skyblue', edgecolor='black')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.title('Distribution of second analogue track lengths')
plt.grid(False)
totalCount = len(secondAnalogueLength)
plt.text(0.8, 0.8, f'Total count: {totalCount}', transform=plt.gca().transAxes, fontsize=12, ha='center')
plt.text(0.8, 0.7, f'Median size: {secondAnalogueN50}', transform=plt.gca().transAxes, fontsize=12, ha='center')
filename = f"{outputDir}Second_Analogue_Distribution.png"
plt.savefig(filename)

plt.clf()

filename = f"{outputDir}LengthRatioScatter.png"
plt.scatter(secondAnalogueLength/firstAnalogueLength, detectLength, c='skyblue', label='Filtered')
plt.scatter(cutSecondAnalogueLength/cutFirstAnalogueLength, cutDetectionLength, c='red', label='Cutoff')
plt.scatter(2, 25000, c='green', label='Optimal')
plt.xlabel('Length Ratio')
plt.ylabel('Length')
plt.title('Track Length vs. Length Ratio')
plt.legend()
plt.savefig(filename)


file = open(bedgraphBrdU, 'r')
keep = False
read = ''
chromosome = []
start = []
stop = []
probability = []

for line in file:
	splitLine = line.rstrip().split()
	if splitLine[0] == 'track':
		if keep:
			makeDF('BrdU', chromosome, start, stop, probability, read)
		chromosome = []
		start = []
		stop = []
		probability = []
		read_id = line[26:62]
		keep = read_id in detectReads
		if keep:
			read = read_id
	elif keep:
		chromosome.append(splitLine[0])
		start.append(splitLine[1])
		stop.append(splitLine[2])
		probability.append(splitLine[3])

# Handle last group of data
if keep:
    makeDF('BrdU', chromosome, start, stop, probability, read)

file.close()

file = open(bedgraphEdU, 'r')
keep = False
read = ''
chromosome = []
start = []
stop = []
probability = []

for line in file:
	splitLine = line.rstrip().split()
	if splitLine[0] == 'track':
		if keep:
			makeDF('EdU', chromosome, start, stop, probability, read)
		chromosome = []
		start = []
		stop = []
		probability = []
		read_id = line[26:62]
		keep = read_id in detectReads
		if keep:
			read = read_id
	elif keep:
		chromosome.append(splitLine[0])
		start.append(splitLine[1])
		stop.append(splitLine[2])
		probability.append(splitLine[3])

# Handle last group of data
if keep:
	makeDF('EdU', chromosome, start, stop, probability, read)

file.close()