#Extracts the number of B.1.1.7 and non-B.1.1.7 genomes and the proportion of B.1.1.7 genomes from each adm2 in each week

import argparse
import pandas as pd
import datetime

#Creates a dictionary containing each date as keys and the corresponding first date in the week as values between the given firstDate and lastDate
#firstDate must be the first day of the first week, subsequent weeks are calculated from this
#Both firstDate and lastDate are given as strings
def getWeekDict(firstDate, lastDate):
    #Will be filled with dates as keys and week start dates as values
    epiWeekDict = {}

    #The first date to be examined
    startDate = datetime.date(int(firstDate.split("-")[0]), int(firstDate.split("-")[1]), int(firstDate.split("-")[2]))
    #The last date to be examined
    endDate = datetime.date(int(lastDate.split("-")[0]), int(lastDate.split("-")[1]), int(lastDate.split("-")[2]))
    #A single day incrementor
    iterator = datetime.timedelta(days = 1)
    #A week incrementor
    weekIterator = datetime.timedelta(days = 7)
    #Will be incremented with each day and when divided by 7 will have remainder 0 at the start of a new week
    epiWeekIterator = 0

    #The starting epi week, should be the first day of that epi week
    epiWeek = startDate

    #Iterate through the dates and add them to epiWeekDict
    while startDate <= endDate:
        epiWeekDict[str(startDate)] = str(epiWeek)
        epiWeekIterator += 1
        if (epiWeekIterator % 7) == 0:
            epiWeek += weekIterator
        startDate += iterator

    return(epiWeekDict)

#Used to check if a value is nan
def isNaN(string):
    return(string != string)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", help = "20210122_Genomes_UTLA_corrected_clustered_adm2.csv containing sequence metadata")
    parser.add_argument("-o", help = "Output csv file")
    args = parser.parse_args()

    #Import the sequence metadata
    metadata = pd.read_csv(args.m)

    #Open output file
    outFile = open(args.o, "w")
    outFile.write("adm2,week_start,B.1.1.7_sequences,other_sequences,proportion_sequences_B.1.1.7\n")

    #Dates as keys, week start as values
    weekDict = getWeekDict("2020-08-30", "2021-01-19")

    #All adm2s and weeks in data, used to write number of sequences
    allAdm2s = list()
    allWeeks = list()

    #adm2:week as keys, number of sequences as values
    b117 = dict()
    nonB117 = dict()

    #Iterate through the sequences and add to b117 and nonB117
    for i in range(metadata.shape[0]):
        a = metadata["aggregated_adm2"][i]
        #print(a, isNaN(a))
        if (isNaN(a) == False) and ("|" not in a):
            d = metadata["sample_date"][i].split("/")
            date = d[2] + "-" + d[1] + "-" + d[0]
            week = weekDict[date]

            if metadata["lineage"][i] == "B.1.1.7":
                if (a + ":" + week) in b117:
                    b117[a + ":" + week] += 1
                else:
                    b117[a + ":" + week] = 1
            else:
                if (a + ":" + week) in nonB117:
                    nonB117[a + ":" + week] += 1
                else:
                    nonB117[a + ":" + week] = 1

            if a not in allAdm2s:
                allAdm2s.append(a)
            
            if week not in allWeeks:
                allWeeks.append(week)
    
    #Write the sequences in each week
    for eachAdm2 in allAdm2s:
        for eachWeek in allWeeks:
            if "," in eachAdm2:
                outFile.write('"' + eachAdm2 + '",')
            else:
                outFile.write(eachAdm2 + ",")
            outFile.write(eachWeek + ",")
            aw = eachAdm2 + ":" + eachWeek
            if (aw in b117) and (aw in nonB117):
                outFile.write(str(b117[aw]) + "," + str(nonB117[aw]) + "," + str(float(b117[aw])/(float(b117[aw]) + float(nonB117[aw]))) + "\n")
            elif aw in b117:
                outFile.write(str(b117[aw]) + ",0,1\n")
            elif aw in nonB117:
                outFile.write("0," + str(nonB117[aw]) + ",0\n")
            else:
                outFile.write("0,0,0\n")
    
    outFile.close()