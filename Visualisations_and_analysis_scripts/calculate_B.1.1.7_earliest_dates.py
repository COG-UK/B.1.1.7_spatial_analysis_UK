#Calculates the earliest B.1.1.7 sequence in each adm2
#Uses sequence data in 20210122_Genomes_UTLA_corrected_clustered_adm2.csv
#Uses the locations in column aggregated_adm2 as the adm2

import argparse
import pandas as pd
import datetime

#Creates a dictionary containing each date as keys and their epi weeks as values between the given firstDate and lastDate
#firstDate must be the first day of the firstEpiWeek, which is the Sunday of that week
#Takes the week to start the count on, which is the week firstDate is in
#firstDate is given as a string, e.g. "2020-09-20" while lastDate is given as a date
def getEpiWeekDict(firstDate, lastDate, firstEpiWeek):
    #Will be filled with dates as keys and epi weeks as values
    epiWeekDict = {}

    #The first date to be examined
    startDate = datetime.date(int(firstDate.split("-")[0]), int(firstDate.split("-")[1]), int(firstDate.split("-")[2]))
    #The last date to be examined
    endDate = lastDate
    #A single day incrementor
    iterator = datetime.timedelta(days = 1)
    #Will be incremented with each day and when divided by 7 will have remainder 0 at the start of a new epi week
    epiWeekIterator = 0

    #The starting epi week, should be the first day of that epi week
    epiWeek = int(firstEpiWeek)

    #Iterate through the dates and add them to epiWeekDict
    while startDate <= endDate:
        epiWeekDict[startDate] = epiWeek
        epiWeekIterator += 1
        if (epiWeekIterator % 7) == 0:
            epiWeek += 1
        startDate += iterator

    return(epiWeekDict)

#Used to check if a value is nan
def isNaN(string):
    return(string != string)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", help = "20210122_Genomes_UTLA_corrected_clustered_adm2.csv containing sequence metadata")
    parser.add_argument("-o", help = "Output csv file name")
    args = parser.parse_args()
    
    outFile = open(args.o, "w")
    outFile.write("adm2,earliest_sequence,earliest_epi_week\n")

    #Import the sequence metadata
    metadata = pd.read_csv(args.s)

    #adm2s as keys, earliest sequence as values
    adm2Seq = dict()

    #Iterate through the sequences and identify the earliest B.1.1.7 sequence from each adm2
    for i in range(metadata.shape[0]):
        if metadata["B117"][i]:
            a = metadata["aggregated_adm2"][i]

            #Check if the adm2 is ambiguous
            if not isNaN(a):
                if ("|" not in a) and (a != "Needs_manual_curation"):
                    d = metadata["sample_date"][i].split("/")
                    date = datetime.date(int(d[2]), int(d[1]), int(d[0]))

                    if a in adm2Seq:
                        if date < adm2Seq[a]:
                            adm2Seq[a] = date
                    else:
                        adm2Seq[a] = date

    #Extract an epi week dictionary for the required dates
    ewDict = getEpiWeekDict("2020-09-20", max(adm2Seq.values()), 39)

    #Write the adm2s and earliest dates
    for adm2 in adm2Seq:
        if "," in adm2:
            outFile.write('"' + adm2 + '"')
        else:
            outFile.write(adm2)
        outFile.write("," + str(adm2Seq[adm2]) + "," + str(ewDict[adm2Seq[adm2]]) + "\n")

    outFile.close()