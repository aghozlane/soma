#!/usr/bin/env python
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html


"""Trim reads."""


#For argument parsing
import os
import sys
import argparse

__author__ = "Mathieu Almeida, Amine Ghozlane"
__copyright__ = "Copyright 2014, INRA"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Mathieu Almeida, Amine Ghozlane"
__email__ = "mathieu.almeida@jouy.inra.fr, amine.ghozlane@jouy.inra.fr"
__status__ = "Developpement"

#python trimReads.py -f <read_file> -s <read_size> -a <adapted_trimming>

def isfile(path):
    """Check if path is an existing file.
      Arguments:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

#=============================================================================
# Parameters parsing
#=============================================================================

"""
@brief Parameters Parsing
@param NULL
@return options.fastaFile, options.cutSize
"""
def config_parameters():
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h".format(sys.argv[0]))
    parser.add_argument("-f", "--readFile", dest="readFile", required=True,
                        type=isfile, help="read file in fasta format")
    parser.add_argument("-q", "--qualityFile", dest="qualityFile", type=isfile,
                        default=None, help="quality read file in fasta format")
    parser.add_argument("-s", "--trimSize", dest="trimSize", type=int,
                        help="read size asked", default = '15')
    parser.add_argument("-a", "--adapted", dest="adapted", default=False,
                        action='store_true', help="if adapted trimming asked")
    parser.add_argument("-t", "--threshold", dest="threshold", default = '9',
                        type=int, help="quality threshold")
    parser.add_argument("-c", "--correct", dest="correct", default=False,
                        action='store_true',
                        help="correct position under quality threshold")
    parser.add_argument("-o", "--output", dest="output", default=None,
                        type=str,
                        help="correct position under quality threshold")
    return parser.parse_args()

#END config_parameters

#==============
#trim reads
#==============
def trimReads(readFile, qualityFile, trimSize, output_file):
    """
    """
    if not output_file:
        output_file = (readFile.split('.')[0] + '_trim_' + str(trimSize) + 
                       '.' + readFile.split('.')[1])
    with open(readFile, 'rt') as filinRead:
        if qualityFile:
            filinQual = open(qualityFile, 'rt')
        with open(output_file, 'wt') as filoutRead:
            if qualityFile:
                filoutQual = open(qualityFile.split('.')[0] + '_trim_' + str(trimSize)
                                  + '.' + qualityFile.split('.')[1], 'wt')
            for lineRead in filinRead:
                if lineRead[0] == '>':
                    filoutRead.write(lineRead)
                else:
                    lineRead = lineRead.rstrip('\n')
                    #warning : add an extra position to include the first T
                    filoutRead.write(lineRead[:-trimSize] + os.linesep)
            if qualityFile:
                for lineQual in filinQual:
                    if lineQual[0] == '>':
                        filoutQual.write(lineQual)
                    else:
                        listQuality = []
                        qual = ''
                        listQuality = lineQual.split(' ')[:-trimSize]
                        qual = ' '.join(listQuality)
                        filoutQual.write(qual + os.linesep)
        if qualityFile:
            filinQual.close()
        if qualityFile:
            filoutQual.close()


#==============
#trim reads
#==============
def trimAdaptedReads(readFile, qualityFile, threshold, output_file):
    if not output_file:
        output_file = (readFile.split('.')[0] + '_trimAdapt_' +
                      readFile.split('.')[1])
    with open(readFile, 'rt') as filinRead:
        with open(qualityFile, 'rt') as filinQual:
            with open(output_file, "wt") as filoutRead:
                with open(qualityFile.split('.')[0] + '_trimAdapt_' +
                                  qualityFile.split('.')[1], 'w') as filoutQual:
                    ReadLines = filinRead.readlines()
                    QualLines = filinQual.readlines()
                    #for iLine in xrange(min(len(filinRead), len(filinQual))):
                    for iLine in range(len(ReadLines)):
                        if not (ReadLines[iLine][0]=='#'):
                            if not (ReadLines[iLine][0]=='>'):
                                #warning : add an extra position to include the first T
                                
                                listQuality = []
                                listQuality = QualLines[iLine].split(' ')[:-1]
                                sizeSeq = 0
                                error=0
                                for aQual in listQuality:
                                    if int(aQual) >= threshold:
                                        sizeSeq = sizeSeq + 1
                                    else:
                                        error = error + 1 
                                    if error > 1:
                                        break
                                if (sizeSeq >= 30 and
                                    ('.' not in ReadLines[iLine][0:sizeSeq])):
                                    #extract header
                                    filoutRead.write(ReadLines[iLine-1])
                                    filoutQual.write(QualLines[iLine-1])
                                    #extract sequence
                                    filoutRead.write(ReadLines[iLine][0:sizeSeq]
                                                     + os.linesep)
                                    qual = ' '.join(listQuality[0:sizeSeq-1])
                                    filoutQual.write(qual + os.linesep)


#==============
#trim reads
#==============
def CorrectReads(readFile, qualityFile, threshold, output_file):
    """
    """
    if not output_file:
        output_file = (readFile.split('.')[0] + '_correct_' + str(threshold)
                      + '.' + readFile.split('.')[1])
    with open(readFile, 'rt') as filinRead:
        with open(qualityFile, 'rt') as filinQual:
            with open(output_file, 'wt') as filoutRead:
                with open(qualityFile.split('.')[0] + '_correct_' + str(threshold)
                          + '.' + qualityFile.split('.')[1], 'wt') as filoutQual:
                    ReadLines = filinRead.readlines()
                    QualLines = filinQual.readlines()
                    #for iLine in xrange(min(len(filinRead), len(filinQual))):
                    for iLine in xrange(len(ReadLines)):
                        if not ReadLines[iLine][0]=='>':
                            #warning : add an extra position to include the first T
                            listQuality = []
                            listColor = []
                            listQuality = QualLines[iLine].split(' ')[:-1]
                            listColor = ReadLines[iLine][:-1]
                            sizeSeq = 0
                            nbWrong = 0
                            aSeq = ['T']
                            for aQualPos in range(len(listQuality)):
                                if int(listQuality[aQualPos]) < threshold:
                                    listQuality[aQualPos] = '-1'
                                    aSeq.append('.')
                                else:
                                    aSeq.append(listColor[aQualPos+1])
                            #write header
                            filoutRead.write(ReadLines[iLine-1])
                            filoutQual.write(QualLines[iLine-1])
                            #write sequence
                            qual = ' '.join(listQuality)
                            filoutRead.write(''.join(aSeq) + os.linesep)
                            filoutQual.write(qual + os.linesep)


#===================================
#MAIN
#===================================
"""
@brief MAIN 
@param BlastFile 
@return NULL
"""
def main():
    """Main program
    """
    fastaFile = ''
    cutSize = 0
    args = config_parameters()
    if not args.adapted:
        print('trimming mode')
        trimReads(args.readFile, args.qualityFile, args.trimSize, args.output)
    if args.adapted and not args.correct:
        print('adapted trimming mode')
        trimAdaptedReads(readFile, args.qualityFile, args.threshold,
                         args.output)
    if not args.adapted and args.correct:
        print('corrected mode')
        CorrectReads(args.readFile, args.qualityFile, args.threshold,
                     args.output)
#END MAIN                


if __name__=="__main__":
    main()
