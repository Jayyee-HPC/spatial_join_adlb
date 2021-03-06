#!/usr/bin/env python

import sys

adlbch = open('adlb.h')
adlbfh = open('adlbf.h','w')
print >>adlbfh,''
print >>adlbfh,'      ! ** do NOT edit this file; it is automatically generated'
print >>adlbfh,''
for line in adlbch:
    line = line.rstrip()
    splitLine = line.split()
    if 'ADLB_VERSION' in line:
        if 'ADLB_VERSION_NUMBER' in line:
            versionNumber = splitLine[2]
        elif 'ADLB_VERSION_DATE' in line:
            versionDate = splitLine[2]
        else:
            version = splitLine[2]
    elif line.startswith('#define')  and  'ADLB_' in line:
        splitLine[2] = splitLine[2].replace('(','')    # fix things like  (-999999999)
        splitLine[2] = splitLine[2].replace(')','')    ##
        if splitLine[2].isdigit()  or  splitLine[2][1:].isdigit():  # strip '-' if nec
            print >>adlbfh,'      integer,  parameter ::      ' + ' '*40 + '&'
        else:
            print >>adlbfh,'      character(4),  parameter :: ' + ' '*40 + '&'
        print >>adlbfh,'     &    %s = %s' % (splitLine[1],splitLine[2])

print >>adlbfh,''
print >>adlbfh,'      character(4),   parameter ::' + ' '*40 + '&'
print >>adlbfh,'     &    ADLB_VERSION      = "%s%s"' % (version[4],versionNumber)
print >>adlbfh,'      character(99),  parameter ::' + ' '*40 + '&'
print >>adlbfh,'     &    ADLB_VERSION_TEXT = "%s%s - %s"' % (version,versionNumber,versionDate)
