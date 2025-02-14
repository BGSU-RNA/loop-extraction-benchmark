#!/usr/bin/python

"""

Download pdb files from a results html page generated by the RLooM webserver.
It's necessary to download the pdbs in order to further process them in matlab.
RLooM uses its own internal indexes to describe loops, so one has to perform
geometric FR3D searches to match them to their original nucleotides.

"""

__author__ = 'Anton Petrov'

import re, os
from urllib2 import Request, urlopen, URLError, HTTPError

rloom_results_html = 'RNA Loop Modeling.html'
output_dir         = 'structures'
pdb_url            = 'http://rloom.mpimp-golm.mpg.de/structures/'

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

f = open(rloom_results_html, 'r')

# example pattern: frid=pdb1s7210.n198-225"
pattern = re.compile('frid=(.+?)"')

matches = problems = downloaded = 0

for line in f:
    if '<th class="cltable">Segment' in line:
        # don't download "Segments", only internal, hairpin and multiloops
        break

    m = re.findall(pattern, line)
    if m:
        matches += 1;

        print m[0]
        ofn = os.path.join(output_dir, m[0] + '.pdb')
        if os.path.exists(ofn): # skip already created files
            continue

        req = Request(''.join([pdb_url, m[0], '.fpdb']))
        try:
            p = urlopen(req)
            downloaded += 1
        except HTTPError, e:
            problems += 1
            print 'Error code: ', e.code
            continue
        except URLError, e:
            problems += 1
            print 'Reason: ', e.reason
            continue

        result = p.read()
        out = open(ofn, 'w')
        out.write(result)
        out.close()
        print 'Downloaded file %s' % ofn

f.close()
print 'Downloaded %i files out of %i' % (downloaded, matches)
print '%i files could not be downloaded' % problems