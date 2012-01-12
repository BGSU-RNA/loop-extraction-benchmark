"""


"""

__author__ = 'Anton Petrov'


import re, os, logging, csv, pdb, collections, itertools, numpy


class LoopBenchmark():
    """
    """

    def __init__(self):
        """
        """
        self.methods = ['fr3d', 'rna3dmotif', 'rnajunction', 'cossmos', 'scor',
                        'rloom']
        self.chains = ['0', '9']
        """input file locations"""
        self.scor_ifn             = 'scor/1s72.html'
        self.fr3d_ifn             = 'fr3d/fr3d_loops.csv'
        self.rloom_ifn            = 'rloom/structures/loops.txt'
        self.cossmos_location     = 'cossmos/results/tmp'
        self.rna3dmotif_location  = 'rna3dmotif/CATALOGUE/DESC'
        self.rnajunction_location = 'rnajunction/junctions'
        """initialize empty defaultdicts"""
        self.scor        = collections.defaultdict(list)
        self.fr3d        = collections.defaultdict(list)
        self.rloom       = collections.defaultdict(list)
        self.cossmos     = collections.defaultdict(list)
        self.rna3dmotif  = collections.defaultdict(list)
        self.rnajunction = collections.defaultdict(list)

    def parse_scor(self):
        """Extract loops from one html file"""
        logging.info('Parsing scor')
        duplicates = 0
        f = open(self.scor_ifn, 'r')
        file = f.read()
        ids = re.findall('&idElement=(.+?)#', file)
        for i, id in enumerate(ids):
            # 1s72:0:78-80,0:94-99
            chain = id[5]
            parts = re.findall('\d+-\d+', id)
            if not parts:
                nt = re.findall('1s72:\d:(\d+)', id)
                loop = [int(nt[0])]
            else:
                loop = []
                for part in parts:
                    nts = part.split('-')
                    loop.extend(xrange(int(nts[0]), int(nts[1])+1))
            loop.sort()
            if loop not in self.scor[chain]:
                self.__append_id('scor', chain, i)
                self.scor[chain].append(loop)
            else:
                duplicates += 1
        if duplicates > 0:
            logging.info('%i duplicates found', duplicates)

    def print_report(self):
        """
        """
        for method in self.methods:
            loops = getattr(self, method)
            logging.info('Chain 9:')
            logging.info(loops['9'])
            logging.info('Found %i loops in chain 0', len(loops['0']))
            logging.info('Found %i loops in chain 9', len(loops['9']))
            logging.info('%s', '='*50)

    def parse_cossmos(self):
        """Parse multiple custom csv files with 2x2, 3x4 etc loops"""
        logging.info('Parsing cossmos')
        files = os.listdir(self.cossmos_location)
        i = 1 # loop counter for assigning ids
        duplicates = 0
        for ifn in files:
            f = open(os.path.join(self.cossmos_location, ifn), 'r')
            lines = f.readlines()
            description = lines.pop(0)
            keys = []
            for k in description.split('#'):
                keys.append(k)
            strand1 = keys.index('Aseq_num')
            try:
                strand2 = keys.index('Bseq_num')
            except: # hairpin
                strand2 = None
            for line in lines:
                try:
                    a = line.index("0'")
                    chain = '0'
                except:
                    chain = '9'
                parts = line.split('#')
                nums = []
                # ['118-G', '119-A', '120-A', '121-U', '122-C', '']
                nts = re.findall('\d+', parts[strand1])
                [nums.append(int(x)) for x in nts]
                if strand2:
                    nts = re.findall('\d+', parts[strand2])
                    [nums.append(int(x)) for x in nts]
                nums.sort()
                if nums not in self.cossmos[chain]:
                    self.__append_id('cossmos', chain, i)
                    self.cossmos[chain].append(sorted(nums))
                    i += 1 # increment cossmos loop counter
                else:
                    duplicates += 1
        if duplicates > 0:
            logging.info('Found %i duplicates', duplicates)

    def parse_rna3dmotif(self):
        """
        """
        logging.info('Parsing rna3dmotif')
        files = os.listdir(self.rna3dmotif_location)
        duplicates = 0
        for i, ifn in enumerate(files):
            chain = ifn[5]
            f = open(os.path.join(self.rna3dmotif_location, ifn), 'r')
            file = f.read()
            #	Bases: 118_G  119_A  120_A  121_U  122_C
            nts = [re.findall('(\d+)',x) for x in re.findall('Bases: (.+?)\s+\n', file)]
            numbers = [int(x) for x in nts[0]]
            numbers.sort()
            if numbers not in self.rna3dmotif[chain]:
                self.__append_id('rna3dmotif', chain, i)
                self.rna3dmotif[chain].append(sorted(numbers))
            else:
                duplicates += 1
        if duplicates > 0:
            logging.info('Found %i duplicates', duplicates)

    def parse_rnajunction(self):
        """Parse rnajunction pdb files to get loop intervals"""
        logging.info('Parsing rnajunction')
        duplicates = 0
        files = os.listdir(self.rnajunction_location)
        for i, ifn in enumerate(files):
            if 'j2' not in ifn and 'j3' not in ifn:
                logging.info('Skipping file %s', ifn)
                continue # skip junction > j3 and kissing loops (k2)
            f = open(os.path.join(self.rnajunction_location, ifn), 'r')
            lines = f.readlines()
            chain = ifn[20]
            nums = []
            for line in lines:
                if line[0:4] == 'ATOM':
                    nums.append(line[23:27].rstrip().lstrip())
            nums = set(nums) # remove duplicates
            nts = []
            [nts.append(int(x)) for x in nums]
            nts.sort()
            if nts not in self.rnajunction[chain]:
                self.__append_id('rnajunction', chain, i)
                self.rnajunction[chain].append(sorted(nts))
            else:
                duplicates += 1
        if duplicates > 0:
            logging.info('Found %i duplicates', duplicates)

    def __append_id(self, method, chain, i):
        """
        """
        ids = getattr(self, method)
        ids['ids%s' % chain].append('%s%i' % (method, i))

    def parse_rloom(self):
        """Parse a csv file that's generated by matlab after geometric FR3D
           searches. It contains hairpin, internal and"""
        logging.info('Parsing rloom')
        duplicates = 0
        f = open(self.rloom_ifn,'r')
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            nts = re.findall('\d+', row[0])
            nums = []
            [nums.append(int(x)) for x in nts]
            nums.sort()
            chain = row[1]
            if nums not in self.rloom[chain]:
                self.__append_id('rloom', chain, i)
                self.rloom[chain].append(nums)
            else:
                duplicates += 1
        if duplicates > 0:
            logging.info('Found %i duplicates', duplicates)

    def parse_fr3d(self):
        """Parse a csv file that's generated from a mysql query from loops_all
           table for all `id` and `loop_name` where pdb = '1S72' """
        logging.info('Parsing fr3d')
        f = open(self.fr3d_ifn,'r')
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            if i == 0:
                i = 1 # skip header row
                continue
            chain = row[1][2]
            self.fr3d['ids%s' % chain].append(row[0])
            parts = re.findall('\d+:\d+', row[1])
            nts = []
            for part in parts:
                nums = part.split(':')
                nts.extend(xrange(int(nums[0]), int(nums[1])+1))
            self.fr3d[chain].append(nts)

    def integrate_results(self):
        """Fix a method, for example, fr3d, compare each loop with every other
           loop. If there is a match, report it and remove the matched loop.
           Then go to the next method and compare it with the remaining ones."""
        for chain in self.chains:
            for i, method in enumerate(self.methods):
                loops = getattr(self, method) # e.g. self.fr3d
                for j, loop in enumerate(loops[chain]):
                    line = [loops['ids%s' % chain][j]] # id
                    line.append('%s' % self.count_chainbreaks(loop)) # loop type
                    line.append(chain) # chain
                    for k in xrange(i):
                        line.append('0') # pad with empty fields for skipped methods
                    line.append(",".join(["%s" % x for x in loop])) # own nts
                    for method2 in xrange(i+1, len(self.methods)):
                        line.append(self.found_in(loop,
                                                  self.methods[method2],
                                                  chain))
                    print ''.join(['"', '","'.join(line), '"'])

    def found_in(self, query, method, chain):
        """
        """
        loops = getattr(self, method) # e.g. self.rna3dmotif
        for i, nt_range in enumerate(loops[chain]):
            common = len(set(nt_range).intersection(query))
            pos = loops[chain].index(nt_range)
            if common == len(query) and common == len(nt_range):
                del loops[chain][pos]
                del loops['ids%s' % chain][pos]
                return '1' # exact match
            elif common >= 1:
                del loops[chain][pos]
                del loops['ids%s' % chain][pos]
                return ",".join(["%s" % x for x in nt_range]) # partial match
        return '0' # no match

#     def remove_duplicates(self):
#         """Rloom lists all internal loops twice."""
#         for method in self.methods:
#             duplicates = 0
#             for chain in self.chains:
#                 m = getattr(self, method)
#                 k = m[chain][:] # make a copy
#                 len1 = len(k)
#                 k.sort()
#                 k = list(k for k,_ in itertools.groupby(k))
#                 """delete from the original list to maintain the order"""
#                 for i, loop in enumerate m[chain]:
#
#                 duplicates += len1 - len(k)
#         if duplicates > 0:
#             logging.info('Found %i duplicates in %s', duplicates, method)
#
    def count_chainbreaks(self, k):
        """Input: a list of nucleotide numbers.
           Chain breaks are places where diff btw neighbors > 1"""
        return len(filter(lambda x: x>1, numpy.diff(k).tolist()))

    def remove_multiloops(self):
        """Since we are only comparing IL, HL and J3, it's necessary to remove
           higher order loops from the analysis."""
        for method in self.methods:
            multiloops = 0
            for chain in self.chains:
                loops = getattr(self, method)
                for loop in loops[chain]:
                    if self.count_chainbreaks(loop) > 2:
                        multiloops += 1
                        """safe only after duplicates have been removed"""
                        logging.info(loop)
                        loops[chain].remove(loop)
            logging.info('Removed %i multiloops found by %s', multiloops, method)


def main():
    """
    """
    logging.basicConfig(level=logging.DEBUG)

    L = LoopBenchmark()
    L.parse_scor()
    L.parse_rna3dmotif()
    L.parse_cossmos()
    L.parse_rnajunction()
    L.parse_fr3d()
    L.parse_rloom()
    L.print_report()
#     L.remove_duplicates()
    L.remove_multiloops()
    L.integrate_results()


if __name__ == "__main__":
    main()