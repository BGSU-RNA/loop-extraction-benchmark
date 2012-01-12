"""


"""

__author__ = 'Anton Petrov'


import re, os, numpy, logging, csv, pdb, collections


class LoopBenchmark():
    """
    """

    def __init__(self):

        self.methods = ['fr3d', 'rna3dmotif', 'rnajunction', 'cossmos', 'scor', 'rloom']
        self.scor_ifn = 'scor/1s72.html'
        self.rna3dmotif_location = 'rna3dmotif/CATALOGUE/DESC'
        self.cossmos_location = 'cossmos/results/tmp'
        self.rnajunction_location = 'rnajunction/junctions'
        self.fr3d_ifn = 'fr3d/fr3d_loops.csv'
        self.rloom_ifn = 'rloom/structures/loops.txt'

        self.scor = collections.defaultdict(list)
        self.rna3dmotif = collections.defaultdict(list)
        self.cossmos = collections.defaultdict(list)
        self.rnajunction = collections.defaultdict(list)
        self.fr3d = collections.defaultdict(list)
        self.rloom = collections.defaultdict(list)

    def parse_scor(self):
        """
        """
        logging.info('Parsing scor')
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
            self.__append_id('scor', chain, i)
            self.scor[chain].append(sorted(loop))

    def print_report(self):
        """
        """
        for method in self.methods:
            loops = getattr(self, method)
            logging.info('Chain 9:')
            print loops['9']
            logging.info('Found %i loops in chain 0', len(loops['0']))
            logging.info('Found %i loops in chain 9', len(loops['9']))

    def parse_cossmos(self):
        """
        """
        logging.info('Parsing cossmos')
        files = os.listdir(self.cossmos_location)
        i = 1
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
            except:
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
                self.__append_id('cossmos', chain, i)
                self.cossmos[chain].append(sorted(nums))
                i += 1

    def parse_rna3dmotif(self):
        """
        """
        logging.info('Parsing rna3dmotif')
        files = os.listdir(self.rna3dmotif_location)
        for i, ifn in enumerate(files):
            chain = ifn[5]
            f = open(os.path.join(self.rna3dmotif_location, ifn), 'r')
            file = f.read()
            #	Bases: 118_G  119_A  120_A  121_U  122_C
            nts = [re.findall('(\d+)',x) for x in re.findall('Bases: (.+?)\s+\n', file)]
            numbers = [int(x) for x in nts[0]]
            self.__append_id('rna3dmotif', chain, i)
            self.rna3dmotif[chain].append(sorted(numbers))

    def parse_rnajunction(self):
        """Parse rnajunction pdb files to get loop intervals"""
        logging.info('Parsing rnajunction')
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
            self.__append_id('rnajunction', chain, i)
            self.rnajunction[chain].append(sorted(nts))

    def __append_id(self, method, chain, i):
        a = getattr(self, method)
        self.rnajunction['ids%s' % chain].append('%s%i' % (method, i))

    def parse_rloom(self):
        """
        """
        logging.info('Parsing rloom')
        f = open(self.rloom_ifn,'r')
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            nts = re.findall('\d+', row[0])
            nums = []
            [nums.append(int(x)) for x in nts]
            chain = row[1]
            self.__append_id('rloom', chain, i)
            self.rloom[chain].append(nums)

    def parse_fr3d(self):
        """
        """
        logging.info('Parsing fr3d')
        f = open(self.fr3d_ifn,'r')
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            if i == 0:
                i = 1 # skip header row
                continue
            chain = row[1][2]
            self.fr3d['ids'+chain].append(row[0])
            parts = re.findall('\d+:\d+', row[1])
            nts = []
            for part in parts:
                nums = part.split(':')
                nts.extend(xrange(int(nums[0]), int(nums[1])+1))
            self.fr3d[chain].append(nts)

    def integrate(self):

        for i, method in enumerate(self.methods):
            logging.info('%s', method)
            loops = getattr(self, method) # e.g. self.fr3d
            for j, loop in enumerate(loops['0']):
                line = [",".join(["%s" % x for x in loop])]
                for method2 in xrange(i+1, len(self.methods)):
                    line.append(self.found_in(loop, self.methods[method2]))
#                 pdb.set_trace()
                print ''.join(['"', '","'.join(line), '"'])

    def found_in(self, query, method):
        """
        """
        chain = '0'
        loops = getattr(self, method) # e.g. self.rna3dmotif
        for i, nt_range in enumerate(loops[chain]):
            common = len(set(nt_range).intersection(query))
            if common == len(query) and common == len(nt_range):
                loops[chain].remove(nt_range)
                return '1' # exact match
            elif common > 1:
                loops[chain].remove(nt_range)
                return ",".join(["%s" % x for x in nt_range]) # partial match
        return '0'

    def process_duplicates(self):
        pass

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
    L.integrate()
#     pdb.set_trace()

if __name__ == "__main__":
    main()