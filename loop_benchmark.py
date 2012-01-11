"""


"""

__author__ = 'Anton Petrov'


import re, os, numpy, logging


class LoopBenchmark():
    """
    """

    def __init__(self):
        self.scor_ifn = 'scor/1s72.html'
        self.rna3dmotif_location = 'rna3dmotif/CATALOGUE/DESC'
        self.cossmos_location = 'cossmos/results/tmp'
        self.rnajunction_location = 'rnajunction/junctions'

        self.scor = []
        self.rna3dmotif = []
        self.cossmos = []
        self.rnajunction = []

    def parse_scor(self):
        """
        """
        f = open(self.scor_ifn, 'r')
        file = f.read()
        ids = re.findall('&idElement=(.+?)#', file)
        for id in ids:
            # 1s72:0:78-80,0:94-99
            parts = id.replace('1s72:','').replace(':','/').replace('-',':').split(',')
            parts = ['/'.join(['1',x]) for x in parts]
            self.scor.append(','.join(parts))
        logging.info('Scor parsed')

    def __cossmos_get_strand(self, nts, chain):
        """
        """
        model_num = '1'
        numbers =  ':'.join([nts[0][0:nts[0].index('-')],
                             nts[-2][0:nts[-2].index('-')]])
        return '/'.join([model_num, chain, numbers])

    def parse_cossmos(self):
        """
        """
        files = os.listdir(self.cossmos_location)
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
                nts = parts[strand1].split(' ')
                lstrand = self.__cossmos_get_strand(nts, chain)
                if strand2:
                    nts = parts[strand2].split(' ')
                    rstrand = self.__cossmos_get_strand(nts, chain)
                    self.cossmos.append(','.join([lstrand, rstrand]))
                else:
                    self.cossmos.append(lstrand)
        print self.cossmos

    def parse_rna3dmotif(self):
        """
        """
        files = os.listdir(self.rna3dmotif_location)
        for ifn in files:
            chain = '/'.join(['1',ifn[5]])
            f = open(os.path.join(self.rna3dmotif_location, ifn), 'r')
            file = f.read()
            #	Bases: 118_G  119_A  120_A  121_U  122_C
            nts = [re.findall('(\d+)',x) for x in re.findall('Bases: (.+?)\s+\n', file)]
            numbers = [int(x) for x in nts[0]]
            d = numpy.diff(numbers) # difference between neighbors
            b = filter(lambda x: x>1, d) # locate chain break
            if b: # internal loop
                chbr = d.tolist().index(b[0]) # chain break
                self.rna3dmotif.append(','.join([
                    '/'.join([chain, ':'.join([nts[0][0], nts[0][chbr]]) ]),
                    '/'.join([chain, ':'.join([nts[0][chbr+1], nts[0][-1]]) ])
                                ]))
            else: # hairpin loop
                self.rna3dmotif.append('/'.join([chain, ':'.join([nts[0][0], nts[0][-1]]) ]))

        print self.rna3dmotif

    def parse_rnajunction(self):
        """Parse rnajunction pdb files to get loop intervals"""
        files = os.listdir(self.rnajunction_location)
        for ifn in files:
            f = open(os.path.join(self.rnajunction_location, ifn), 'r')
            lines = f.readlines()
            nums = []
            for line in lines:
                if line[0:4] == 'ATOM':
                    nums.append(line[23:27].rstrip().lstrip())
                    chain = '/'.join(['1',line[21:22]])
            nums = set(nums) # remove duplicates
            numbers = []
            nts = []
            [numbers.append(int(x)) for x in nums] # convert to list
            numbers.sort()
            [nts.append(str(x)) for x in numbers] # convert to strings
            d = numpy.diff(numbers) # difference between neighbors
            b = filter(lambda x: x>1, d) # locate chain breaks
            if b: # internal loop or a junction
                chbr = []
                d = d.tolist()
                [chbr.append(d.index(x)) for x in b]
                start = 0
                fragments = []
                for c in chbr: # loop over chainbreaks
                    fragments.append('/'.join([chain, ':'.join([nts[start], nts[c]]) ]))
                    start = c + 1
                self.rnajunction.append(','.join(fragments))

        print self.rnajunction
        print len(self.rnajunction)


    def parse_rloom(self):
        pass

    def parse_fr3d(self):
        pass


def main():
    """
    """
    logging.basicConfig(level=logging.DEBUG)

    L = LoopBenchmark()
#     L.parse_scor()
#     L.parse_rna3dmotif()
#     L.parse_cossmos()
    L.parse_rnajunction()


if __name__ == "__main__":
    main()