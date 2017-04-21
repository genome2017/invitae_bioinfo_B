import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))

from nose import with_setup
from nose.tools import nottest
from translate_coordinate import Mappings

class TestTranscriptToGenomicInvitaeInput:
    mappings={}

    def __init__ (self):
        self.initialized = True

    def setup(self):
        #example1
        cigar = '8M7D6M2I2M11D7M'
        chr = 'CHR1'
        genomic_pos = 3
        mapping_orientation = '+'
        TR1 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR1']=TR1

        #example2
        cigar = '20M'
        chr = 'CHR2'
        genomic_pos = 10
        mapping_orientation = '+'
        TR2 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR2']=TR2

    def test(self):

        # plus strand tests
        coord1 = Mappings.transcript_to_genomic_pos(4,self.mappings['TR1'])
        coord2 = Mappings.transcript_to_genomic_pos(0,self.mappings['TR2'])
        coord3 = Mappings.transcript_to_genomic_pos(13,self.mappings['TR1'])
        coord4 = Mappings.transcript_to_genomic_pos(10,self.mappings['TR2'])
        assert coord1 == (7,7)
        assert coord2 == (10,10)
        assert coord3 == (23,23)
        assert coord4 == (20,20)


class TestTranscriptToGenomic:
    mappings={}

    def __init__ (self):
        self.initialized = True

    def setup(self):
        #example1
        cigar = '8M7D6M2I2M11D7M'
        chr = 'CHR1'
        genomic_pos = 3
        mapping_orientation = '+'
        TR1 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR1']=TR1

        #example2
        cigar = '20M'
        chr = 'CHR2'
        genomic_pos = 10
        mapping_orientation = '+'
        TR2 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR2']=TR2

        #example3
        cigar = '8M7D6M2I2M11D7M'
        chr = 'CHR1'
        genomic_pos = 3
        mapping_orientation = '-'
        TR1 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR3']=TR1

        #example4
        cigar = '20M'
        chr = 'CHR2'
        genomic_pos = 10
        mapping_orientation = '-'
        TR2 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR4']=TR2

    def test(self):

        # plus strand tests
        coord1 = Mappings.transcript_to_genomic_pos(4,self.mappings['TR1'])
        coord2 = Mappings.transcript_to_genomic_pos(0,self.mappings['TR2'])
        coord3 = Mappings.transcript_to_genomic_pos(13,self.mappings['TR1'])
        coord4 = Mappings.transcript_to_genomic_pos(10,self.mappings['TR2'])

        str1 = "TR1\t4\tCHR1\t7"
        str2 = "TR2\t0\tCHR2\t10"
        str3 = "TR1\t13\tCHR1\t23"
        str4 = "TR2\t10\tCHR2\t20"
        assert coord1 == (7,7)
        assert coord2 == (10,10)
        assert coord3 == (23,23)
        assert coord4 == (20,20)

        # minus strand tests
        coord5 = Mappings.transcript_to_genomic_pos(4,self.mappings['TR3']) #input example
        coord6 = Mappings.transcript_to_genomic_pos(0,self.mappings['TR3']) #input example
        coord7 = Mappings.transcript_to_genomic_pos(13,self.mappings['TR3']) #input example
        coord8 = Mappings.transcript_to_genomic_pos(24,self.mappings['TR3']) # boundary

        str5 = "TR3\t4\tCHR1\t39"
        str6 = "TR3\t0\tCHR2\43"
        str7 = "TR3\t13\tCHR1\t25"
        str8 = "TR3\t24\tCHR2\t3"
        print(str(coord7))
        assert coord5 == (39,39)
        assert coord6 == (43,43)
        assert coord7 == (21,21)
        assert coord8 == (3,3)

class TestIndelsTranscriptToGenomic:
    mappings={}

    def __init__ (self):
        self.initialized = True

    def setup(self):
        #example1
        cigar = '8M7D6M2I2M11D7M'
        chr = 'CHR1'
        genomic_pos = 3
        mapping_orientation = '+'
        TR1 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR1']=TR1

        #example2
        cigar = '20M'
        chr = 'CHR2'
        genomic_pos = 10
        mapping_orientation = '+'
        TR2 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR2']=TR2

        #example3
        cigar = '8M7D6M2I2M11D7M'
        chr = 'CHR1'
        genomic_pos = 3
        mapping_orientation = '-'
        TR3 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR3']=TR3

        #example4
        cigar = '20M'
        chr = 'CHR2'
        genomic_pos = 10
        mapping_orientation = '-'
        TR4 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR4']=TR4

    def test(self):

        # plus strand tests
        coord1 = Mappings.transcript_to_genomic_pos(17,self.mappings['TR1'])	# deleltion
        coord2 = Mappings.transcript_to_genomic_pos(18,self.mappings['TR1']) # boundary to an indel
        coord3 = Mappings.transcript_to_genomic_pos(15,self.mappings['TR1'])	#insertion
        coord4 = Mappings.transcript_to_genomic_pos(7,self.mappings['TR1'])
        assert coord1 == (25,25)
        assert coord2 == (37,37)
        assert coord3 == (23,24)
        assert coord4 == (10,10)


        # minus strand tests
        coord5 = Mappings.transcript_to_genomic_pos(6,self.mappings['TR3']) #boundary of indel
        coord6 = Mappings.transcript_to_genomic_pos(7,self.mappings['TR3']) #boundary of indel
        coord7 = Mappings.transcript_to_genomic_pos(9,self.mappings['TR3']) #insertion
        coord8 = Mappings.transcript_to_genomic_pos(10,self.mappings['TR3']) #insertion
        coord9 = Mappings.transcript_to_genomic_pos(17,self.mappings['TR3']) # boundary

        #PROBLEMS HERE
        print(str(coord5))
        print(str(coord6))
        print(str(coord7))
        print(str(coord8))
        print(str(coord9))

        assert coord5 == (37,37)
        assert coord6 == (25,25)
        assert coord7 == (23,24)
        assert coord8 == (23,24)
        assert coord9 == (10,10)


class TestGenomicToTranscript:
    mappings={}

    def __init__ (self):
        self.initialized = True

    def setup(self):
        #example1
        cigar = '8M7D6M2I2M11D7M'
        chr = 'CHR1'
        genomic_pos = 3
        mapping_orientation = '+'
        TR1 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR1']=TR1

        #example2
        cigar = '20M'
        chr = 'CHR2'
        genomic_pos = 10
        mapping_orientation = '+'
        TR2 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR2']=TR2

        #example3
        cigar = '8M7D6M2I2M11D7M'
        chr = 'CHR1'
        genomic_pos = 3
        mapping_orientation = '-'
        TR3 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR3']=TR3

        #example4
        cigar = '20M'
        chr = 'CHR2'
        genomic_pos = 10
        mapping_orientation = '-'
        TR4 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR4']=TR4

    def test(self):
        coord1 = Mappings.genomic_to_transcript_pos(7,self.mappings['TR1'])
        coord2 = Mappings.genomic_to_transcript_pos(10,self.mappings['TR2'])
        coord3 = Mappings.genomic_to_transcript_pos(23,self.mappings['TR1'])
        coord4 = Mappings.genomic_to_transcript_pos(20,self.mappings['TR2'])

        # test the same coordinates as in the inputfiles provided
        # if we supply the same transcript coordinates, do we get back the original genome coordinates?
        str1 = "TR1\t4\tCHR1\t7"
        str2 = "TR2\t0\tCHR2\t10"
        str3 = "TR1\t13\tCHR1\t23"
        str4 = "TR2\t10\tCHR2\t20"
        assert coord1 == (4,4)
        assert coord2 == (0,0)
        assert coord3 == (13,13)
        assert coord4 == (10,10)

        # minus strand examples
        coord5 = Mappings.genomic_to_transcript_pos(43,self.mappings['TR3']) # boundary
        coord6 = Mappings.genomic_to_transcript_pos(3,self.mappings['TR3']) # boundar
        coord7 = Mappings.genomic_to_transcript_pos(20,self.mappings['TR3'])
        coord8 = Mappings.genomic_to_transcript_pos(24,self.mappings['TR3'])

        # test the same coordinates as in the inputfiles provided
        # if we supply the same transcript coordinates, do we get back the original genome coordinates?
        str5 = "TR1\t4\tCHR1\t7"
        str6 = "TR2\t0\tCHR2\t10"
        str7 = "TR1\t13\tCHR1\t23"
        str8 = "TR2\t10\tCHR2\t20"
        assert coord5 == (0,0)
        assert coord6 == (24,24)
        assert coord7 == (14,14)
        assert coord8 == (8,8)

class TestIndelsGenomicToTranscript:
    mappings={}

    def __init__ (self):
        self.initialized = True

    def setup(self):
        #example1
        cigar = '8M7D6M2I2M11D7M'
        chr = 'CHR1'
        genomic_pos = 3
        mapping_orientation = '+'
        TR1 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR1']=TR1

        #example2
        cigar = '20M'
        chr = 'CHR2'
        genomic_pos = 10
        mapping_orientation = '+'
        TR2 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR2']=TR2

        #example3
        cigar = '8M7D6M2I2M11D7M'
        chr = 'CHR1'
        genomic_pos = 3
        mapping_orientation = '-'
        TR3 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR3']=TR3

        #example4
        cigar = '20M'
        chr = 'CHR2'
        genomic_pos = 10
        mapping_orientation = '-'
        TR4 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR4']=TR4

    def test(self):
        coord1 = Mappings.genomic_to_transcript_pos(10,self.mappings['TR1'])
        coord2 = Mappings.genomic_to_transcript_pos(18,self.mappings['TR1'])
        coord3 = Mappings.genomic_to_transcript_pos(12,self.mappings['TR1'])
        coord4 = Mappings.genomic_to_transcript_pos(23,self.mappings['TR1'])
        coord5 = Mappings.genomic_to_transcript_pos(24,self.mappings['TR1'])
        coord6 = Mappings.genomic_to_transcript_pos(30,self.mappings['TR1'])

        # test the same coordinates as in the inputfiles provided
        # if we supply the same transcript coordinates, do we get back the original genome coordinates?
        assert coord1 == (7,7)
        assert coord2 == (8,8)
        assert coord3 == (7,8)
        assert coord4 == (13,13)
        assert coord5 == (16,16)
        assert coord6 == (17,18)

        # minus strand examples
        coord7 = Mappings.genomic_to_transcript_pos(43,self.mappings['TR3']) # boundary
        coord8 = Mappings.genomic_to_transcript_pos(3,self.mappings['TR3']) # boundary
        coord9 = Mappings.genomic_to_transcript_pos(15,self.mappings['TR3'])
        coord10 = Mappings.genomic_to_transcript_pos(24,self.mappings['TR3'])
        coord11 = Mappings.genomic_to_transcript_pos(23,self.mappings['TR3'])
        coord12 = Mappings.genomic_to_transcript_pos(6,self.mappings['TR3'])
        coord13 = Mappings.genomic_to_transcript_pos(30,self.mappings['TR3'])
        print(coord12)
        # test the same coordinates as in the inputfiles provided
        # if we supply the same transcript coordinates, do we get back the original genome coordinates?
        assert coord7 == (0,0)
        assert coord8 == (24,24)
        assert coord9 == (16,17)
        assert coord10 == (8,8)
        assert coord11== (11,11)
        assert coord12 == (21,21)
        assert coord13 == (6,7)


@nottest
class TestInvalidInput:
    mappings={}

    def __init__ (self):
        self.initialized = True

    def setup(self):
        #example1
        cigar = '8M7D6M2I2M11D7M'
        chr = 'CHR1'
        genomic_pos = 3
        mapping_orientation = '+'
        TR1 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR1']=TR1

        #example2
        cigar = '20M'
        chr = 'CHR2'
        genomic_pos = 10
        mapping_orientation = '+'
        TR2 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR2']=TR2

        #example3 - invalid length
        cigar = '0M'
        chr = 'CHR2'
        genomic_pos = 10
        mapping_orientation = '+'
        TR2 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR2']=TR2

        #example4 - invalid operation
        cigar = '10Q'
        chr = 'CHR2'
        genomic_pos = 10
        mapping_orientation = '+'
        TR2 = Mappings(cigar, chr, genomic_pos, mapping_orientation)
        self.mappings['TR2']=TR2

    def test(self):
        coord1 = Mappings.transcript_to_genomic_pos(27,self.mappings['TR1']) # too long transcript
        coord2 = Mappings.genomic_to_transcript_pos(50,self.mappings['TR1']) # too long genomic
        coord3 = Mappings.transcript_to_genomic_pos(3,self.mappings['TR3']) # invalid cigar length
        coord4 = Mappings.transcript_to_genomic_pos(0,self.mappings['TR4']) # invalid cigar op

        str1 = "TR1\t4\tCHR1\t7"
        str2 = "TR2\t0\tCHR2\t10"
        str3 = "TR1\t13\tCHR1\t23"
        str4 = "TR2\t10\tCHR2\t20"
        print(str(coord4))
        assert coord1 == 7
        assert coord2 == 10
        assert coord3 == 23
        assert coord4 == 20

