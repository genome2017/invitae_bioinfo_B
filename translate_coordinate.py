import os
import sys
import re
from argparse import ArgumentParser
from mappings import Mappings


class GenomicMapping:
    """
    Class to hold the input information about the mappings (contents from inputfile1.txt)
    :param transcriptName: string, transcript name
    :param: chromosome: string, genomic chromsome
    :param: pos: int, genomic mapping position (presumably from the sam/bam file). According to spec, this represents the
                position of the first match.
    :param: cigar: string, cigar string of alignment of transcript to genome.
    :param: orientation: string (+/-) which represents if transcript maps 5'->3' or 3'->5'
    """

    def __init__(self, transcript_name, chromosome, pos, cigar, orientation):
        self.transcript_name = transcript_name
        self.chromosome = chromosome

        try:
            self.pos = int(pos)
        except:
            raise ValueError("Invalid alignment position for "+transcript_name+" "+chromosome+" "+pos+"\n")

        if pos < 0:
            raise ValueError("Alignment position is negative for "+transcript_name+" "+chromosome+" "+pos+"\n")

        self.cigar_string = cigar
        if orientation != "+" and orientation != "-":
            raise ValueError("Invalid alignment orientation " + str(orientation))
        self.orientation = orientation



def is_valid_cigar(cigar_string):
    """
    :param cigar_string, string: Input cigar string
    :return: boolean if valid or not
    """

    matches = re.findall(r'(\d+)([A-Z]{1})', cigar_string)
    matches_2 = re.findall(r'([A-Z]{1})', cigar_string)
    is_valid = True

    if len(matches) != len(matches_2):
        # sanity check that we when we parse only the operations, we get the same number
        return False

    for m in matches:
        op_len = m[0]
        operation = m[1]
        try:
            int_test = int(op_len)
        except:
            sys.stderr.write("Operation length is not an integer: "+m[0]+" "+m[1]+"\n")
            is_valid = False

        if operation not in 'MID':
            is_valid = False
            sys.stderr.write("Invalid cigar operation: "+operation+".\n")
    return is_valid


def validate_input(args):
    if not os.path.isfile( args.genome_mapping_file):
        # TODO: at this stage can also verify format of the inputs
        return (False, "Provided genome mapping file does not exist" )

    if not os.path.isfile( args.transcript_processing_file):
        # TODO: at this stage can also verify format of the inputs
        return (False, "Provided processing file does not exist" )
    if not os.path.isdir(os.path.dirname(os.path.abspath(args.output_file))):
        return (False, "Specified parent directory for output file location does not exist")

    return (True, "")


def translate_coordinates(genome_mapping_file, processing_file, output_file):
    """
    Translate coordinates specified in processing_file based on alignments specified in genome_mapping_file. 
    Write translations to output_file
    
    :param genome_mapping_file, string: Name of file specifying alignment of transcript to genome.  See documentation 
    for file spec.
    :param processing_file, string: Name of file specifying transcripts and positions to process.  See documentation 
    for file spec.
    :param output_file, string:  Name of output file translations will be written t..
    :return: void
    """

    mappings = {}
    with open(genome_mapping_file) as in_handle:
        for line in in_handle:
            data = line.rstrip().split("\t")

            if len(data) < 4:
                sys.stderr.write("Line in genome mapping file does have atleast 4 columns\n")
                sys.stderr.write(line)
                continue

            # if not specified in the input file, assume the mapping orientation is 5'=>3'

            mapping_orientation = "+"
            if len(data) == 5:
                mapping_orientation = data[4]

            if not is_valid_cigar(data[3]):
                sys.stderr.write("Input cigar string is not valid. Skipping "+data[0]+" "+data[3])
                continue

            try:
                GM = GenomicMapping(data[0], data[1], data[2], data[3], mapping_orientation)
            except:
                sys.stderr.write("Excluding: "+data[0]+" from analysis - invalid input data.\n")
                continue
            mappings[data[0]] = GM


    o_handle = open(output_file,'w')
    with open(processing_file) as in_handle:
        for line in in_handle:
            data = line.rstrip().split("\t")

            if len(data)<2:
                sys.stderr.write("Line in processing file does have atleast 4 columns\n")
                sys.stderr.write(line)
                continue

            # TODO: check types of input
            transcript = data[0]
            query_position = int(data[1])

            # default mapping is from transcript -> genome
            mapping_direction = "TRANSCRIPT"

            if len(data) == 3:
                mapping_direction = data[2]

            if transcript not in mappings:
                sys.stderr.write("Can't find mappings for : " + data[0] + "\n")
                continue

            genome_mapping_info = mappings[transcript]
            mapping_direction = data[2]

            if mapping_direction != "TRANSCRIPT" and mapping_direction!="GENOMIC":
                sys.stderr.write ("Specification of mapping direction is not TRANSCRIPT or GENOMIC.Skipping\n")
                continue

            try:
                query_mapping = Mappings(genome_mapping_info.cigar_string, genome_mapping_info.chromosome,
                                    genome_mapping_info.pos, genome_mapping_info.orientation)
            except:
                sys.stderr.write("Could not process this mapping.  Skipping "+transcript+".\n")
                continue

            # in cases where mapped coordinate is not a tuple (for cigar operations which are not insertions),
            # return only 1 int instead of a tuple

            #TODO: move printing to it's own method
            output_line = None
            print_array = []
            print_array.append(genome_mapping_info.transcript_name)
            print_array.append(genome_mapping_info.pos)
            print_array.append(genome_mapping_info.chromosome)
            print_array.append(None) # placeholder for position
            genome_pos = genome_mapping_info.pos
            transcript_pos = None

            if mapping_direction == "GENOMIC":
                output_coordinate = Mappings.genomic_to_transcript_pos(query_position, query_mapping )

                #TODO : handle this better
                if output_coordinate is None:
                     genome_pos = "ERROR"
                else:
                    genome_pos = output_coordinate[0]

                    if output_coordinate[0] != output_coordinate[1]:
                        genome_pos = str(output_coordinate[0]) + "-" + str(output_coordinate[1])
                transcript_pos = query_position

            else:
                output_coordinate = Mappings.transcript_to_genomic_pos(query_position, query_mapping )
                if output_coordinate is None:
                    transcript_pos = "ERROR"
                else:
                    transcript_pos = output_coordinate[0]
                    if output_coordinate[0] != output_coordinate[1]:
                        transcript_pos = str(output_coordinate[0]) + "-" + str(output_coordinate[1])

            print_array[1] = genome_pos
            print_array[3] = transcript_pos

            # map all to string
            output_line = "\t".join(map(str,print_array))
            print(output_line)
            o_handle.write(output_line)


######## MAIN ###############

if __name__ == "__main__":

    parser = ArgumentParser(
        "Translate coordinates from transcripts->genome (or vice versa) based on input mapping information")

    parser.add_argument("--genome-mapping-file", required=True, dest="genome_mapping_file", help="File specifying mappings (inputfile1.txt in exercise specifications) ")
    parser.add_argument("--transcript-processing-file", required=True, dest="transcript_processing_file", help="File specifying transcripts to process (inputfile2.txt in exercise specifications) ")
    parser.add_argument("--output_file", dest="output_file", required=False, default='output.txt',
                        help="Name of output file to write results to.  Default is output.txt)")

    args = parser.parse_args()

    (is_input_valid, msg) = validate_input(args)
    if not is_input_valid:
        sys.stderr.write(msg + "\n")
        sys.exit(-1)

    translate_coordinates(args.genome_mapping_file, args.transcript_processing_file, args.output_file)