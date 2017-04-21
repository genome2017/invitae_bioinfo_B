import os
import sys
import re


class Mappings:
    """
    Class which holds information about how transcript aligns to the reference align to each other based on an input cigar string.  Also called to convert transcript to genomic coordinates and vice versa.
    :param: cigar_string, string.  Cigar string which represents the alignment of transcript to reference
    :param: genomic_chr, string.  Genomic chromosome on which the mapping is found
    :param: genomic_mapping_pos, int.  1-based leftmost mapping POSition of the first matching base (from SAM spec)
    :param: alignment_orientation, string.  Indicates if alignment is on forward or reverse strand. Valid values are + and -.
    
    """
    def __init__(self, cigar_string, genomic_chr, genomic_mapping_pos, alignment_orientation):

        self.cigar_string = cigar_string
        self.genomic_chr = genomic_chr
        self.genomic_mapping_pos = genomic_mapping_pos

        self.reference_ranges = []  # array of SequenceRange objects representing ranges in the reference,

        # one SequenceRange for each cigar operation
        self.query_ranges = []  # same as reference_ranges but for ranges in the transcript
        self.cigar_operations = []  # array of CigarOperations objects for every cigar operation

        self.is_transcript_forward = True

        if alignment_orientation == "-":
            self.is_transcript_forward = False

        if not cigar_string:
            raise ValueError("Cigar string not available - cannot process.")

        # populate arrays
        self.populate_cigar_operations()
        self.populate_ranges()


    def populate_cigar_operations(self):
        """
        Parses the input cigar string and generates a CigarOperation object for every operation.
        :return: void.  
        """

        matches = re.findall(r'(\d+)([A-Z]{1})', self.cigar_string)

        for m in matches:
            operation = m[1]
            try:
                op_len = int(m[0])
            except:
                raise ValueError("Cigar length is not an int: " + m[0])

            CO = CigarOperation(op_len, operation)
            self.cigar_operations.append(CO)


    def increment_indices(self, cigar_operation, query_start, reference_start):
        """
        Given a cigar string and ints representing the start of a new range, find stop coordinate for 
        the range.   
    
        :param cigar_operation, string: CigarOperations object
        :param query_start: Start position of the next range for query.  
        :param reference_start: Start position of the next range for reference.
        :return:  integers representing start/stop positions for query/reference range
        """

        # note about query/range start positions. We assume a new cigar operation extends the
        # last range. e.g. if last range was [6 10], we assume the start coordinate of the next range starts at 11.
        # This is not true in the case of deletions, in which case we adjust the coordinate in the code below.

        query_end = None
        reference_end = None

        if cigar_operation.operation in 'M':
            query_end = query_start + cigar_operation.op_length - 1
            reference_end = reference_start + cigar_operation.op_length - 1
        elif cigar_operation.operation == 'I':
            reference_start -= 1
            query_end = query_start + cigar_operation.op_length - 1
            reference_end = reference_start
        elif cigar_operation.operation in "D":
            query_start -= 1
            query_end = query_start
            reference_end = reference_start + cigar_operation.op_length - 1
        else:
            raise Exception("Cigar operation not supported: " + cigar_operation.operation + "\n")

        return query_start, query_end, reference_start, reference_end


    def populate_ranges(self):
        """
        Based on CIGAR operations delineate the corresponding SequenceRanges for query and reference 
        :return: void
        """

        assert len(self.cigar_operations) > 0

        current_query_start = 0
        current_reference_start = self.genomic_mapping_pos

        # when we have a match or mismatch (M), both query and ref position indices are moved up by the cigar len
        # For an insertion is relative to the reference.  Query index is incremented but reference index remains the
        # same
        # For a deletion, we increment the reference index but keep the query index the same

        for i in self.cigar_operations:
            current_query_start, current_query_end, current_reference_start, current_reference_end = self.increment_indices(
            i, current_query_start, current_reference_start)
            QR = SequenceRange(current_query_start, current_query_end, i)
            RR = SequenceRange(current_reference_start, current_reference_end, i)

            # assume the next range starts in the next base over.  If not, adjust in increment_indices
            current_reference_start = current_reference_end + 1
            current_query_start = current_query_end + 1
            self.query_ranges.append(QR)
            self.reference_ranges.append(RR)

        # adjust the coordinates if the mapping is on the reverse strand
        if not self.is_transcript_forward:
            # mapping is 3'->5'
            # mapping ranges stay the same, but the actual coordinates are modified

            current_pos = 0
            for i in reversed(self.query_ranges):
                range_len = i.stop_pos - i.start_pos
                if range_len == 0:
                    current_pos -= 1
                current_end = current_pos + range_len
                i.start_pos = current_pos
                i.stop_pos = current_end
                # assume the next range will start adjacent to the end
                current_pos = current_pos + range_len + 1

            # something amiss here if this is not true
            assert self.query_ranges[-1].start_pos == 0


    @staticmethod
    def get_pos(SR1, SR2, query_coordinate, is_forward_SR1, is_forward_SR2):
        """
    
        :param SR1: SequenceRange list in which the query_coordinate is location
        :param SR2: SequenceRange lis.  'Other' list in which the query coordinate is to be translated.
        :param query_coordinate: int representing the position to query
        :param is_forward_SR1: boolean.  Is the sequence range SR1 mapping 5'->3'
        :param is_forward_SR2: boolean, Is the sequence range SR2 mapping 5'->3'
        :return: tuple representing the translated coordinate  (min_pos,max_pos).  For a match, min_pos==max_pos.  
        For an insertion, the min_pos is the coordinate in SR2 that is immediately before the insertion and max_pos 
        is the coordinate in SR2 immediately after the insertion.  
        """

        if not is_forward_SR1 and not is_forward_SR2:
            sys.stderr.write("Query and reference mappings cannot both be on reverse strand. Skipping\n")
            return None

        if query_coordinate < 0:
            sys.stderr.write("Query position is negative.  Cannot process. Skipping\n")
            return None

        # get the last coordinate of SR1 to make sure that the queried coordinate is within range
        final_index = -1
        if not is_forward_SR1:
            final_index = 0

        if query_coordinate > SR1[final_index].stop_pos:
            sys.stderr.write("Requested position is greater than the length of the alignment of query on ref.\n")
            return None

        genomic_pos = None
        matching_positions = None

        # iterate over SR1 5'->3'
        range_indices = range(0, len(SR1))
        if not is_forward_SR1:
            range_indices = reversed(range_indices)

        for i in range_indices:
            sequence_range = SR1[i]

            # check if this range contains the query position
            if query_coordinate >= sequence_range.start_pos and query_coordinate <= sequence_range.stop_pos:
                # what's the offset of the position in the range
                offset = query_coordinate - sequence_range.start_pos

                # check that we've gotten the right offset and range
                t_pos = sequence_range.start_pos + offset
                assert t_pos == query_coordinate

                # now find the position in SR2.  The translated position will be in the same range as i
                if sequence_range.cigar_operation.operation == "M":
                    genomic_pos = SR2[i].start_pos + offset
                    if is_forward_SR1 != is_forward_SR2:
                        genomic_pos = SR2[i].stop_pos - offset
                    matching_positions = (genomic_pos, genomic_pos)
                else:
                    # Insertion. Find the coordinates immediately before and after insertion

                    prev_bin = i - 1  # assume that this bin is atleast 2 since we don't start with indel cigar ops
                    next_bin = i + 1
                    prev_coord = SR2[prev_bin].stop_pos
                    next_coord = SR2[next_bin].start_pos
                    if not is_forward_SR2:
                        next_coord = SR2[prev_bin].start_pos
                        prev_coord = SR2[next_bin].stop_pos
                    matching_positions = (prev_coord, next_coord)

                if matching_positions is not None:
                    break

        if matching_positions is None:
            sys.stderr.write("Could not locate position in query sequence\n")
        return matching_positions


    @staticmethod
    def transcript_to_genomic_pos(input_position, M):
        genomic_pos = M.get_pos(M.query_ranges, M.reference_ranges, input_position, M.is_transcript_forward, True)
        return genomic_pos


    @staticmethod
    def genomic_to_transcript_pos(input_position, M):
        transcript_pos = M.get_pos(M.reference_ranges, M.query_ranges, input_position, True, M.is_transcript_forward)
        return transcript_pos

class CigarOperation:
    """
    Class to hold information related to a one cigar operation
    :param: op_length, int.  The length of the cigar string
    :param: operation, string.  The type of cigar operation (MDI)
    """

    def __init__(self, op_length, operation):
        # assume that the input is well formed
        self.op_length = op_length
        self.operation = operation


class SequenceRange:
    """
     Class which represents one contiguous segment of sequence(query or reference) which matches the cigar string
    :param: start_pos, int.  The starting position of this range.  Start is lower than top, regardless of mapping 
    orientation
    :param: stop_pos, int.  The last position in this range. 
    :param: cigar_operation, CigarOperation object.  Represents the cigar operation which corresponds to this range. 
    """
    def __init__(self, start_pos, stop_pos, cigar_operation):
        assert stop_pos >= start_pos
        self.start_pos = start_pos
        self.stop_pos = stop_pos
        self.cigar_operation = cigar_operation

