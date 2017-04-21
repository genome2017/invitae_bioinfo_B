Invitae Exercise B

Extended documentation available in doc/README.docx

**Usage**


python translate_coordinates.py  
--genome-mapping-file GENOME_MAPPING_FILE
                        REQUIRED. File specifying mappings (inputfile1.txt in exercise
                        specifications)
  --transcript-processing-file TRANSCRIPT_PROCESSING_FILE, REQUIRED
                        File specifying transcripts to process (inputfile2.txt
                        in exercise specifications)
  --output_file OUTPUT_FILE, OPTIONAL
                        Name of output file to write results to. Default is
                        output.txt if none provided.

**Unit Tests**

Unit tests on the method which performs coordinate translation is found in tests/test_class.py.  There are several flavors of tests here.  To add a new test, create a Mappings object (see examples in setup classes), and you can run tests in the ‘test’ method.

To run the tests, please enter:
nosetests test/test_class.py

