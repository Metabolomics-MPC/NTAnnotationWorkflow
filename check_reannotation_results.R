# load required libraries
library(tools)
library(tidyverse)

# check if MS1 annotation results are the same ---------------------------------
# external databases, negative mode
md5sum("Demo/new/test_output/Annotation_MS1_external/neg_hmdb_sample_ms1annotation.tsv") == 
  md5sum("Demo/new/test_output_reanno/Annotation_MS1_external/neg_hmdb_sample_ms1annotation.tsv")

# external databases, positive mode
md5sum("Demo/new/test_output/Annotation_MS1_external/pos_hmdb_sample_ms1annotation.tsv") == 
  md5sum("Demo/new/test_output_reanno/Annotation_MS1_external/pos_hmdb_sample_ms1annotation.tsv")

# internal databases, negative mode
md5sum("Demo/new/test_output/Annotation_MS1_inhouse/neg_F5_ms1annotation.tsv") == 
  md5sum("Demo/new/test_output_reanno/Annotation_MS1_inhouse/neg_F5_ms1annotation.tsv")

# internal databases, positive mode
md5sum("Demo/new/test_output/Annotation_MS1_inhouse/pos_F5_ms1annotation.tsv") == 
  md5sum("Demo/new/test_output_reanno/Annotation_MS1_inhouse/pos_F5_ms1annotation.tsv")


# check if MS2 annotation results are the same ---------------------------------
# internal databases, negative mode
md5sum("Demo/new/test_output/Annotation_MS2_inhouse/neg_MassbankRecord_Neg_ms2annotation.tsv") ==
  md5sum("Demo/new/test_output_reanno/Annotation_MS2_inhouse/neg_MassbankRecord_Neg_ms2annotation.tsv")

# internal database, positive mode
md5sum("Demo/new/test_output/Annotation_MS2_inhouse/pos_MassbankRecord_Pos_ms2annotation.tsv") ==
  md5sum("Demo/new/test_output_reanno/Annotation_MS2_inhouse/pos_MassbankRecord_Pos_ms2annotation.tsv")

# external databases, negative mode
md5sum("Demo/new/test_output/Annotation_MS2_external/") ==
  md5sum("Demo/new/test_output_reanno/Annotation_MS2_inhouse/neg_MassbankRecord_Neg_ms2annotation.tsv")

