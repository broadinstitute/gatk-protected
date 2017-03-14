#!/bin/bash -l
set -e
#cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

ln -s /home/travis/build/broadinstitute/gatk-protected/scripts/mutect2_wdl/mutect2.wdl
java -jar ~/cromwell-0.25.jar run /home/travis/build/broadinstitute/gatk-protected/scripts/mutect2_wdl/mutect2_multi_sample.wdl /home/travis/build/broadinstitute/test_m2_wdl_multi_mod.json - test_m2_wdl.metadata
