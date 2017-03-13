#!/bin/bash -l
set -e
#cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

java -jar ~/cromwell-0.25.jar run /home/travis/build/broadinstitute/gatk-protected/scripts/mutect2_wdl/mutect2_multi_sample.wdl pon_gatk_workflow_exome.json
