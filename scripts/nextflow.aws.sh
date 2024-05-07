#!/bin/bash

NEXTFLOW_PROJECT=$1  # first argument is a project
shift
NEXTFLOW_PARAMS="$@"  # additional arguments are nextflow options (e.g. -resume) or workflow parameters

# Create the default config using environment variables
# passed into the container by AWS Batch
NF_CONFIG=~/.nextflow/config

cat << EOF > $NF_CONFIG
workDir = "$NF_WORKDIR"
process.executor = "awsbatch"
process.queue = "$NF_JOB_QUEUE"
aws.batch.cliPath = "/home/ec2-user/miniconda/bin/aws"
EOF

echo "=== CONFIGURATION ==="
cat ~/.nextflow/config

# stage in session cache
# .nextflow directory holds all session information for the current and past runs.
# it should be `sync`'d with an s3 uri, so that runs from previous sessions can be
# resumed
echo "== Restoring Session Cache =="
aws s3 sync --only-show-errors $NF_LOGSDIR/.nextflow .nextflow

function preserve_session() {
    # stage out session cache
    if [ -d .nextflow ]; then
        echo "== Preserving Session Cache =="
        aws s3 sync --only-show-errors .nextflow $NF_LOGSDIR/.nextflow
    fi

    # .nextflow.log file has more detailed logging from the workflow session and is
    # nominally unique per session.
    #
    # when run locally, .nextflow.logs are automatically rotated
    # when syncing to S3 uniquely identify logs by the AWS Batch Job ID
    if [ -f .nextflow.log ]; then
        echo "== Preserving Session Log =="
        aws s3 cp --only-show-errors .nextflow.log $NF_LOGSDIR/.nextflow.log.$AWS_BATCH_JOB_ID
    fi
}

# set a trap so that session information is preserved when the container exits
trap preserve_session EXIT

echo "== Running Workflow =="
echo "nextflow run $NEXTFLOW_PROJECT $NEXTFLOW_PARAMS"
nextflow run $NEXTFLOW_PROJECT $NEXTFLOW_PARAMS
