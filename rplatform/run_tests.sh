#/bin/sh

echo "Executing $0"
echo "Environment: ${rp_env}"
echo "Working directory: `pwd`"
echo "Working directory contains: `ls | tr '\n' ' '`"

# Define how you run tests below
# exit when any command fails
set -e

echo ">>>>> RUNNING UNIT TESTS"
Rscript -e "devtools::test(pkg = '/mnt/vol/MiDAS', stop_on_failure = TRUE)"
