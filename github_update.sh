#!/bin/bash

cd /data/home/s2215768/fly_annotation/scripts
git add .
git commit -m "$1"
PASS="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"

xyz=$(expect -c "
spawn git push origin main 
expect \"Username for 'https://github.com':\"
send \"pd16\r\"
expect \"Password for 'https://pd16@github.com':\"
send \"$PASS\r\"
")

echo "$xyz"
echo "github repository updated!"

cd $(pwd)
