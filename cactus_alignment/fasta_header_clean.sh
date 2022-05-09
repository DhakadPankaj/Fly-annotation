#! /bin/bash

fa_name_cleaner() {
    filename=$1
    cat ${filename} | sed -E 's/(>[0-9A-Z\._]+).+/\1/' > $(echo ${filename} | sed -E 's/(\.f.a)/\.renamed\1/')
    cat ${filename} | grep ">" | sed -E 's/>([0-9A-Z\._]+)(.+)/\1\t\2/' > $(echo ${filename} | sed -E 's/\.(f.a)/\.\1info/')
} 
export -f fa_name_cleaner


fnafile=$1

fa_name_cleaner ${fnafile}

#mv ${fnafile%.fna}.renamed.fna ${fnafile}
