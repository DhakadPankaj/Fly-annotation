# /bin/bash

main() {
  # Check the number of arguments
  if (( $# != 1 )); then
    echo "Usage: cat_setup.sh cactus_config_file"
    exit 1
  fi

  local cactusfile=$1
  local clade=$1

  local bampath="/data/home/s2215768/fly_annotation/data/rnaseq_all/bams"
  local halpath="/data/home/s2215768/fly_annotation/cat_annot/hal_GenomeAlignments"
  local gffpath="/data/home/s2215768/fly_annotation/cat_annot/refseq_annotation"
  local proteinpath="/data/home/s2215768/fly_annotation/cat_annot/miniprot_hints"

  # Copy cactus_config file
  local root_anc=$(awk 'NR==3' "$clade" | sed 's/#Root: //')
  local hal_name="${clade%_cactus.config}.hal"
  local refsp=$(awk -F " " '{if($1=="#Reference_species:") print $2}' "$clade")
  local out_sp=`cat ${cactusfile}|awk -F " " '{if($1=="#Outgroup_species:") print $2}'`
  cat ${cactusfile}|tail -n +8|grep -v "^$"|cut -d " " -f1|sed 's/*//'|grep -v "${out_sp}" > ${cactusfile%_cactus.config}.tmp
  # Process BAM and INTRONBAM files
  local bam_count=0
  local intronbam_count=0
  local tmp_file="${cactusfile%_cactus.config}.tmp"
  while IFS= read -r species; do
    if [[ -s "${bampath}/${species}_rnaseq.bam" ]]; then
      bam_count=$((bam_count+1))
    elif [[ -s "${bampath}/${species}_rnaseq.intron.bam" ]]; then
      intronbam_count=$((intronbam_count+1))
    fi
  done < "$tmp_file"

  # Generate config file
  generate_config_file "$cactusfile" "$refsp" "$gffpath" "$bampath" "$proteinpath" "$tmp_file"

  # Create directories
  create_directories "$clade"

  # Get hints data if already exists
  copy_hints_data "$clade" "$tmp_file"

  # Run CAT commands
  run_cat_commands "$bam_count" "$intronbam_count" "$halpath" "$hal_name" "$refsp" "$cactusfile" "$clade" "$wd"

  # Cleanup
  cleanup "$clade" "$tmp_file"
}

generate_config_file() {
  local cactusfile=$1
  local refsp=$2
  local gffpath=$3
  local bampath=$4
  local proteinpath=$5
  local tmp_file=$6

  (
    echo -e "[ANNOTATION]"
    echo -e "${refsp} = ${gffpath}/${refsp}.cleaned.gff3"

    if [[ $bam_count -ge 1 ]]; then
      echo -e "\n[BAM]"
      while IFS= read -r species; do
        if [[ -s "${bampath}/${species}_rnaseq.bam" ]]; then
          echo -e "${species} = ${bampath}/${species}_rnaseq.bam"
        fi
      done < "$tmp_file"
    fi

    if [[ $intronbam_count -ge 1 ]]; then
      echo -e "\n[INTRONBAM]"
      while IFS= read -r species; do
        if [[ -s "${bampath}/${species}_rnaseq.intron.bam" ]]; then
          echo -e "${species} = ${bampath}/${species}_rnaseq.intron.bam"
        fi
      done < "$tmp_file"
    fi

    echo -e "\n[PROTEIN_FASTA]"
    while IFS= read -r species; do
      if [[ -s "${proteinpath}/${species}_miniprothints.gff" ]]; then
        echo -e "${species} = ${proteinpath}/${species}_miniprothints.gff"
      fi
    done < "$tmp_file"
  ) > "${cactusfile%_cactus.config}.cat.config"
}

create_directories() {
  local clade=$1

  mkdir -p "${clade%_cactus.config}/${clade%_cactus.config}_work_dir" "${clade%_cactus.config}/${clade%_cactus.config}_out_dir"
  mkdir -p "${clade%_cactus.config}/${clade%_cactus.config}_work_dir/hints_database"
}

copy_hints_data() {
  local clade=$1
  local tmp_file=$2

  mkdir -p "${clade%_cactus.config}/${clade%_cactus.config}_work_dir/hints_database"
  while IFS= read -r species; do
    if [[ -s "CAT_hintsdata/${species}.extrinsic_hints.gff" ]]; then
      cp "CAT_hintsdata/${species}.extrinsic_hints.gff" "${clade%_cactus.config}/${clade%_cactus.config}_work_dir/hints_database/"
    fi
  done < "$tmp_file"
}

run_cat_commands() {
  local bam_count=$1
  local intronbam_count=$2
  local halpath=$3
  local hal_name=$4
  local refsp=$5
  local cactusfile=$6
  local clade=$7
  local wd=$8

  luigid --background --logdir luigi_logs

  if [[ $bam_count -ge 1 ]] || [[ $intronbam_count -ge 1 ]]; then
    luigi --module cat RunCat --hal="${halpath}/${hal_name}" --ref-genome="$refsp" --config="${cactusfile%_cactus.config}.cat.config" \
    --work-dir "${clade%_cactus.config}/${clade%_cactus.config}_work_dir" \
    --workDir "${clade%_cactus.config}/toil_tmpDIR" --out-dir "${clade%_cactus.config}/${clade%_cactus.config}_out_dir" --augustus --local-scheduler \
    --augustus-species fly --augustus-cgp --cgp-param "${wd}/augustus_cfgs/cgp_parameters.cfg" \
    --augustus-cgp-cfg-template "${wd}/augustus_cfgs/cgp_extrinsic_template.cfg" \
    --tm-cfg "${wd}/augustus_cfgs/extrinsic.ETM1.cfg" \
    --tmr-cfg "${wd}/augustus_cfgs/extrinsic.ETM2.cfg" --assembly-hub --binary-mode local --workers=45 --maxCores=9
  else
    luigi --module cat RunCat --hal="${halpath}/${hal_name}" --ref-genome="$refsp" --config="${cactusfile%_cactus.config}.cat.config" \
    --work-dir "${clade%_cactus.config}/${clade%_cactus.config}_work_dir" \
    --workDir "${clade%_cactus.config}/toil_tmpDIR" --out-dir "${clade%_cactus.config}/${clade%_cactus.config}_out_dir" --augustus --local-scheduler \
    --augustus-species fly --assembly-hub --binary-mode local --workers=40 --maxCores=8
  fi
}

cleanup() {
  local clade=$1
  local tmp_file=$2

  mkdir -p CAT_workdir
  mkdir -p CAT_hintsdata

  if [[ -d "${clade%_cactus.config}/${clade%_cactus.config}_out_dir/consensus_gene_set" ]]; then
    cp "${clade%_cactus.config}/${clade%_cactus.config}_work_dir/hints_database/"*extrinsic_hints.gff CAT_hintsdata/
    tar cf - "${clade%_cactus.config}/${clade%_cactus.config}_work_dir" | pigz -p40 > "${clade%_cactus.config}_work_dir.tar.gz"
    mv "${clade%_cactus.config}_work_dir.tar.gz" CAT_workdir/
    rm -r "${clade%_cactus.config}/toil_tmpDIR" "${clade%_cactus.config}/${clade%_cactus.config}_work_dir"
  fi

  rm -r "$tmp_file"
}

main "$@"
