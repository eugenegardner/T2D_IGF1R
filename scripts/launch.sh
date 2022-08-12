#!/bin/bash

TAR=$1
PHENO=$2
BIN=$3
OUT=$4
SEX=$5
COVARFILE=$6
COVARNAMES=$7

dx run mrcepid-runassociationtesting --brief --priority normal --destination results/ -iassociation_tarballs=$TAR -iphenofile=$PHENO -iis_binary=$BIN -ioutput_prefix="${OUT}.bolt" -isex=$SEX -itool=bolt -imode=burden -iinclusion_list=file-G6qXvjjJ2vfQGPp04ZGf6ygj --yes --name "${OUT}.bolt"

dx run mrcepid-runassociationtesting --brief --priority normal --destination results/ -iassociation_tarballs=$TAR -iphenofile=$PHENO -iis_binary=$BIN -ioutput_prefix="${OUT}.staar" -isex=$SEX -itool=staar -imode=burden -iinclusion_list=file-G6qXvjjJ2vfQGPp04ZGf6ygj --yes --name "${OUT}.staar"

dx run mrcepid-runassociationtesting --brief --priority normal --destination results/ -iassociation_tarballs=$TAR -iphenofile=$PHENO -iis_binary=$BIN -ioutput_prefix="${OUT}.glm" -isex=$SEX -itool=glm -imode=burden -iinclusion_list=file-G6qXvjjJ2vfQGPp04ZGf6ygj --yes --name "${OUT}.glm"






