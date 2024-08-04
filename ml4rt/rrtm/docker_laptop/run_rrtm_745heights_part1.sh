#!/usr/bin/bash

# Runs the full RRTM with 745 heights.
# Argument 1 is the string ID for the Docker container.
# Argument 2 is the sudo password.

SEPARATOR_STRING="\r\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\r\n"

DATE_STRINGS=("20191201" "20191202" "20191203" "20191204" "20191205" "20191206" "20191207" "20191208" "20191209" "20191210" "20191211" "20191212" "20191213" "20191214" "20191215" "20191216" "20191217" "20191218" "20191219" "20191220" "20191221" "20191222" "20191223" "20191224" "20191225" "20191226" "20191227" "20191228" "20191229" "20191230" "20191231" "20200101" "20200102" "20200103" "20200104" "20200105" "20200106" "20200107" "20200108" "20200109" "20200110" "20200111" "20200112" "20200113" "20200114" "20200115" "20200116" "20200117" "20200118" "20200119" "20200120" "20200121" "20200122" "20200123" "20200124" "20200125" "20200126" "20200127" "20200128" "20200129" "20200130" "20200131" "20200201" "20200202" "20200203" "20200204" "20200205" "20200206" "20200207" "20200208" "20200209" "20200210" "20200211" "20200212" "20200213" "20200214" "20200215" "20200216" "20200217" "20200218" "20200219" "20200220" "20200221" "20200222" "20200223" "20200224" "20200225" "20200226" "20200227" "20200228" "20200229" "20200301" "20200302" "20200303" "20200304" "20200305" "20200306" "20200307" "20200308" "20200309" "20200310" "20200311" "20200312" "20200313" "20200314" "20200315" "20200316" "20200317" "20200318" "20200319" "20200320" "20200321" "20200322" "20200323" "20200324" "20200325" "20200326" "20200327" "20200328" "20200329" "20200330" "20200331" "20200401" "20200402" "20200403" "20200404" "20200405" "20200406" "20200407" "20200408" "20200409" "20200410" "20200411" "20200412" "20200413" "20200414" "20200415" "20200416" "20200417" "20200418" "20200419" "20200420" "20200421" "20200422" "20200423" "20200424" "20200425" "20200426" "20200427" "20200428" "20200429" "20200430" "20200501" "20200502" "20200503" "20200504" "20200505" "20200506" "20200507" "20200508" "20200509" "20200510" "20200511" "20200512" "20200513" "20200514" "20200515" "20200516" "20200517" "20200518" "20200519")

container_id_string=$1
sudo_password=$2

for((i=0; i<22; i++)); do
    echo ${DATE_STRINGS[$i]}
    
    bash run_rrtm_one_day_745heights.sh "${container_id_string}" "${DATE_STRINGS[$i]}" "${sudo_password}"
    echo $SEPARATOR_STRING
done