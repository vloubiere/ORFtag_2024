#!/bin/bash
set -e

function display_help {
    echo "Usage: $0 -l <layout> -b <BC12_regexpr> -B <BC14_regexpr> -o <output_prefix> -u <UMI_flag>"
    echo "Options:"
    echo "  -l <layout>            Layout type: 'PE' for paired-end or 'SE' for single-end"
    echo "  -b <BC12_regexpr>      Regular expression for BC12"
    echo "  -B <BC14_regexpr>      Regular expression for BC14"
    echo "  -o <output_prefix>     Output file prefix"
    echo "  -u <UMI_flag>          UMI flag: 'TRUE' or 'FALSE'"
    exit 1
}

layout=""
BC12_regexpr=""
BC14_regexpr=""
output_prefix=""
UMI_flag=""

while getopts "l:b:B:o:u:" opt; do
  case $opt in
    l)
      layout="$OPTARG"
      ;;
    b)
      BC12_regexpr="$OPTARG"
      ;;
    B)
      BC14_regexpr="$OPTARG"
      ;;
    o)
      output_prefix="$OPTARG"
      ;;
    u)
      UMI_flag="$OPTARG"
      ;;
    \?)
      display_help
      ;;
  esac
done

if [ -z "$layout" ] || [ -z "$BC12_regexpr" ] || [ -z "$BC14_regexpr" ] || [ -z "$output_prefix" ] || [ -z "$UMI_flag" ]; then
    display_help
fi

if [ "$UMI_flag" == "TRUE" ]; then
    if [ -n "$BC12_regexpr" ] && [ -n "$BC14_regexpr" ]; then
        echo "Error: When UMI_flag is set to TRUE, only one of BC12 or BC14 can be specified. The other one will be appended to read name."
        display_help
    fi
fi

output1="${output_prefix}_1.fq.gz"
output2="${output_prefix}_2.fq.gz"
output_single="${output_prefix}_single.fq.gz"

if [ -e "$output1" ]; then
    rm "$output1"
fi
if [ -e "$output2" ]; then
    rm "$output2"
fi
if [ -e "$output_single" ]; then
    rm "$output_single"
fi

function process_reads {
    local mate=0
    local line_count=0
    declare -A BC12_counts
    declare -A BC14_counts

    while read -r line; do
        IFS=$'\t' read -ra cur <<< "$line"
        ((line_count++))

        if [ "$UMI_flag" == "TRUE" ]; then
            if [ -n "$BC12_regexpr" ]; then
                BC="${cur[11]:5:${#cur[11]}}"
            elif [ -n "$BC14_regexpr" ]; then
                BC="${cur[13]:5:${#cur[13]}}"
            fi
            read_id="${cur[0]}_${BC}"
        else
            read_id="${cur[0]}"
        fi

        if [ "$layout" == "PE" ]; then
            if [ "$mate" -eq 0 ]; then
                if [[ "${cur[11]}" =~ $BC12_regexpr && "${cur[13]}" =~ $BC14_regexpr ]]; then
                    echo "@${read_id}"
                    echo "${cur[9]}"
                    echo "+"
                    echo "${cur[10]}" | gzip -c >> "$output1"
                    mate=1
                fi
            else
                echo "@${read_id}"
                echo "${cur[9]}"
                echo "+"
                echo "${cur[10]}" | gzip -c >> "$output2"
                mate=0
            fi
        elif [ "$layout" == "SE" ]; then
            if [[ "${cur[11]}" =~ $BC12_regexpr || "${cur[13]}" =~ $BC14_regexpr ]]; then
                echo "@${read_id}"
                echo "${cur[9]}"
                echo "+"
                echo "${cur[10]}" | gzip -c >> "$output_single"
            fi
        fi

        # Count the occurrences of BC12 and BC14 in the first 1000 lines
        if [ "$line_count" -le 1000 ]; then
            BC12_counts["${cur[11]:5:${#cur[11]}}"]=$((BC12_counts["${cur[11]:5:${#cur[11]}}"] + 1))
            BC14_counts["${cur[13]:5:${#cur[13]}}"]=$((BC14_counts["${cur[13]:5:${#cur[13]}}"] + 1))
        fi
    done

    # Display the top 10 values for BC12 and BC14
    echo "Top 10 values for BC12:"
    printf "%s\n" "${!BC12_counts[@]}" | sort -rn -k2 | head -n 10
    echo "Top 10 values for BC14:"
    printf "%s\n" "${!BC14_counts[@]}" | sort -rn -k2 | head -n 10
}

process_reads
