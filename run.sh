#!/bin/bash -e

dumpname=
dumpana_flags=
elements=
logfile=">&1"
name=
path="$(pwd)"

# Text to be shown when help is requested
show_help() {
cat << EOF

Usage: ${0##*/} --scripts=<name> --dump=<dumpname> --elements="<elements>"
    [-h | --help] [--logfile=<path>] [--path=<path>] [-- <extra dumpana options>]

Run dumpana binary multiple times with <dumpname> file loaded and <elements> mapped, and commands
parsed from files with name <name>*.commands (in lexographical order).

**IMPORTANT NOTICE**
<name>*.commands files must have all commands listed in separate new line each. If multiple values
are being read from single standard input buffer flush, then list them space-separeted.

    -h, --help                          display help & exit
    ________________________________________________________________________________________________
    mandatory:

    --dump=<dumpname>                   dump file to be loaded name/names
                                        e.g. (--dump=frames1.dump, --dump="frames1.dump frames2.dump")
    --elements="<elements>""            define elements to map atomic types to; must be valid element
                                        chemical symbol
                                        e.g. (--elements="Zr Cu", --elements="Na Cl")
    --scripts=<name>                    scripts with commands base filename (this needs to be equivalent
                                        to key for 'find . -name <name>*.commands')
                                        e.g. (--scripts=auto: will execute all commands in files
                                        auto1.commands, auto2.commands, auto3.commands, etc.
                                        --scripts=automatic_scripts: will execute all comands in files
                                        automatic_scripts_a.commands, automatic_scripts_b.commands, etc.)
    ________________________________________________________________________________________________
    optional:

    --logfile=<path>                    path to file to which output from dumpana binary will be redirected;
                                        script will add time signature after '_' character
                                        e.g. (--logfile=log1, --logfile=userslog)
    --path=<path>                       path to dumpana binary; use when script is run not from the
                                        directory containing binary
                                        e.g. (--path=/path/to/binary)
    -- <extra dumpana options>          pass extra paremeters to dumpana binary
                                        e.g. (-- -1, -- -os -mm)

EOF
}

# trap function for all exit signals
on_exit() {
    rv=$?
    if [[ $rv -eq 0 ]]; then
        echo -e "\e[32mFinished successfully"
    else
        echo -e "\e[1m\e[91mFailed with code: $rv"
    fi
    exit $rv
}

# options parsing
OPTS=`getopt -n "dumpana-script-runner" -a -o h -l dump:,elements:,help,logfile:,path:,scripts: -- "$@"`
eval set -- "$OPTS"
while :; do
    case "$1" in
        --dump)
            dumpname="$2"
            shift 2 ;;
        --elements)
            elements="$2"
            shift 2 ;;
        --logfile)
            logfile="$2"
            shift 2 ;;
        --path)
            path="$2"
            shift 2 ;;
        --scripts)
            name="$2"
            shift 2 ;;
        --)
            shift
            dumpana_flags="$*"
            if ! [[ -z $dumpname || -z $elements || -z $name ]]; then
                echo "Dump file(s) to load: \"$dumpname\""
                echo "Matched commands files: \"$name*.commands\""
                if [[ ! -z $dumpana_flags ]]; then
                    echo "Extra dumpana flags: \"$dumpana_flags\""
                fi
                echo "Mapped elements: \"$elements\""
                if [[ $logfile != ">&1" ]]; then
                    echo "Logfile with output: $logfile"
                fi
                break
            fi ;;
        --h|--help|*)
            show_help >&2
            exit 1 ;;
    esac
done

# set trap on exit signals
trap "on_exit" EXIT

# parse relative path
if [[ ${path:0:2} == "./" ]]; then
	path="$(pwd)/${path:2}"
fi
# change possible "//" to "/"
path=${path//\/\///}

#check if binary exists
if ! [[ -e $path/dumpana ]]; then
	echo "Could not find dumpana binary file in path $path"
	exit 1
fi

# reset single commands file
if [[ -e tmp.commands ]]; then
    rm tmp.commands
fi

# feed the elements to commands file
echo "$elements" > tmp.commands

# find all scripts with commands
# execute find & split it's output by newline delimiter
IFS=$'\n' cmds_array=(`find . -path "*$name*.commands" | sort`)
echo "Found ${#cmds_array[@]} commands file(s)."
for file in "${cmds_array[@]}"; do
    cat "$file" >> tmp.commands
    echo -e "" >> tmp.commands
done

# clear all empty lines
sed -i '/^$/d' tmp.commands
# push exit command to file
echo -e "0\n" >> tmp.commands

if [[ $logfile != ">&1" ]]; then
    rm -rf "$logfile"
    echo "Dump file(s) loaded: \"$dumpname\"" &>> $logfile
    echo "Matched ${#cmds_array[@]} commands files: \"$name*.commands\"" &>> $logfile
    if [[ ! -z $dumpana_flags ]]; then
        echo "Extra dumpana flags: \"$dumpana_flags\"" &>> $logfile
    fi
    echo "Mapped elements: \"$elements\"" &>> $logfile
    logfile="&>> $logfile"
fi

# call program & redirect output to log file
cmd="$path/dumpana $dumpana_flags $dumpname < tmp.commands $logfile"
echo "Command called: $cmd" | tr -s ' '
eval $cmd

# cleanup
rm tmp.commands
