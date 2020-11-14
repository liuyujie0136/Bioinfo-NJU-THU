#!/bin/bash

### Bioinformatics Shell Scripting - scripts in book Bioinformatics Data Skills (chapter 12)

## basic bash scripts
results_dir="results/"
echo $results_dir

sample="CNTRL01A"
mkdir ${sample}_aln/ # $sample_aln/ may fail, because bash may think "sample_aln" as a whole
mkdir "${sample}_aln/" # more robust! recommended!

# Bash handles command-line arguments (which are assigned to the value $1, $2, $3, etc.). The variable $0 stores the name of the script.
bash args.sh arg1 arg2 arg3

if [ "$#" -lt 3 ] # $#: argument numbers
then
    echo "error: too few arguments, you provided $#, 3 required"
    echo "usage: script.sh arg1 arg2 arg3"
    exit 1
fi

test "ATG" = "atg" ; echo "$?" # echo "$?" to print the exit status

# -a as logical AND, -o as logical OR, ! as negation
if [ "$#" -ne 1 -o ! -r "$1" ]
then
    echo "usage: script.sh file_in.txt"
    exit 1
fi

# Bash arrays (0-indexed)
sample_names=(zmaysA zmaysB zmaysC)
echo ${sample_names[0]} # note: {} is used
echo ${sample_names[@]} # extract all elements

# use "basename" to strip(remove) a suffix - useful in creating output filenames
# usage: basename NAME [SUFFIX]; or: basename OPTION... NAME...
# Print NAME with any leading directory components removed. If specified, also remove a trailing SUFFIX.
basename seqs/zmaysA_R1.fastq .fastq
basename -s .fastq seqs/zmaysA_R1.fastq # these two are equal

# directly loop over files
#!/bin/bash
set -e
set -u
set -o pipefail

for file in *.fastq
do
    echo "$file: " $(bioawk -c fastx 'END {print NR}' $file) # count how many FASTQ records are in a file
done


## Automating File-Processing with find and xargs
# Finding Files with find
find zmays-snps # running find on a directory (without other arguments) can be a quick way to see a project directory’s structure

# basic syntax: find <path> <expression>
find data/seqs -name "zmaysB*fastq"
find data/seqs -name "zmaysB*fastq" -type f # only files; syntax: f for files, d for directories, and l for links
find data/seqs -name "zmaysB*fastq" -and -type f # equivalent to no-and version
find data/seqs -name "zmaysA*fastq" -or -name "zmaysC*fastq" # two patterns
find seqs -type f "!" -name "zmaysC*fastq" -and "!" -name "*-temp*" # not and not

# Common find expressions:
-name <pattern>; -type <x>
-iname <pattern> # case insensitive
-empty # empty files or directories
-size <filesize> # match file size
-regex # regular expression (extended: -E)
-print0 # separate results with null byte, not newline
-and; -or; -not; "!" # and, or, not, not

# find’s -exec: Running Commands on find’s Results
# passing each file that matches find’s expressions to the command specified by -exec
touch zmays{A,C}_R{1,2}-temp.fastq # create temp files
find . -name "*-temp.fastq" -exec rm -i {} \; # path ".": this dir; end ";": must be added to indicate the end of command; rm -i: prompt before every removal

# xargs (similar to R’s do.call())
# xargs works by taking input from standard in and splitting it by spaces, tabs, and newlines into arguments. Then, these arguments are passed to the command supplied.
find . -name "*-temp.fastq" | xargs rm # same as fine .. -exec rm

# using the null byte delimiter
find . -name "samples [AB].txt" -print0 | xargs -0 rm

# set how many arguments are passed to "each" command call (here is one "rm")
# important!
find . -name "*-temp.fastq" | xargs -n 1 rm

# write to a simple Bash script
find . -name "*-temp.fastq" | xargs -n 1 echo "rm -i" > delete-temp.sh # simply places arguments after the command you provide

# Using xargs with Replacement Strings to Apply Commands to Files
find . -name "*.fastq" | xargs basename -s ".fastq"
find . -name "*.fastq" | xargs basename -s ".fastq" | xargs -I{} fastq_stat --in {}.fastq --out ../summaries/{}.txt # -I{}: more fine-grained placement of arguments, replace {} with args; OR: --replace={}; fastq_stat: a foo program or shell script

# xargs and Parallelization
xargs -P <num> # <num> is the number of processes to run simultaneously

# xargs can run shell script, pragrams containing "|" or ">" must use this
find . -name "*.fastq" | xargs -n 1 -P 4 bash myscript.sh # -n 1 forces xargs to process one input argument at a time

# GNU Parallel
find . -name "*.fastq" | parallel --max-procs=6 '<program> {/.} > {/.}-out.txt' # {/.}: extract base filename without basename


## Make and Makefiles (skipped)
