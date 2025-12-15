echo "Cleaning fastq reads in $1"

for i in $1/*R1*.fastq.gz;do
	R1="$i"
	R2=$(echo $R1 | sed -e 's/R1/R2/')
	SAMPLE=$(echo $R1 | cut -d_ -f1 | cut -d/ -f3)
	OUT1=${R1##*/}
	OUT2=${R2##*/}


	echo "Using: "
	echo $R1
	echo $R2
	echo $SAMPLE

	fastp --in1 $R1 --in2 $R2 --out1 ${OUT1%%.*}.clean.fastq.gz --out2 ${OUT2%%.*}.clean.fastq.gz --detect_adapter_for_pe -g -x -q 28 --length_required 150 -h $SAMPLE.html -j $SAMPLE.json -R $SAMPLE  --thread 16
done;

