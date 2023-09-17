# Source: https://github.com/biomonika/HPP/blob/main/assembly/wdl/workflows/assessAsemblyCompletness.wdl
# Download date: 17 Sept 2023

version 1.0


workflow assessAssemblyCompletness{
    input {
        File assembly
        String assembly_name=basename(sub(sub(sub(assembly, "\\.gz$", ""), "\\.fasta$", ""), "\\.fa$", "")) #remove the file extension
        File reference
        Int threadCount=1
        Int preemptible = 1
    }

    call formatAssembly {
        input:
            assembly=assembly,
            assembly_name=assembly_name,
            preemptible=preemptible
    }

    call generateAssemblyEdges {
        input:
            assembly=formatAssembly.formattedAssembly,
            assembly_name=assembly_name,
            preemptible=preemptible
    }

    call runBioawk {
        input:
            assembly=assembly,
            assembly_name=assembly_name,
            edges=generateAssemblyEdges.edges,
            preemptible=preemptible
    }
    
    call combineIntoFasta {
        input:
            assembly_name=assembly_name,
            headers=runBioawk.headers,
            edges=generateAssemblyEdges.edges,
            preemptible=preemptible
    }
    
    call runNCRF {
        input:
            edgesFasta=combineIntoFasta.edgesFasta,
            assembly_name=assembly_name,
            preemptible=preemptible
    }


    call runMashMap {
        input:
            assembly=assembly,
            assembly_name=assembly_name,
            reference=reference,
            threadCount=threadCount,
            preemptible=preemptible
    }

    call assessCompletness {
        input:
            assembly_name=assembly_name,
            telomericEnds=runNCRF.ends,
            lengths=runBioawk.lengths,
            unknown=runBioawk.unknown,
            mashmap=runMashMap.mashmap,
            preemptible=preemptible
    }

    output {
        File T2Tcontigs = assessCompletness.contigs
        File T2Tscaffolds = assessCompletness.scaffolds
    }

    parameter_meta {
        assembly: " Assembly for which to assess completness in fasta format"
    }
    meta {
        author: "Monika Cechova"
        email: "mcechova@ucsc.edu"
    }
}

task formatAssembly {
    input{
        File assembly
        String assembly_name
        Int memSizeGB = 32
        Int preemptible
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        seqtk seq -l0 ~{assembly} >~{assembly_name}.formatted.fa

        
    >>>

    output {
        File formattedAssembly="~{assembly_name}.formatted.fa"
    }

    runtime {
    }
}

task generateAssemblyEdges {
    input{
        File assembly
        String assembly_name
        Int memSizeGB = 32
        Int preemptible
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        bases=1000; 
        #extract first and last bases from each fasta entry
        
        #extract bases from the start and end of each entry in fasta
        cat ~{assembly} | gawk -v len="$bases" -F "" '{if (NR % 2 == 0) {for (i=1; i<=NF; i++) {if (i <= len || (NF-i) < len) {printf $(i)}}; printf "\n"}}' >~{assembly_name}.edges.tmp.txt
        
        #split start and end into two rows
        cat ~{assembly_name}.edges.tmp.txt | sed -e "s/.\{$bases\}/&\n/g" | sed '/^$/d' >~{assembly_name}.edges.txt
        #rm ~{assembly_name}.edges.tmp.txt

        
    >>>

    output {
        File edges="${assembly_name}.edges.txt"
    }

    runtime {
    }
}

task runBioawk {
    input{
        File assembly
        File edges
        String assembly_name
        Int memSizeGB = 32
        Int preemptible

    }

    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        #create file with the number of unknown bases
        bioawk -c fastx '{n=gsub(/N/, "", $seq); print $name "\t" n}' ~{assembly} >~{assembly_name}.unknown.txt

        #create headers with the contig names and annotating starts and ends
        bioawk -c fastx '{print ">"$name"_start";print ">"$name"_end";}' ~{assembly} >~{assembly_name}.headers.txt

        #create file with sequence lengths
        bioawk -c fastx '{ print $name, length($seq) }' < ~{assembly} >~{assembly_name}.lengths.txt


        
    >>>

    output {
        File unknown = "${assembly_name}.unknown.txt"
        File headers = "${assembly_name}.headers.txt"
        File lengths = "${assembly_name}.lengths.txt"
    }

    runtime {
    }
}

task combineIntoFasta {
    input{
        File headers
        File edges
        String assembly_name
        Int preemptible
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        #combine extracted sequences with headers to create .fasta file
        paste -d '\n' ~{headers} ~{edges} > ~{assembly_name}.combined.fasta
        
    >>>

    output {
        File edgesFasta = "${assembly_name}.combined.fasta"
    }

    runtime {
    }

}

task runNCRF {
    input{
        File edgesFasta
        String assembly_name
        Int memSizeGB = 32
        Int preemptible
    }
    command <<<

    #handle potential errors and quit early
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    #search this fasta for the telomeric repeat
    cat ~{edgesFasta} | NCRF --stats=events TTAGGG | ncrf_summary >~{assembly_name}.telomeric.summary.txt
    cat ~{assembly_name}.telomeric.summary.txt | tail -n +2 | cut -f3 | sort | uniq >~{assembly_name}.telomeric.ends.txt
        
    >>>

    output {
        File ends = "${assembly_name}.telomeric.ends.txt"
    }

    runtime {
    }

}

task runMashMap {
    input{
        File assembly
        String assembly_name
        File reference
        Int memSizeGB = 32
        Int threadCount
        Int preemptible
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mashmap --threads ~{threadCount} --perc_identity 95 --noSplit -r ~{reference} -q ~{assembly} -o ~{assembly_name}.mashmap.txt
        
    >>>

    output {
        File mashmap = "${assembly_name}.mashmap.txt"
    }

    runtime {
    }

}

task assessCompletness {
    input{
        String assembly_name
        File telomericEnds
        File lengths
        File unknown
        File mashmap
        Int memSizeGB = 32
        Int preemptible
    }
    command <<<

    #handle potential errors and quit early
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    #CREATE A SUMMARY FILE
    cat ~{telomericEnds} | while read line; do contig=`echo $line | sed 's/_.*//'`; egrep --color ${contig} ~{lengths} | tr -d '\n'; echo -ne '\t'; egrep --color ${contig} ~{unknown} | tr -d '\n'; echo -ne '\t'; egrep --color ${contig} ~{mashmap} | tr -s "[:blank:]" "\t" ; done | cut -f1,2,4,10 | sort | uniq -c | sort -rgk1 | tr -s "[:blank:]" "\t" | sed 's/^\t//' >~{assembly_name}.SUMMARY.txt

    #write headers
    echo -e "name\tlength\tNs\tchromosome" >~{assembly_name}.T2T.scaffolds.txt
    echo -e "name\tlength\tNs\tchromosome" >~{assembly_name}.T2T.contigs.txt

    #check if there are telomeres on both sides
    if grep -q "^2" ~{assembly_name}.SUMMARY.txt; then
        #write T2T scaffolds, the number of Ns must be 0
        cat ~{assembly_name}.SUMMARY.txt | grep "^2" | awk '{if ($4>=0) print;}'  | cut -f2- | sort -k4,4 -V -s >>~{assembly_name}.T2T.scaffolds.txt || true

        #write T2T contigs, the number of Ns is greater than 0
        cat ~{assembly_name}.SUMMARY.txt | grep "^2" | awk '{if ($4==0) print;}'  | cut -f2- | sort -k4,4 -V -s >>~{assembly_name}.T2T.contigs.txt || true
    else
        echo "None of the contigs/scaffolds contained telomeres on both ends."
    fi
    
    echo "Done."    
    >>>

    output {
        File scaffolds = "${assembly_name}.T2T.scaffolds.txt"
        File contigs = "${assembly_name}.T2T.contigs.txt"
    }

    runtime {
    }
}

