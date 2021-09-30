load "pipeline-config.groovy"

align_pbmm2 = {
    transform("bam") {
        exec """
            ~/tools/miniconda3/envs/pbmm2/bin/pbmm2 align --preset HIFI
                $REF $input.fasta $output.bam
        """
    }
}

bam_to_cram = {
    transform("cram") {
        exec """
            samtools sort -O CRAM --reference $REF 
                -o $output.cram -T $output.cram.prefix
                $input.bam
        """
    }
}

@preserve("*.crai")
index_cram = {
    output.dir=file(input.cram).absoluteFile.parentFile.absolutePath
    transform("cram") to("cram.crai") {
        exec "samtools index $input.cram", "local"
    }
    forward input
}

trf = {
    from('fasta') transform('trf.dat') {
        exec """
            ~/tools/trf409.linux64 $input.fasta 2 7 7 80 10 50 6 -d -h -ngs > $output.dat
        """
    }
}

parse_trf = {
    from('dat', 'bam') transform('bed') {
        exec """
            $python $longstr/trf.py
                --dat $input.dat
                --out $output.bed
                $input.bam
        """
    }
}

cov = {
//    doc = "Get PacBio contigs in bed format - exclude secondary and supplementary alignments"
    from('cram') transform('cov') {
        exec """
            samtools view -h  -F 2304 $input.cram | bedtools bamtobed -i - > $output.cov
        """
    }
}

// Call variants with bcftools mpileup and filter to insertions great than 10bp
call_variants = {
    from('cram') transform('vcf') {
        exec """
            bcftools mpileup --min-ireads 1 -f $REF $input.cram | bioawk -Hc vcf '(length($alt)-length($ref))>=10' > $output.vcf
        """
    }
}

// Report most frequent kmer in all insertions >= 10 bp in a vcf
count_kmers = {
    from('vcf') transform('bed') {
        exec """
            $python $longstr/kinsertions.py $input.vcf $output.bed
        """
    }
}
compare = {
    from ('-genotype.txt', 'trf.bed') produce('.strling_pacbio.bed') {
        exec """
            $python $longstr/compare.py --strling $input.txt --trf $input.bed --cov $input.cov --out $output.bed
        """
    }
}


//run {
//    "%.fasta" * [align_pbmm2 + bam_to_cram + index_cram +
//        trf + parse_trf + cov ]//+ compare]
//}

run {
    "%.fasta" * [align_pbmm2 + bam_to_cram + index_cram +
                cov + call_variants + count_kmers]//+ compare]
}
