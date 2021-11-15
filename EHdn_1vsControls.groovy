EHdn_path="~/tools/ExpansionHunterDenovo-v0.9.0-linux_x86_64"
PYTHON="/uufs/chpc.utah.edu/common/HIPAA/u6026198/tools/miniconda3/envs/STR/bin/python"
REF="/uufs/chpc.utah.edu/common/HIPAA/u6026198/storage/ref-data/GATK_Bundle_Build38/Homo_sapiens_assembly38.fasta"
//disease_loci="~/storage/ref-data/hg38.STR_disease_loci.2000slop.bed"
disease_loci="~/storage/ref-data/hg38.STR_disease_loci.500slop.bed"

//load "EHdn.groovy"
load "STR_tp_samples_config.groovy"
// control_samples
REF="/uufs/chpc.utah.edu/common/HIPAA/u6026198/storage/ref-data/GATK_Bundle_Build38/Homo_sapiens_assembly38.fasta"

def get_sample(filename) {
    return(filename.split("/")[-1].split("\\.")[0])
}

if(args.any { it.endsWith('.cram') })
    input_type = 'cram'
else
    input_type='bam'

set_test_sample = {
    branch.test_sample = get_sample(branch.name)
    println "test sample: "+branch.test_sample
}

profile = {
    //output.dir = "profiles"
    branch.sample = get_sample(branch.name)
    println "profile sample: "+branch.sample

    produce(branch.sample+".str_profile.json") {
        exec """
            $EHdn_path/bin/ExpansionHunterDenovo profile
                --reads $input
                --reference $REF
                --output-prefix $output.prefix.prefix
        """
    }
}

manifest_case = {
    //output.dir = "profiles"
    produce(branch.test_sample+".manifest") {
        exec """
            echo -e "${branch.test_sample}\tcase\t${branch.test_sample}.str_profile.json" > $output.manifest
        """
    }
}

manifest_control = {
    //output.dir = "profiles"
    produce(branch.sample+".manifest") {
        exec """
            echo -e "${branch.sample}\tcontrol\t$input.json" > $output.manifest
        """
    }
}

manifest_controls = {
    from("*.manifest") produce("control.manifest") {
        exec """
            cat $inputs.manifest > $output.manifest
        """
    }
}

manifest_all = {
    //output.dir = "profiles"
    from("*.manifest") produce(branch.test_sample + ".all.manifest") {
        exec """
            cat $inputs.manifest control.manifest > $output.manifest
        """
    }
}

merge = {
    //output.dir = "profiles"
    produce(branch.test_sample + ".all.multisample_profile.json") {
        exec """
            $EHdn_path/bin/ExpansionHunterDenovo merge
                --reference $REF
                --manifest $input.manifest
                --output-prefix ${branch.test_sample}.all
        """
    }
}

outlier = {
    produce(branch.test_sample + ".all.outlier_locus.tsv") {
        exec """
            $PYTHON $EHdn_path/scripts/outlier.py locus
                --manifest $input.manifest
                --multisample-profile $input.json
                --output $output.tsv
        """
    }
}

sort_tsv = {
    filter("sorted") {
        exec """
            bash ~/storage/runEHdn/sort.sh $input.tsv > $output.tsv
        """
    }
}

// Filter to pathogenic loci and also limit to rows matching the test sample
path_loci = {
    filter("path_loci") {
        exec """
            bedtools intersect -wo -a $input.tsv -b $disease_loci | grep ${branch.test_sample} > $output.tsv
        """
    }
}

// Filter to pathogenic loci and also limit to rows matching the test sample
filter_path = {
    filter("filtered") {
        exec """
            ~/tools/miniconda3/bin/python ~/storage/git/STRling/scripts/match-repeatunit.py EHdn $input.tsv $output.tsv
        """
    }
}

run {
    control_samples * [profile + manifest_control] + manifest_controls +
    "%.${input_type}" * [profile] +
    "%.${input_type}" * [set_test_sample + manifest_case +
        manifest_all + merge +
        outlier +
        sort_tsv + path_loci + filter_path]
}

