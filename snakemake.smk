input = open('samples.txt', 'r')
samples = input.read()
samples = samples.split('\n')
samples = samples[:-1]

rule all:
    input:
        expand('06_MG_PROFILES/{sample}/PROFILE.db', sample = samples),
        expand('05_CONTIGS_DB/{sample}/contigs.db', sample = samples)

rule simplifyContigs:
    input:
        fasta='02_ASSEMBLY/{sample}/final.contigs.fa'
    output:
        reformattted = '02_ASSEMBLY/{sample}/contigs.simplified.fa',
        report='02_ASSEMBLY/{sample}/contig_Rename_Report-txt'
    log:
        err ='02_ASSEMBLY/{sample}/reformat.err',
        out='02_ASSEMBLY/{sample}/reformat.out'
    resources:
        cpus_per_task=1,
        mem_mb=8000,
        tasks=1,
        time='5h',
        nodes=1,
        account='pi-blekhman',
        partition='blekhman'
    conda:
        'anvio-dev-no-update'
    shell:
        """
        anvi-script-reformat-fasta {input.fasta} \
        -o {output.reformattted} \
        --simplify-names \
        -l 750 \
        --report-file {output.report} > {log.out} 2> {log.err}
        """
rule makeIndex:
    resources:
        cpus_per_task = 4,
        mem_mb = 55000,
        tasks = 1,
        account = 'pi-blekhman',
        partition = 'blekhman'
    conda:
        'biobakery3'
    input:
        fasta='02_ASSEMBLY/{sample}/contigs.simplified.fa'
    output:
        indexLog='03_indexes/index{sample}.log',
        index = '03_indexes/{sample}.1.bt2'
    params:
        index='03_indexes/{sample}'
    threads: 4
    shell:
        """
        bowtie2-build {input.fasta} 03_indexes/{params.index} --threads {threads} > {output.indexLog}
        """

rule align:
    input:
        R1='/project/blekhman/jjcolgan/mgDenovoAssembly/populationAssembly/01_QC/{sample}_R1_dehosted.fastq.gz',
        R2='/project/blekhman/jjcolgan/mgDenovoAssembly/populationAssembly/01_QC/{sample}_R2_dehosted.fastq.gz',
        index = '03_indexes/{sample}.1.bt2'
    log:
        bowtieAlignment='04_MG_ALIGNED/logs/{sample}Alignment.err',
        bowtieAlignmentOut='04_MG_ALIGNED/logs/{sample}Alignment.out'
    output:
        sam=temp('04_MG_ALIGNED/{sample},[\d\w]+_[\d\wy]+}MgAligned.sam'),
        bam='04_MG_ALIGNED/{sample},[\d\w]+_[\d\wy]+}MgAligned.bam'
    params:
        index = '03_indexes/{sample}'
    threads: 4
    resources:
        cpus_per_task=4,
        mem_mb=12000,
        tasks=1,
        time='5h',
        nodes=1,
        account='pi-blekhman',
        partition='blekhman'
    conda:
        'biobakery3'
    shell:
        """
        bowtie2 -x {params.index} -1 {input.R1} -2 {input.R2} \
        -p {threads} -S {output.sam}  --very-sensitive > {log.bowtieAlignment} \
        2> {log.bowtieAlignmentOut}
        samtools view -bS {output.sam} > {output.bam}
        """
rule filter:
    threads: 1
    resources:
        cpus_per_task=1,
        mem_mb=2000,
        tasks=1,
        time='5h',
        nodes=1,
        account='pi-blekhman',
        partition='blekhman'
    conda:
        'biobakery3'
    input:
        bam=temporary('04_MG_ALIGNED/{sample},[\d\w]+_[\d\wy]+}MgAligned.bam')
    output:
        filtered='04_MG_ALIGNED/{sample,[\d\w]+_[\d\wy]+}mapped.bam'
    log:
        filterLog='04_MG_ALIGNED/logs/{sample,[\d\w]+_[\d\wy]+}Filter.log',
        filterError='04_MG_ALIGNED/logs/{sample,[\d\w]+_[\d\wy]+}Filter.err'
    threads: 1
    shell:
        """
        samtools view -b -F 4 -q 1 {input.bam} -o {output.filtered} > {log.filterLog} 2> {log.filterError}
        """
rule sort:
    threads: 1
    resources:
        cpus_per_task=1,
        mem_mb=2000,
        tasks=1,
        time='2h',
        nodes=1,
        account= 'pi-blekhman',
        partition='blekhman'
    conda:
        'biobakery3'
    input:
        filtered=temp('04_MG_ALIGNED/{sample}mapped.bam')
    output:
        sorted=temp('04_MG_ALIGNED/{sample,[\d\w]+_[\d\wy]+}Sorted.bam')
    log:
        out = '04_MG_ALIGNED/{sample,[\d\w]+_[\d\wy]+}Sorted.log',
        err ='04_MG_ALIGNED/{sample,[\d\w]+_[\d\wy]+}Sorted.err'
    threads: 1
    shell:
        """
        samtools sort {input.filtered} -o {output.sorted} > {log.out} 2> {log.err}
        """

rule index_bam_konzo:
    resources:
        cpus_per_task = 5,
        mem_mb = 15000,
        tasks = 1,
        time = '15h',
        nodes = 1,
        account = 'pi-blekhman'
    conda:
        'biobakery3'
    input:
        sorted='04_MG_ALIGNED/{sample,[\d\w]+_[\d\wy]+}Sorted.bam'
    output:
        indexedDone='04_MG_ALIGNED/{sample}_index.done'
    log:
        out="04_MG_ALIGNED/{sample}_index.out",
        err='04_MG_ALIGNED/{sample}_index.err'
    threads:
        5
    shell:
        """
        samtools index -@ {threads} {input.sorted} > {log.out} 2> {log.err}
        touch {output.indexedDone}
        """
rule genContigs:
    resources:
        cpus_per_task=4,
        mem_mb=50000,
        tasks=1,
        time='10h',
        nodes=1,
        account='pi-blekhman',
        partition='blekhman'
    input:
        fasta = '02_ASSEMBLY/{sample}/contigs.simplified.fa'
    conda:
        'anvio-dev-no-update'
    output:
        contigs_db = '05_CONTIGS_DB/{sample}/contigs.db'
    log:
        contigsDbLog ='05_CONTIGS_DB/{sample}/konzoDb.log',
        contigsDbError='05_CONTIGS_DB/{sample}/konzoDb.err'
    threads: 8
    shell:
        """
        anvi-gen-contigs-database -f {input.fasta} -o {output.contigs_db} -T {threads} > {log.contigsDbLog} 2> {log.contigsDbError}
        """
rule annotate:
    resources:
        mem_mb=25000,
        tasks=4,
        time='10h',
        nodes=1,
        account='pi-blekhman',
        partition='blekhman'
    conda:
        'anvio-dev-no-update'
    threads:
        4
    log:
        sampleHmms='05_CONTIGS_DB/{sample}/hmm.err',
        hmmsLog='05_CONTIGS_DB/{sample}/hmm.log',
        keggOut = '05_CONTIGS_DB/{sample}/kegg.log',
        keggErr = '05_CONTIGS_DB/{sample}/kegg.err',
        scgOut ='05_CONTIGS_DB/{sample}/scg.log',
        scgErr='05_CONTIGS_DB/{sample}/scg.err',
        pfamsOut = '05_CONTIGS_DB/{sample}/pfams.log',
        pfamsErr = '05_CONTIGS_DB/{sample}/pfams.err',
        cazymesout ='05_CONTIGS_DB/{sample}/caz.log',
        cazymesErr = '05_CONTIGS_DB/{sample}/caz.log'
    input:
        contigs_db = '05_CONTIGS_DB/{sample}/contigs.db'
    output:
        '05_CONTIGS_DB/{sample}/annotations.done'
    shell:
        '''
        anvi-run-hmms -c {input.contigs_db} \
        -T {threads} \
        --just-do-it \
        > {log.hmmsLog} 2> {sampleHmms}
        
        anvi-run-kegg-kofams -c {input.contigs_db} \
        -T {threads} \
        --just-do-it \
        --kegg-data-dir kegg \
        > {log.keggOut}
        2> {log.keggOut}
        
        anvi-run-scg-taxonomy -c {input.contigs_db} \
        -T {threads} \
        --scgs-taxonomy-data-dir /project/blekhman/shared/anvioDbs/scg \
        > {log.pfamsErr} 2> {log.pfamsErr}
        
        anvi-run-pfams -c {input.contigs_db} \
        -T {threads} \
        --pfam-data-dir /project/blekhman/shared/anvioDbs/pfam \
        > {log.pfamsErr} 2> {log.pfamsErr}
        
        anvi-run-cazymes -c {input.contigs_db} \
        -T {threads} \
        --cazyme-data-dir /project/blekhman/shared/anvioDbs/cazyme \
        > {log.cazymesout} 2> {log.cazymesErr}
        touch {output}
        '''
rule profile:
    conda: 'anvio-dev-no-update'
    resources:
        cpus_per_task=10,
        mem_mb=50000,
        tasks=1,
        time='15h',
        nodes=1,
        account='pi-blekhman'
    input:
        indexedDone = '04_MG_ALIGNED/{sample}_index.done',
        contigs_db='05_CONTIGS_DB/{sample}/contigs.db',
        bam='04_MG_ALIGNED/{sample}Sorted.bam'

    log:
        out='06_MG_PROFILES/{sample}Profile.log',
        error='06_MG_PROFILES/{sample}Profile.err'
    output:
        done='06_MG_PROFILES/{sample}/PROFILE.db'
    params:
        dir = '06_MG_PROFILES/{sample}'
    threads: 10
    shell:
        """
        anvi-profile -i {input.bam} \
        -c {input.contigs_db} \
        -T {threads} \
        -o {params.dir} --force-overwrite > {log.out} 2> {log.error}
        """
