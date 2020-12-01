datasets=["ERR1050075_1.fastq", "ERR992657_2.fastq", "SRR12801265.fastq", "SRR12873993.fastq", "SRR12989394.fastq", "SRR8653791.fastq"]

size=["100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", "1100", "1200", "1300", "1400", "1500", "1600", "1700", "1800", "1900", "2000"]

rule all:
  input: expand("size.{dataset}.{size}M",dataset=datasets,size=size)

rule size:
    input:
        "{dataset}"
    output:
        "size.{dataset}.{size}M"
    log:
        "size.{dataset}.{size}M.log"
    shell:
        """
		./sizeTest -d kmers  -k {input}  -n {wildcards.size}000000 > {output} 2> {log}
	"""

 
