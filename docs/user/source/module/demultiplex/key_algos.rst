
Key algorithm(s)
++++++++++++++++++++++++

The main algorithm of this code involves using a grep function from seqkit (https://bioinf.shenwei.me/seqkit/) which is essentially just grep for fastq. 
For those that do not know how grep works, please read this article (https://en.wikipedia.org/wiki/Grep). The demultiplexing code here relies on this tool to query
the fastq for specific patterns that connect it to a reference sequence. This usually means querying for a sequnece of nucleotides unique to each reference.

For example:

    