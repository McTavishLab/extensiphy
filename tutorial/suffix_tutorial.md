---

title: Suffix Tutorial

author: Jasper Toscani Field

output: html_document

---

## Short read paired end file suffixes

The files for paired end reads come in paris - usally signified by an R1 and R2 before the file extention.

e.g.

```
taxon_name_R1.fastq
taxon_name_R2.fastq
```
in order to recognize and pair up those files, Extensiphy need to know that the pairs are denoted by:

```
_R1.fastq
_R2.fastq
```

Depending on how the data was generated or processed, the suffixes can be thinge like 'r1_001.fq' and 'r2_001.fq'.
So for paired end reads, pass in the suffix and the extention to -1 and -2.


