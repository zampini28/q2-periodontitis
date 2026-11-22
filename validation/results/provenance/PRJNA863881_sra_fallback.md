# Deviation: PRJNA863881 run SRR20721137, reverse reads

On 2026-11-21 the ENA file ftp.sra.ebi.ac.uk/vol1/fastq/SRR207/037/SRR20721137/SRR20721137_2.fastq.gz
was an empty directory (HTTP 301 to an empty index; FTP listing: drwxr-xr-x),
although the ENA file report lists it (121779007 bytes, md5 3b25ab79cbfaad68123febe0d68666a6).
The run was fetched from NCBI SRA with fasterq-dump 3.1.1 (--split-files).
SRA forward reads were identical in sequence and quality to the ENA forward file
(610040 reads). The SRA reverse reads (610040) were used, with headers rewritten to
the ENA form "@SRR20721137.<n> <n>/2"; sequences and qualities unchanged.
PRJNA863881.runs.tsv: md5 of that file replaced by 1f547799390f2aadf0cb5214c6cd7637;
the original ENA report is PRJNA863881.runs.ena-original.tsv.
