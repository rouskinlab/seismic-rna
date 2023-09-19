
.. _sam-format:

SAM: Sequence Alignment Maps
------------------------------------------------------------------------

Aligned reads are stored as sequence alignment map (SAM) files. These
files come in three formats -- SAM, BAM, and CRAM -- documented on the
`Samtools website <https://samtools.github.io/hts-specs/>`_:

======================== ========= ====== =========== ========= ========= ===================
Format                   Extension Type   Readability I/O Speed File Size Uses in SEISMIC-RNA
======================== ========= ====== =========== ========= ========= ===================
Sequence Alignment Map   ``.sam``  text   \+\+\+      \+\+      \-\-\-    parsing and editing
Binary Alignment Map     ``.bam``  binary \+\+        \+\+\+    \-\-      short-term storage
CompRessed Alignment Map ``.cram`` binary \+          \+        \-        long-term storage
======================== ========= ====== =========== ========= ========= ===================
