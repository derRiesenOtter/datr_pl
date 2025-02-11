library(rtracklayer)
library(Biostrings)
library(stringr)

# report results gffread

print("Report results of gffread")

transcriptome <- readDNAStringSet("../material/1_gffread_transcriptome.fa")

number_of_transcripts <- length(transcriptome)

print("Number of reference Transcripts generated:")
print(number_of_transcripts)

gff <- import("../material/reference.gff", format = "gff")

transcribed_types <- c(
  "mRNA",
  "tRNA",
  "snoRNA",
  "rRNA",
  "snRNA",
  "ncRNA",
  "telomerase_RNA",
  "transposable_element_gene"
)

number_of_mrna_seq <- sum(gff$type %in% transcribed_types)

print("Number of mRNA sequences:")
print(number_of_mrna_seq)

print("All Transcripts have been translated. Transcribable elements are:")
print(transcribed_types)

number_of_genes <- sum(gff$type == "gene")

sprintf("The reference.gff file contained %d genes.", number_of_genes)
print("Splice variants are responsible for the difference in transcripts and genes.")

genes_without_cds <- sub(" .*", "", names(transcriptome))
genes <- sub("_.*", "", genes_without_cds)
number_of_splice_vairants <- table(genes)

genes_with_one_splice <- sum(number_of_splice_vairants == 1)
genes_with_two_splice <- sum(number_of_splice_vairants == 2)

number_transcribed_elements <- sum(gff$type %in% transcribed_types)

sprintf(
  "Of the %d transcribed elements in the gff, %d have only one splice variant, while %d have two.",
  number_transcribed_elements,
  genes_with_one_splice,
  genes_with_two_splice
)
