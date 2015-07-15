SET(SOURCES
  sequence/Nucleotide.C sequence/AminoAcid.C sequence/NTSequence.C
  sequence/AASequence.C sequence/Codon.C sequence/Mutation.C
  sequence/CodingSequence.C evolution/NucleotideSubstitutionModel.C
)  

ADD_LIBRARY(seq ${SOURCES})

INCLUDE_DIRECTORIES(${SEQ_SOURCE_DIR}/src/sequence ${SEQ_SOURCE_DIR}/src/evolution)
INSTALL_TARGETS(/lib seq)
SUBDIRS(sequence evolution)
