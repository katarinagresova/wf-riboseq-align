# star_transcript's alignment settings. contaminant_compete_align uses the same
# ones: it compares a read's contaminant AS with the AS star_transcript would
# give it.
RIBO_TRANSCRIPTOME_STAR_ARGS = (
    "--seedSearchLmax 10 "
    "--outFilterMultimapNmax 255 "
    "--outFilterMismatchNmax 2 "
    "--outFilterMultimapScoreRange 0 "
    "--outFilterIntronMotifs RemoveNoncanonical"
)
