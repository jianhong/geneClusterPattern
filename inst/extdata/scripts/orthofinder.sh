for f in *fa ; do python orthofinder/bin/primary_transcript.py $f ; done
orthofinder -a 12 -f primary_transcripts/
