# Identify rescued contained reads by hifiasm (?)

grep -P "^A\t" output.bp.r_utg.gfa | grep -o -P "read=[^\t]*" | sort > output.bp.r_utg.gfa.sorted.headers
comm -12 output.bp.r_utg.gfa.sorted.headers ../map_mm_noncontained/contained.100.txt > output.bp.r_utg.gfa.contained.headers