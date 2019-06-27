awk 'BEGIN{OFMT="%06.3f"} {if (NR%2==1)print $7+0.0,$8+0.0,$9+0.0,$2}' | ./trisort.sh 1 2 3
