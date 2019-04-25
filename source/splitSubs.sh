#!/bin/sh

mkdir split
cd split

curDir=$(pwd)

nowDir=$(basename $curDir)

if [ $nowDir != "split" ]; then
	echo "Couldn't change to split dir. Exiting"
	exit 1
fi

for i in $*
do
	echo "Splitting $i"
	awk -v name=$i '
		/^!/ {
			next
		}

		/end *subroutine/ {
			#print "fine", $3
			inSub = 0
			next
		}

		/ *subroutine/ {
			split($2,nn,"(")
			#print "inizio", nn[1]
			inSub = 1
			fname = name "." nn[1]
			next
		}

		{
			if (inSub) print $0 > fname
	}' ../$i
done
