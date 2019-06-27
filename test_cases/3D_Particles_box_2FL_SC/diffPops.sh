
it=$2

echo "Pops:$1 ---------------------"
./maxDiff_1_32.sh "$1""$it".000000.dat  master/"$1""$it".000000.dat 5 6

echo "Vars:$1 ---------------------"
./maxDiff_1_32.sh "$1""$it".000000.dat1 master/"$1""$it".000000.dat1 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18

echo ""
