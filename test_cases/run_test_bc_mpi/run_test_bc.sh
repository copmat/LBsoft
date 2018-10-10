#!/bin/bash

  #! 0 F F F
  #! 1 T F F
  #! 2 F T F
  #! 3 T T F
  #! 4 F F T
  #! 5 T F T
  #! 6 F T T
  #! 7 T T T

declare -a pbcx=('0' '1' '0' '1' '0' '1' '0' '1')

declare -a pbcy=('0' '0' '1' '1' '0' '0' '1' '1')

declare -a pbcz=('0' '0' '0' '0' '1' '1' '1' '1')

declare -a dtype=('1' '2' '3' '4' '5' '6' '7')

declare -a ddime=('1 1 4' '1 4 1' '4 1 1' '1 2 2' '2 1 2' '2 2 1' '2 2 2')

declare -a ncpu=('4' '4' '4' '4' '4' '4' '8')

#declare -a ncpu=('1' '1' '1' '1' '1' '1' '1')

itest=()

i=0
ifine=7

WDR=`pwd`

cd $WDR

make

myprog="../../../execute/main.x"
myprogp="../../../execute/main_mpi.x"
compare="../compare.x"

k=1

while (( $i <= $ifine )); do
  cd $WDR

  
  
  j=0
  jfine=6
  while (( $j <= $jfine )); do
    cd $WDR
    mio3=$(printf "%03d\n" $k)
    echo  "case: $k boundary condition ${pbcx[$i]} ${pbcy[$i]} ${pbcz[$i]}"
    echo  "case: $k decomposition type ${dtype[$j]} decomposition dimension ${ddime[$j]}"
    echo "mkdir -p $mio3"
    mkdir -p $mio3
    cd $WDR/$mio3
    cp -f ../input.dat .
    sed -i "/boundary condition/ c\boundary condition ${pbcx[$i]} ${pbcy[$i]} ${pbcz[$i]}" input.dat
    sed -i "/decomposition type/ c\decomposition type ${dtype[$j]}" input.dat
    sed -i "/decomposition dimension/ c\decomposition dimension ${ddime[$j]}" input.dat
    echo "cp -f $myprog ."
    echo "cp -f $myprogp ."
    cp -f $myprog .
    cp -f $myprogp .
    echo "execute in serial case: $k"
    echo "./main.x | tee test.serial.txt"
    ./main.x | tee test.serial.txt
    mv output000000.map output000000.map.orig
    mv global.map global.map.orig
    echo "execute in parallel case: $k"
    echo "mpirun -np ${ncpu[$j]} ./main_mpi.x | tee test.mpi.txt"
    mpirun -np ${ncpu[$j]} ./main_mpi.x | tee test.mpi.txt
    rm main.x main_mpi.x time.dat
    cp -f $compare .
    echo "comparison case: $k"
    ./compare.x
    echo "end case    : $k"
    mioout=`echo ${textArray[$i]} | awk '{print $1}' test.out`
    itest[$k]=$mioout
    let j++
    let k++
  done

  let i++
  
done

let k--
kfine=k

cd $WDR

if [ -f "number_case.txt" ]; then
  rm -f number_case.txt
fi
echo "$k" >> number_case.txt


if [ -f "myresults.out" ]; then
  rm -f myresults.out
fi

k=1
i=0
while (( $i <= $ifine )); do
  cd $WDR
  j=0
  jfine=6
  while (( $j <= $jfine )); do
    cd $WDR
    echo "case: $k boundary condition ${pbcx[$i]} ${pbcy[$i]} ${pbcz[$i]} on ${ncpu[$j]} CPUs"
    echo "case: $k decomposition type ${dtype[$j]} decomposition dimension ${ddime[$j]}"
    echo "case: $k boundary condition ${pbcx[$i]} ${pbcy[$i]} ${pbcz[$i]} on ${ncpu[$j]} CPUs" >> myresults.out
    echo "case: $k decomposition type ${dtype[$j]} decomposition dimension ${ddime[$j]}" >> myresults.out
    abm=${itest[$k]}
    if [ $abm == "1" ]; then
      echo "case: $k RESULT OK"
      echo "case: $k RESULT OK" >> myresults.out
    else 
      echo "case: $k RESULT NO"
      echo "case: $k RESULT NO" >> myresults.out
    fi
    let j++
    let k++
  done

  let i++
  
done


exit 0
