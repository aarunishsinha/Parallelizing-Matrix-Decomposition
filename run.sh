n=$1
filename=$2
threads=$3
strt=$4
if [[ ${strt} == "0" ]];
    then ./omp_strats $n $filename $threads $strt
fi
if [[ ${strt} == "1" ]];
    then ./omp_strats $n $filename $threads $strt
fi
if [[ ${strt} == "2" ]];
    then ./omp_strats $n $filename $threads $strt
fi
if [[ ${strt} == "3" ]];
    then ./omp_strats $n $filename $threads $strt
fi
