# Relatedness-Filtering
Script for greedily remove related samples

# Installation
```
g++ -std=c++11 main.cpp misc.cpp -o greedy_related
```

# Instruction
You can just type
```
greedy_related -h
```
To get all the instructions

To use with plink relatedness file, you can use the following awk script

```{bash}
awk -v NUM=1 'NR==1{print $2" "$4" "$10} NR!=1{print $2" "NUM" "$10; print $4" "NUM" "$10; NUM=NUM+1}'
```

Please note that **greedy_related** assume the input has a header line
