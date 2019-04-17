# MPI_Kmeans
An MPI implementation of K Means Clustering Algorithm to run KMeans in parallel.

## MPI Installation 
```sudo apt install mpich```

Once MPI is installed, compile the program using 
```mpicc kmeans.c -o obj```

Run the program using num processes with the command : 
```mpiexec -n $num ./obj```
