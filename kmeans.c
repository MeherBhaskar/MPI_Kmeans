#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

int n;
int k;
int d;


float getDistance(float *v1, float *v2, int d)
{
	float dist = 0.0;
	for(int i = 0; i < d ; i++)
	{
		dist += (v1[i] - v2[i])*(v1[i] - v2[i]);
	}
	return dist;
}

int main(int argc, char *argv[])
{
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	if(rank == 0)
	{
		printf("Enter number of instances per process : \n");
		scanf("%d", &n);
		printf("Number of instances :::::: %d \n",n*size);
		printf("Enter the number of means : \n");
		scanf("%d", &k);
		printf("Enter the number of attributes : \n");
		scanf("%d", &d);
	}
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&d, 1, MPI_INT, 0, MPI_COMM_WORLD);

	float entries[n*size][d];
	float means[k][d];
	float procentries[n][d];
	int   procClosestMean[n];
	int   closestMean[n*size];
	int   procMeanCount[k];		//Reduce at the root for getting new centers
	int   meanCount[k];
	
	for(int i = 0; i < k ;i++)
	{
		procMeanCount[i] = 0;
	}
	if(rank == 0)
	{
		printf("Enter the instances : \n");
		for(int i = 0; i < n*size; i++)
		{
			for(int j = 0; j < d; j++)
			{
				scanf("%f", &entries[i][j]);
			}
		}

		printf("Enter the initial means : \n");
		for(int i = 0; i < k; i++)
		{
			for(int j = 0; j < d; j++)
			{
				scanf("%f", &means[i][j]);
			}
		}

		printf("\n\n\n\n");
		printf("Initial Means : \n");
		for(int i = 0; i < k ; i++)
		{
			for(int j = 0; j < d ; j++)
			{
				printf(" %f  ",means[i][j]);
			}
			printf("\n");
		}

	}

	double t1, t2;

	t1 = MPI_Wtime();
	MPI_Scatter(entries, n*d , MPI_FLOAT, procentries, n*d, MPI_FLOAT,0, MPI_COMM_WORLD);	
	int flag2 = 1;
	int flag = 1;
	int itr = 0;

	while(flag2 >= 1)
	{	
		flag2 = 0;	
		flag = 0;
		itr++;
		MPI_Barrier(MPI_COMM_WORLD);
		flag = 0;
		MPI_Bcast(means, k*d, MPI_FLOAT, 0, MPI_COMM_WORLD);
			
		for(int i = 0; i < k ; i++)
		{
			meanCount[i] = 0;
			procMeanCount[i] = 0;
		}
			
		for (int i = 0; i < n; i++)
		{
			float currDist = 99999999;
			int currClosest = -1;
			
			for (int j = 0; j < k; j++)
			{
				if(getDistance(procentries[i], means[j], d) < currDist)
				{
					currDist = getDistance(procentries[i], means[j], d);
					currClosest = j;
				}
			}	
			if(currClosest != procClosestMean[i])
			{
				flag = 1;
			}

					procClosestMean[i] = currClosest;
				procMeanCount[currClosest]++;
		}
		
		MPI_Allreduce(&flag, &flag2, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(procClosestMean, n, MPI_INT, closestMean, n, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Reduce(procMeanCount, meanCount, k, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

		
		if(rank == 0)
		{
			float tempsums[k][d];

			for(int i =0; i < k; i++)
			{
				for(int j = 0; j < d ; j++)
				{
					tempsums[i][j] = 0;
				}
			}
			for(int i = 0; i < n*size;  i++)
			{
				for(int j = 0; j < d; j++)
				{
					tempsums[closestMean[i]][j] += entries[i][j];
				}
			}

			for(int i = 0; i < k; i++)
			{
				for(int j = 0; j < d ; j++)
				{
					if(meanCount[i] != 0)
						means[i][j] = tempsums[i][j]/meanCount[i];
				}
			}

		}
		
		MPI_Barrier(MPI_COMM_WORLD);
	}


	t2 = MPI_Wtime();
	float time = t2 - t1;
	float max_time;
	MPI_Reduce(&time, &max_time, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD );
	


	if(rank == 0)
	{
		printf("\nUpdated Means : \n");
		
			for(int i = 0; i < k; i++ )
			{
				for(int j = 0; j < d; j++)
				{
					printf(" %f  ",means[i][j]);
				}
				printf("\n");
			}
		
		printf("Time Taken : %f \n",max_time );
	}


	MPI_Finalize();
	return 0;
}