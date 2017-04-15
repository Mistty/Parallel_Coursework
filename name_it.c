#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
void cmd_input(int *n, int *t, int *c, int *MAX_ITRS, int *tile);
void move_red(int **sub_grid, int tile, int n, int k);
void move_blue(int** sub_grid, int tile, int n, int k);
void count_color(int** sub_grid, int** redcount, int** bluecount, int n, int t, int tile, int k);
int** get_memory(int r, int c);
void delete_memory(int **grid);
void free_memory(int myid, int **grid, int **sub_grid, int **redcount, int **bluecount, int *k_size);


void cmd_input(int *n, int *t, int *c, int *MAX_ITRS, int *tile){
	printf("Please input parameters below: \n");
	printf("Grid Size: ");
	fflush(stdout);
	scanf("%d", n);
	printf("Tile Size: ");
	fflush(stdout);
	scanf("%d", t);
	printf("Percentage: ");
	fflush(stdout);
	scanf("%d", c);
	printf("Max iterations: ");
	fflush(stdout);
	scanf("%d", MAX_ITRS);
	*tile = *n/ *t;
	printf("Input Finished, Computation Start...\n");
}
/////////////////////////////////////////////////////////////////////////////////////////////
// Red and Blue movement procedure
/*
white = 0, red = 1, blue = 2,
red or blue just moved in = 3 and
red or blue (in the first row or column) just moved out = 4
*/

//grid initialization
/*
void grid_init(int n, int **sub_grid){
    int i, j, r, t[3] = {0,0,0}, m[3];
    switch(n * n % 3){
        case 0:
            m[0] = n * n / 3;
            m[1] = n * n / 3;
            m[2] = n * n / 3;    
            break;
        case 1:
            m[0] = n * n / 3;
            m[1] = n * n / 3;
            m[2] = n * n / 3 + 1;    
            break;
        case 2:
            m[0] = n * n / 3 + 1;
            m[1] = n * n / 3 + 1;
            m[2] = n * n / 3;    
            break;
    }
    
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++){
            r = rand() % 3;
            while(t[r] + 1 > m[r])
                r = rand() % 3;
            t[r]++;
            sub_grid[i][j] = r;
        }
}
*/

void move_red(int **sub_grid, int tile, int n, int k){
	int i, j;
	for(i=1;i <tile*k+1;i++){
       if (sub_grid[i][0]== 1 && sub_grid[i][1]== 0){
            sub_grid[i][0]= 0;
            sub_grid[i][1]= 3;
        }
        for (j=1;j<n;j++){
            if (sub_grid[i][j]== 1 && sub_grid[i][(j+1)%n]== 0){
                 sub_grid[i][j]= 0;
               sub_grid[i][(j+1)%n]= 3;
            }
            else if (sub_grid[i][j]== 3)
                sub_grid[i][j]= 1;
        }
        if (sub_grid[i][0]== 3)
            sub_grid[i][0]= 1;
    }
}
void move_blue(int** sub_grid, int tile, int n, int k){
	int i, j;
	for(j= 0;j<n;j++){
		if (sub_grid[0][j]== 2 && sub_grid[1][j]== 0){
			sub_grid[0][j]= 0;
			sub_grid[1][j]= 3;
		}
		for(i =1; i<tile*k+2; i++){
			if (sub_grid[i][j] == 2 && sub_grid[(i+1)%(tile * k + 2)][j] == 0){
                sub_grid[i][j] = 0;
                sub_grid[(i+1)%(tile * k + 2)][j] = 3;
            }
			else if (sub_grid[i][j]== 3){
                sub_grid[i][j]= 2;
			}
		}	
        if (sub_grid[0][j]== 3){
            sub_grid[0][j]= 2;
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////
void count_color(int** sub_grid, int** redcount, int** bluecount, int n, int t, int tile, int k){
	int i, j;
	for (i = 0; i < k; i++){
		for (j = 0; j < t; j++){
			redcount[i][j] = 0, bluecount[i][j] = 0;
		}
	}
	for (i = 1; i < tile*k+1; i++){
		for(j=0; j< t; j++){
			if(sub_grid[i][j] == 1){
				redcount[(i-1)/tile][j/tile]++;
			}
			if(sub_grid[i][j] == 2){
				bluecount[(i-1)/tile][j/tile]++;
			}
		}
	}
}
// get memory
int** get_memory(int r, int c){
    int i;
    int **grid = (int **) malloc(r * sizeof(int*));
    grid[0] = malloc(r * c * sizeof(int));
    for(i = 0; i < r; ++i)
        grid[i] = grid[0] + i * c;
    return grid;
}

// release memory for the grids
void delete_memory(int **grid){
    free(grid[0]);
    free(grid);
}

// release memory
void free_memory(int myid, int **grid, int **sub_grid, int **redcount, int **bluecount, int *k_size){
    if(myid == 0){
        delete_memory(grid);
        free(k_size);
    }
    delete_memory(sub_grid);
    delete_memory(redcount);
    delete_memory(bluecount);
    
}


int main(int argc, char** argv){
	int n, t, c, MAX_ITRS;
	int tile;
	int k;						//counts of tiles allocated to each processor
	int *k_size;					//total tile numbers
	int myid, numprocs;			//processor id, numbers of processors
	int **redcount, **bluecount;//the numbers of red and blue
	int **grid; 				//the main grid
	int **sub_grid; 			//sub grid allocated to each processor
	int end = 0, n_itrs = 0;	//end signal and number of iterations
	int currentTotal = 0;
	int i, j, u;
	int sr;
	int used_pro;
	
	MPI_Status status;
	
	//initialize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	
	    /*
    *  initialize grid and get user input
    */
    
    if(myid == 0){        
	
        cmd_input(&n,&t,&c,&MAX_ITRS,&tile);		// get input parameters
		
        k_size = (int*)malloc(numprocs * sizeof(int));
		
		for(u = 0; u < numprocs; u++){
			if(numprocs > t){
				if(u < t){
					k_size[u] = 1;
				}
				else{
					k_size[u] = 0;
				}
			}
			else{ 
				if(u < t % numprocs){
					k_size[u] = t / numprocs + 1;
					}
				else{
					k_size[u] = t / numprocs;
				}
			}
		}
        //getk_size(pn, k_size, t); 		// number of tiles for each process
		
        used_pro = (numprocs > t)? t:numprocs; 			// number of used procossor
			
        grid = get_memory(n, n); 						// allocate memory for grid
		
		//grid_init(n, grid);
		
		srand((unsigned)time(NULL));					//use time as a seed to generate random numbers
		for (i = 0; i < n; ++i) {
			for (j = 0; j < n; ++j) {
				grid[i][j] = rand() % 3;
			}
		}
		
		for (i =0; i<n; i++){
			for(j=0;j<n;j++){
				printf(" %d -", grid[i][j]);
			}
			printf("\n");
		}
	
        //board_init(n, board); // init grid
        printf("The grid above is the initialized grid for this computation.\n");
        fflush( stdout );
    }

    // broadcast the input parameters to other processors
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&c, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&MAX_ITRS, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
    MPI_Bcast(&tile, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
    MPI_Bcast(&used_pro, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    //send the sub grid size
    if(myid == 0){
		
        k = k_size[0]; //get k
        sr = 0; // start tile row    
        sub_grid = get_memory(tile * k + 2, n);  			// allocate memory for sub grid
		

        // get the sub grid for first processor
        for(i = 1; i < tile * k + 1; i++)
            for(j =0; j < n; j++)
                sub_grid[i][j] = grid[i - 1][j];

        // send parameters to other processors
        currentTotal += k_size[0];
        for(i = 1; i < numprocs; i++){
            MPI_Send(&currentTotal, 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD);
            MPI_Send(&k_size[i], 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD);
            MPI_Send(grid[currentTotal * tile], tile * k_size[i] * n , MPI_INTEGER, i, 1, MPI_COMM_WORLD);
            currentTotal += k_size[i];
        }     
          
    }else{
        // get parameters from processor 0
        MPI_Recv(&sr, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&k, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, &status);
        sub_grid = get_memory(tile * k + 2, n);
        MPI_Recv(sub_grid[1], tile * k * n, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, &status);
    }
    
	// get memory for red
	redcount = get_memory(k, t);
	
	//get memory for blue
    bluecount = get_memory(k, t);
  
    /*
    * red/blue movement
    * check termination standard
    */
    while (!end && n_itrs < MAX_ITRS){
        n_itrs++;

        if(k == 0){
            // ignore idle process
            MPI_Allreduce(&end, &end, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD); 
            continue;
        }
            
        move_red(sub_grid, tile, n, k); // red movement

        // exchange rows
        MPI_Sendrecv(sub_grid[1], n, MPI_INT, (used_pro + myid - 1) % used_pro, 1, sub_grid[tile * k + 1], n, MPI_INT, (myid + 1) % used_pro, 1, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(sub_grid[tile * k], n, MPI_INT, (myid + 1) % used_pro, 1, sub_grid[0], n, MPI_INT, (used_pro + myid - 1) % used_pro, 1, MPI_COMM_WORLD, &status);

        move_blue(sub_grid, tile, n, k); // blue movement

        // detect termination status
        count_color(sub_grid, redcount, bluecount, n, t, tile, k);
        for(i = 0; i < k; i++)
            for(j = 0; j < t; j++){
                if((float)redcount[i][j] / (tile * tile) > (float)c / 100 || (float)bluecount[i][j] / (tile * tile) > (float)c / 100){
                    printf("(%d, %d) has the colored squares more than %d%%\n", sr + i, j, c);
                    fflush( stdout );
                    end = 1;
                }
            }

        // reduce termination status to all
        MPI_Allreduce(&end, &end, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD); 
    }

    /*
    * print result
    */
    // print if termination with maximum iteration number reached
    if(myid == 0)
        if(!end){
            printf("Maximum number of iteration reached\n");
            fflush( stdout );
        }
    
    // gather all sub grids and print current grid
    if(myid == 0){
        currentTotal = k;
        for(i = 1; i < tile * k + 1; i++){
            for(j =0; j < n; j++){
                grid[i - 1][j] = sub_grid[i][j];
			}
		}
        for(i = 1; i < numprocs; i++){
            MPI_Recv(grid[currentTotal * tile], tile * k_size[i] * n , MPI_INTEGER, i, 1, MPI_COMM_WORLD, &status);
            currentTotal += k_size[i];
        }
        
        fflush( stdout );
		
        for(i = 0; i < n; i++){
			for(j = 0; j < n; j++){
				printf("%d ", grid[i][j]);
				if((j+1) % t == 0)
					printf(" ");  
			}
			printf("\n");
			if((i+1) % t == 0)
				printf("\n");
    }
		
    }else{
		
        MPI_Send(sub_grid[1], tile * k * n, MPI_INTEGER, 0, 1, MPI_COMM_WORLD);
    }
    free_memory(myid, grid, sub_grid, redcount, bluecount, k_size); 

    MPI_Finalize();

    return 0;
	
}
