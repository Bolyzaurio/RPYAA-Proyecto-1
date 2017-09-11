#include "nasa_shuttle.h"
/**
   Naive Bayes Formula
   P( Y = y | X = x1, x2, ... , xm) = 
   PRODUCTO( P (xi|y) P(Y=y) /
   SUMA ( PRODUCTO ( P (xi | y) , P( Y = yk)))
**/
double variance [CLASSES] [COLUMNS - 1] = {0};
double mean [CLASSES] [COLUMNS - 1]= {0};
double totals [CLASSES] [COLUMNS - 1]= {0};
int main (int argc, char *arcs[] )
{
  gsl_matrix *data = read_data();
  read_new_set();
  gsl_matrix_free(data);
  return 0;
}
void read_new_set()
{
  FILE *data;
  size_t bytes;
  int read;
  char *token = " ";
  char * line = NULL;
  char *atom = NULL;
  data = fopen("shuttle.tst","rb");
  double auxiliary_vector[COLUMNS-1] ={0};
  int j, class, aux, candidata , i , true_class,incorrecto,correcto;
  int coso = 25;
  long double probability = 0;
  long double biggest = 0;
  long double zero = 0.0;
  while( (read = getline(&line, &bytes , data)) != -1 && coso)
    {
      biggest = 0;
      candidata = 0;
      class = 0;
      probability = 0.0;
      for (atom = strtok( line , token ), j = 0; atom != NULL ; atom = strtok( NULL , token ) , j++ )
	if ( atom != NULL)
	  {
	    auxiliary_vector[j] = atof(atom);
	    true_class = atoi(atom);
	  }
      for ( ;class < CLASSES ; class ++ , probability = 0.0)
	{//Reset Candidata, Probability, 
	  for (aux = 0 ; aux < COLUMNS - 1 ; aux ++ )
	    {
	      if (probability == 0.0)
		probability = calculate_probability( mean [class][aux] , variance [class][aux],  auxiliary_vector[aux], 0);
	      probability += calculate_probability( mean [class][aux] , variance [class][aux],  auxiliary_vector[aux], 0);
	    }
	  
	  printf("Iteracion  %d Biggest %Lf Probability %Lf\n",class, biggest , probability);
	  
	  if (biggest == zero)
	    {
	      //printf("\t\tBiggest = %Lf , Prob = %Lf candidata %d\n", biggest, probability, candidata);
	      biggest = probability;
	      candidata = class;
	      //printf("\t\tBiggest = %Lf\n",biggest);
	    }
	  if (biggest < probability)
	    {
	      biggest = probability;
	      candidata = class;
	    }
	}
      printf("Ejemplo %d : Probabilidad mas alta es de pertenecer a la clase  %d con :  %Lf\n", i , candidata+1 , biggest);
      if ( candidata + 1 == true_class)
	{
	  printf("Correcto\n");
	  correcto++;
	}
      else
	{
	  printf("Incorrecto\n");
	  incorrecto++;
	}
      //coso--;
      i++;
    }
  printf("\tCorrrectos %d \t\tIncorrectos%d\n", correcto,incorrecto);
  fclose(data);
}
// Gaussian_Distribution
long double calculate_probability ( double mean , double stdev, double x , int characteristic)
{
  double exponent = exp( - ( pow ( x - mean , 2 )/
			     ( 2 * pow( stdev , 2 ) )));
  double returnable = ( 1 / (sqrt( 2 * M_PI ) * stdev)) * exponent;
  return returnable ;
}
gsl_matrix* read_data()
{
  FILE *data;
  size_t bytes;
  int read, count_lines;
  char *token = " ";
  char * line = NULL;
  char *atom = NULL;
  data = fopen("shuttle.trn","rb");

  double * count_aparitions = calloc ( sizeof(double) , CLASSES );
  
  while ( (read = getline(&line, &bytes , data)) != -1 )
    count_lines ++;
  fclose(data);
  gsl_matrix * data_m = gsl_matrix_alloc (COLUMNS,count_lines);
  data = fopen("shuttle.trn","rb");
  /**
     Read the Matrix from the file
  **/
  for ( int class , j, i = 0 ; (read = getline(&line, &bytes , data)) != -1 ; i ++ )
    {
      for (atom = strtok( line , token ), j = 0; atom != NULL ; atom = strtok( NULL , token ) , j++ )
	if ( atom != NULL)
	    {
	      gsl_matrix_set (data_m, j, i, atof(atom));
	      class = atoi(atom);
	    }
      for (int aux = 0 ; aux < COLUMNS - 1 ; aux ++ )
	totals[class-1][aux] += gsl_matrix_get (data_m , aux , i) ;
      count_aparitions[class-1]++;
    }
  for ( int i = 0 ; i < CLASSES ; i ++)
    printf("CLASSE %d aparece %lf \n", i + 1 , count_aparitions[i]);
  fclose(data);

  // Calculate mean of each caracteristic and CLASS
  // The table should be CLASS X CHARACTERISTIC

  for ( int k = 0 ; k < CLASSES ; k ++){
    for ( int j = 0 ; j < COLUMNS - 1 ; j ++)
      mean [k] [j] = ( totals [k] [j] / count_aparitions[k] ) ;
  }
  // Calculate Variance of each caracteristic and CLASS
  // The table should be CLASS X CHARACTERISTIC
  int index_class;
  for ( size_t j = 0 ; j < data_m-> size2 ; j ++ ) // Variance (CLASS,CHARACTERISTIC);
    {
      index_class = gsl_matrix_get(data_m,COLUMNS-1,j) - 1; // Value from 1 - 7
      for ( int i = 0 ; i < (COLUMNS - 1); i ++ )
	variance[index_class][i] = variance[index_class][i] + pow((gsl_matrix_get(data_m,i,j) - mean[index_class][i]), 2);
    }
  for ( int k = 0 ; k < CLASSES ; k ++)
    for ( int j = 0 ; j < COLUMNS - 1 ; j ++)
      variance [k] [j] = ( variance [k] [j] / count_aparitions[k] ) ;
  for ( int k = 0 ; k < CLASSES ; k ++)
    {
      printf("\t\tClass %d\n", k+1);
      for ( int j = 0 ; j < COLUMNS - 1 ; j ++)
	printf("Mean : %lf \t STD : %lf\n", mean[k][j] , variance[k][j] );
    }
  return data_m;
}
