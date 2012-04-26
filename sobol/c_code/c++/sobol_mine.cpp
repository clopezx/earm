# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

// sample compilation: g++ -c -g sobol.cpp sobol_mine.cpp                     

using namespace std;

# include "sobol.hpp"

int main ( );
void test04 ( );
//void test08 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    SOBOL_PRB calls a set of problems for SOBOL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  //  test04 ( );
  test08 ( );
  return 0;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//  Licensing:
//    This code is distributed under the GNU LGPL license. 
//  Author: John Burkardt
{
# define DIMS 137 // Dimensions
# define SEED_SKIP 500 // beginning seed point
# define SEED_MAX  1500 // end seed point. SEED_MAX - SEED_SKIP is the number of points

  int dim_num;
  int i;
  int j;
  float r[DIMS];
  int seed;
  int seed_in;
  int seed_out;

  cout << "# Using I4_SOBOL algorithm to sample " << DIMS << " dimensions with " << (SEED_MAX - SEED_SKIP) << " samples \n";
  seed = SEED_SKIP;
  dim_num = DIMS;

  for ( i = SEED_SKIP; i < SEED_MAX; i++ )
    {
      seed_in = seed;
      i4_sobol ( dim_num, &seed, r );
      seed_out = seed;
      
      for ( j = 0; j < dim_num; j++ )
	{
	  cout << setw(6) << r[j];
	  if (j+1 == dim_num){
	    cout << "\n";
	  }
	  else
	    {
	      cout << ", ";
		}
	}
    }
  return;
# undef DIMS
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//    This code is distributed under the GNU LGPL license. 
//  Author: John Burkardt
{
# define DIMS 2
# define SEED_SKIP 64 
# define SEED_MAX 1064

  int dim_num;
  int i;
  int j;
  double r[DIMS];
  long long int seed;
  long long int seed_in;
  long long int seed_out;

  cout << "# Using I8_SOBOL algorithm to sample " << DIMS << " dimensions with " << (SEED_MAX - SEED_SKIP) << " samples \n";

  for ( i = SEED_SKIP; i < SEED_MAX; i++ )
    {
      seed_in = seed;
      i8_sobol ( dim_num, &seed, r );
      seed_out = seed;
      
      for ( j = 0; j < dim_num; j++ )
        {
          cout << setw(14) << r[j];
	  if (j+1 == dim_num){
	    cout << "\n"
	      }
	  else
	    {
	      cout << ",";
	    }

	}
    }
  
}

  return;
# undef DIM_MAX
}
//****************************************************************************80
