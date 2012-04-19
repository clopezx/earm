# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "sobol.hpp"

int main ( );
void test04 ( );
//void test05 ( );

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
  //  timestamp ( );

  // cout << "\n";
  // cout << "SOBOL_PRB\n";
  // cout << "  C++ version\n";
  // cout << "  Test the SOBOL library.\n";
  // cout << "-----------------------------------------\n";

  test04 ( );
//  test05 ( );
//
//  Terminate.
//
  // cout << "\n";
  // cout << "SOBOL_PRB\n";
  // cout << "  Normal end of execution.\n";

  // cout << "\n";
  // timestamp ( );

  return 0;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests I4_SOBOL.
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
# define DIMS 137
# define SEED_SKIP 500 // I have seen many possibilities for this but no strong reasons to choose one... 
# define SEED_MAX  1500

  int dim_num;
  int i;
  int j;
  float r[DIMS];
  int seed;
  int seed_in;
  int seed_out;

  //cout << "\n";
  //cout << "TEST04\n";
  //cout << "  I4_SOBOL computes the next element of a Sobol sequence.\n";
  //cout << "\n";
  //cout << "  In this test, we call I4_SOBOL repeatedly.\n";
  cout << "# Using I4_SOBOL algorithm to sample " << DIMS << " dimensions with " << (SEED_MAX - SEED_SKIP) << " samples \n";

  seed = SEED_SKIP;

  dim_num = DIMS;
  // cout << "\n";
  // cout << "  Using dimension DIM_NUM =   " << dim_num << "\n";
  // cout << "\n";
  //  cout << "  Seed  Seed   I4_SOBOL\n";
  //  cout << "  In    Out\n";
  // cout << "\n";

  for ( i = SEED_SKIP; i <= SEED_MAX; i++ )
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
//
//  Purpose:
//
//    TEST08 tests I8_SOBOL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_MAX 4

  int dim_num;
  int i;
  int j;
  double r[DIM_MAX];
  long long int seed;
  long long int seed_in;
  long long int seed_out;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  I8_SOBOL computes the next element of a Sobol sequence.\n";
  cout << "\n";
  cout << "  In this test, we call I8_SOBOL repeatedly.\n";

  for ( dim_num = 2; dim_num <= DIM_MAX; dim_num++ )
  {

    seed = 0;

    cout << "\n";
    cout << "  Using dimension DIM_NUM =   " << dim_num << "\n";
    cout << "\n";
    cout << "  Seed  Seed   I8_SOBOL\n";
    cout << "  In    Out\n";
    cout << "\n";

    for ( i = 0; i <= 110; i++ )
    {
      seed_in = seed;
      i8_sobol ( dim_num, &seed, r );
      seed_out = seed;

      if ( i <= 11 || 95 <= i )
      {
        cout << setw(6) << seed_in << "  "
             << setw(6) << seed_out << "  ";
        for ( j = 0; j < dim_num; j++ )
        {
          cout << setw(14) << r[j] << "  ";
        }
        cout << "\n";
      }
      else if ( i == 12 )
      {
        cout << "....................\n";
      }
    }

  }

  return;
# undef DIM_MAX
}
//****************************************************************************80
